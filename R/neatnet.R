#' @import geos
#' @import sf
#' @import wk
NULL

#' Buffer and Union Lines
#'
#' Buffers a set of lines and unions them into a single polygon.
#'
#' @param geometry A geos_geometry vector of lines.
#' @param dist Buffer distance.
#' @return A geos_geometry vector (polygons).
#' @export
neat_geometry_buffer <- function(geometry, dist) {
  # Buffer each line
  buffered <- geos_buffer(geometry, dist, params = geos_buffer_params(quad_segs = 16))
  
  # Union all buffers
  geos_unary_union(geos_make_collection(buffered))
}

#' Iteratively Prune Short Dangles
#'
#' Iteratively removes lines that are shorter than a threshold and are dangles.
#' A line is considered a dangle if at least one of its endpoints has a degree of 1
#' (i.e., it is connected to only one other segment) within the current network.
#'
#' @param lines A geos_geometry vector of lines.
#' @param min_length Threshold length for dangles.
#' @return A geos_geometry vector of filtered lines.
#' @noRd
prune_short_dangles_iterative <- function(lines, min_length) {
  if (length(lines) == 0) return(geos_empty())
  
  repeat {
    initial_line_count <- length(lines)
    
    # 1. Calculate lengths
    lens <- geos_length(lines)
    
    # 2. Build topology to find node degrees
    # We need all endpoints of ALL lines to determine connectivity
    starts <- geos_point_start(lines)
    ends <- geos_point_end(lines)
    
    # Convert to comparable representation (Hex is safer for binary equality)
    starts_hex <- geos_write_hex(starts)
    ends_hex <- geos_write_hex(ends)
    
    all_points_hex <- c(starts_hex, ends_hex)
    
    # Count occurrences of points
    # A point with count 1 is a terminal (degree 1)
    point_counts <- table(all_points_hex)
    
    # 3. Determine if each line is a short dangle
    # Look up degrees for start and end of each line
    degree_start <- point_counts[starts_hex]
    degree_end <- point_counts[ends_hex]
    
    # Handle cases where point_counts doesn't contain a key (e.g., if a line becomes isolated after previous removals)
    # This can happen if all other lines connected to a point were removed in previous iteration
    degree_start[is.na(degree_start)] <- 0
    degree_end[is.na(degree_end)] <- 0
    
    # A line is a dangle if EITHER end has degree 1 AND it is short
    is_dangle_candidate <- (degree_start == 1 | degree_end == 1) & (lens < min_length)
    
    # If no dangles to remove, break loop
    if (!any(is_dangle_candidate)) break
    
    # Remove dangles
    lines <- lines[!is_dangle_candidate]
    
    # If no lines left, break
    if (length(lines) == 0) break
    
    # If no change in count, break (shouldn't happen if is_dangle_candidate is true and was not empty)
    # This also acts as a safety against infinite loops if the condition for `any(is_dangle_candidate)` remains true but `lines` is not actually changing due to other issues.
    if (length(lines) == initial_line_count) break 
  }
  
  lines
}

#' Iteratively Prune Short Branches (Preserve Terminals)
#
#' Removes short side-branches (dangles attached to higher-degree junctions)
#' without eroding the main backbone (terminal chains).
#
#' @param lines A geos_geometry vector of lines.
#' @param min_length Threshold length for branch removal.
#' @return A geos_geometry vector of filtered lines.
#' @noRd
prune_short_branches_iterative <- function(lines, min_length) {
  if (length(lines) == 0) return(geos_empty())

  repeat {
    initial_line_count <- length(lines)

    lens <- geos_length(lines)

    starts <- geos_point_start(lines)
    ends <- geos_point_end(lines)

    starts_hex <- geos_write_hex(starts)
    ends_hex <- geos_write_hex(ends)
    all_points_hex <- c(starts_hex, ends_hex)
    point_counts <- table(all_points_hex)

    degree_start <- point_counts[starts_hex]
    degree_end <- point_counts[ends_hex]
    degree_start[is.na(degree_start)] <- 0
    degree_end[is.na(degree_end)] <- 0

    is_branch_tip <- (degree_start == 1 & degree_end >= 3) | (degree_end == 1 & degree_start >= 3)
    is_dangle_candidate <- is_branch_tip & (lens < min_length)

    if (!any(is_dangle_candidate)) break
    lines <- lines[!is_dangle_candidate]

    if (length(lines) == 0) break
    if (length(lines) == initial_line_count) break
  }

  lines
}

#' Skeletonize Polygon
#'
#' Creates a skeleton of a polygon using Voronoi diagrams.
#'
#' @param geometry A geos_geometry vector (polygons).
#' @param dist The buffer distance used for the input geometry (needed for final pruning).
#' @param max_segment_length Distance for densifying the boundary to ensure good Voronoi generation.
#' @param final_min_length Minimum length for final segment filtering.
#' @return A geos_geometry vector (lines).
#' @export
neat_skeletonize <- function(geometry, dist, max_segment_length, final_min_length, simplify = FALSE) {
  # 1. Get Boundary
  boundary <- geos_boundary(geometry)
  
  # 2. Densify Boundary
  # geos_densify adds vertices so that no segment is longer than tolerance
  boundary_densified <- geos_densify(boundary, tolerance = max_segment_length)
  
  # 3. Voronoi Diagram (Edges)
  # Uses nodes of the input geometry
  voronoi_edges <- geos_voronoi_edges(boundary_densified)
  
  # 4. Clip Voronoi to original geometry
  skeleton_raw <- geos_intersection(voronoi_edges, geometry)
  
  # 5. Filter spokes (edges touching the boundary)
  # Unnest to individual segments
  # IMPORTANT: keep_multi = FALSE to ensure we get LINESTRINGS not MULTILINESTRINGS
  edges <- geos_unnest(skeleton_raw, keep_multi = FALSE)
  
  if (length(edges) == 0) return(geos_empty())
  
  # Filter boundary spokes using *midpoint distance to boundary*.
  #
  # Voronoi skeletonization produces lots of boundary-directed "spokes".
  # Using endpoint- or whole-segment distance is too aggressive: valid centreline
  # segments may touch the boundary at their ends (e.g., at polygon ends) but are
  # still part of the spine.
  midpoints <- geos_interpolate(edges, geos_length(edges) / 2)
  dist_mid_to_boundary <- geos_distance(midpoints, boundary)
  spoke_midpoint_tol <- dist * 0.75
  is_spoke <- dist_mid_to_boundary < spoke_midpoint_tol
  
  # Keep only edges that do NOT touch the boundary
  spine_edges <- edges[!is_spoke]
  
  if (length(spine_edges) == 0) {
    # If everything is a spoke (e.g. narrow polygon), return empty or center
    return(geos_empty())
  }
  
  # 6. Pre-process spine: prune very short segments (numerical noise) from initial Voronoi segments
  # This initial pruning helps geos_line_merge to create cleaner longer lines
  # Use a much smaller threshold here, e.g., 0.1, to remove numerical noise without breaking main structure.
  spine_edges_pruned_initial <- prune_short_dangles_iterative(spine_edges, min_length = 0.1)
  
  if (length(spine_edges_pruned_initial) == 0) {
    return(geos_empty())
  }
  
  # Snap endpoints to ensure connectivity before merging.
  # Densify first so endpoints can snap to mid-segment vertices (not just existing
  # endpoints), which helps close small gaps at junctions.
  spine_edges_dense <- geos_densify(spine_edges_pruned_initial, tolerance = dist)
  snapped_edges <- geos_snap(spine_edges_dense, spine_edges_dense, tolerance = dist)
  
  # 7. Merge and Simplify
  # Collect all edges into one collection, then union (to node them), then merge
  collection <- geos_make_collection(snapped_edges) # Use snapped edges here
  spine_noded <- geos_unary_union_prec(collection, grid_size = dist / 5)
  spine_merged <- geos_line_merge(spine_noded)
  
  if (isTRUE(simplify)) {
    # Simplify to smooth artifacts.
    # NOTE: Too much simplification can break junction connectivity.
    simplified <- geos_simplify(spine_merged, tolerance = max_segment_length / 10)
    lines <- geos_unnest(simplified, keep_multi = FALSE)
  } else {
    lines <- geos_unnest(spine_merged, keep_multi = FALSE)
  }

  if (length(lines) == 0) {
    return(geos_empty())
  }

  # Attempt an additional merge pass without forcing re-noding (which can fragment).
  lines <- geos_unnest(geos_merge_lines(geos_make_collection(lines)), keep_multi = FALSE)

  # 8. Final Healing
  # Heal small gaps using densify + snap + union/merge, then prune short dangles.
  heal_densify_tol <- max(max_segment_length, dist * 2)
  lines_dense <- geos_densify(lines, tolerance = heal_densify_tol)
  snapped_final <- geos_snap(lines_dense, lines_dense, tolerance = dist)
  final_noded <- geos_unary_union_prec(geos_make_collection(snapped_final), grid_size = dist / 5)
  final_merged <- geos_merge_lines(final_noded)
  result_lines_merged <- geos_unnest(final_merged, keep_multi = FALSE)

  # Remove short side-branches created by the merge/heal step, while preserving
  # backbone terminals (important for gap-free results).
  result_lines <- prune_short_branches_iterative(result_lines_merged, min_length = final_min_length)

  # If pruning removes everything (can happen if the merged result is still highly
  # segmented), fall back to returning the merged result.
  if (length(result_lines) == 0) {
    result_lines <- result_lines_merged
  }
  
  if (length(result_lines) == 0) {
    message("DEBUG: result_lines empty after final filter.")
    return(geos_empty())
  }

  # Final merge pass after branch pruning.
  result_lines <- geos_unnest(geos_merge_lines(geos_make_collection(result_lines)), keep_multi = FALSE)

  if (isTRUE(simplify)) {
    # 9. Final smoothing (reduce waviness)
    # Simplify geometries, then heal to ensure we don't introduce gaps.
    # The tolerance is tied to `dist` so it scales with the network buffer width.
    smooth_tol <- dist / 2
    result_lines_smooth <- geos_simplify(result_lines, tolerance = smooth_tol)
    result_lines_smooth <- geos_unnest(result_lines_smooth, keep_multi = FALSE)

    if (length(result_lines_smooth) > 0) {
      # Drop degenerate geometries that can trigger GEOS densify errors
      ok <- (geos_num_coordinates(result_lines_smooth) > 1) & (geos_length(result_lines_smooth) > 0)
      ok[is.na(ok)] <- FALSE
      result_lines_smooth <- result_lines_smooth[ok]

      if (length(result_lines_smooth) == 0) {
        return(result_lines)
      }

      result_lines_smooth <- geos_densify(result_lines_smooth, tolerance = dist / 2)
      result_lines_smooth <- geos_snap(result_lines_smooth, result_lines_smooth, tolerance = dist)
      result_lines_smooth <- geos_unary_union_prec(
        geos_make_collection(result_lines_smooth),
        grid_size = dist / 5
      )
      result_lines_smooth <- geos_merge_lines(result_lines_smooth)
      result_lines_smooth <- geos_unnest(result_lines_smooth, keep_multi = FALSE)

      # One more pass to remove short branches created by smoothing.
      result_lines_smooth <- prune_short_branches_iterative(result_lines_smooth, min_length = final_min_length)
      if (length(result_lines_smooth) > 0) {
        result_lines_smooth <- geos_unnest(
          geos_merge_lines(geos_make_collection(result_lines_smooth)),
          keep_multi = FALSE
        )
        result_lines <- result_lines_smooth
      }
    }
  }
  
  result_lines
}

#' Neat Network Simplification
#'
#' Simplifies a road network using a buffer-based skeletonization approach.
#'
#' @param x An sf object representing the network.
#' @param dist Buffer distance (half width of road).
#' @param max_segment_factor Multiplier for dist to determine max_segment_length (default 2).
#' @param final_min_factor Multiplier for dist to determine final_min_length (default 3).
#' @param simplify Logical; if TRUE, applies line simplification and a final smoothing/healing
#'   pass to reduce waviness. Default FALSE.
#' @return An sf object of the simplified network.
#' @export
neatnet <- function(x, dist = 10, max_segment_factor = 2, final_min_factor = 3, simplify = FALSE) {
  # Calculate internal parameters
  max_segment_length_internal <- dist * max_segment_factor
  final_min_length_internal <- dist * final_min_factor
  
  # Convert to geos
  g <- as_geos_geometry(x)
  
  # Step 1: Buffer and Union
  g_buff <- neat_geometry_buffer(g, dist)
  
  # Step 2: Skeletonize
  g_skel <- neat_skeletonize(
    g_buff,
    dist,
    max_segment_length_internal,
    final_min_length_internal,
    simplify = simplify
  )
  
  # Step 3: Convert back to sf
  res <- sf::st_as_sf(g_skel)
  
  # Restore CRS if lost
  if (!is.na(sf::st_crs(x))) {
    sf::st_crs(res) <- sf::st_crs(x)
  }
  
  res
}