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

#' Skeletonize Polygon
#'
#' Creates a skeleton of a polygon using Voronoi diagrams.
#'
#' @param geometry A geos_geometry vector (polygons).
#' @param max_segment_length Distance for densifying the boundary to ensure good Voronoi generation.
#' @return A geos_geometry vector (lines).
#' @export
neat_skeletonize <- function(geometry, max_segment_length = 1) {
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
  
  # Check intersection with boundary
  # Use distance check for robustness against precision issues
  # Edges starting/ending at boundary will have distance ~ 0
  # We use a liberal tolerance (e.g. 0.1) to catch edges that are slightly detached due to artifacts,
  # while ensuring we don't remove the spine (which is typically 'dist' away, e.g. >1m).
  dist_to_boundary <- geos_distance(edges, boundary)
  is_spoke <- dist_to_boundary < 0.1
  
  # Keep only edges that do NOT touch the boundary
  spine_edges <- edges[!is_spoke]
  
  if (length(spine_edges) == 0) {
    # If everything is a spoke (e.g. narrow polygon), return empty or center
    return(geos_empty())
  }
  
  # 6. Merge and Simplify
  # Collect all edges into one collection, then union (to node them), then merge
  collection <- geos_make_collection(spine_edges)
  spine_noded <- geos_unary_union(collection)
  spine_merged <- geos_line_merge(spine_noded)
  
  # Simplify to smooth artifacts
  simplified <- geos_simplify(spine_merged, tolerance = max_segment_length / 2)
  
  # Unnest to return individual LineStrings
  lines <- geos_unnest(simplified, keep_multi = FALSE)
  
  # Filter short lines (artifacts/dangles)
  # A reasonable threshold is max_segment_length * 2
  lines[geos_length(lines) > (max_segment_length * 2)]
}

#' Neat Network Simplification
#'
#' Simplifies a road network using a buffer-based skeletonization approach.
#'
#' @param x An sf object representing the network.
#' @param dist Buffer distance (half width of road).
#' @param max_segment_length Densification distance.
#' @return An sf object of the simplified network.
#' @export
neatnet <- function(x, dist = 10, max_segment_length = 1) {
  # Convert to geos
  g <- as_geos_geometry(x)
  
  # Step 1: Buffer and Union
  g_buff <- neat_geometry_buffer(g, dist)
  
  # Step 2: Skeletonize
  g_skel <- neat_skeletonize(g_buff, max_segment_length)
  
  # Step 3: Convert back to sf
  res <- sf::st_as_sf(g_skel)
  
  # Restore CRS if lost
  if (!is.na(sf::st_crs(x))) {
    sf::st_crs(res) <- sf::st_crs(x)
  }
  
  res
}