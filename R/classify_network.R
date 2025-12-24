
#' Classify Network Features (Graph Theory Approach)
#'
#' Classifies network edges based on a graph-theoretical model where 'complex intersections'
#' are clustered nodes, and 'parallel' edges are multiple edges connecting the same two clusters.
#'
#' @param x An sf object containing LINESTRING geometries.
#' @param dist Buffer distance (in map units) to define the size of a 'node cluster'.
#'   Nodes within `dist` of each other will be merged into a single intersection centroid.
#' @param frechet_threshold Threshold for Frechet distance to confirm parallelism.
#'   Parallel candidates (edges connecting the same node clusters) must be within this
#'   distance of the shortest edge in the group to be confirmed as parallel.
#'
#' @return An sf object with added columns:
#'   \itemize{
#'     \item \code{net_type}: Character. One of "loop" (internal to intersection), "parallel" (part of a parallel group), or "simple" (single link).
#'     \item \code{cluster_u}: ID of the start node cluster.
#'     \item \code{cluster_v}: ID of the end node cluster.
#'     \item \code{parallel_group}: ID for the parallel group (if parallel).
#'   }
#' @importFrom geos as_geos_geometry geos_unnest geos_point_start geos_point_end geos_buffer geos_unary_union geos_make_collection geos_idx_strtree geos_strtree_query geos_length geos_distance_frechet
#' @importFrom sf st_as_sf st_geometry
#' @importFrom dplyr mutate filter group_by ungroup rowwise n
#' @export
classify_network = function(x, dist = 10, frechet_threshold = 29) {
  
  # Ensure input is suitable
  # We assume projected CRS for distance operations
  
  # 1. Geometry Prep (GEOS)
  g <- geos::as_geos_geometry(x)
  # Strip CRS from geos object to avoid mismatches during internal ops (sf/geos interaction)
  # Note: This affects only the geos object 'g'
  wk::wk_crs(g) <- NULL
  
  g_lines <- geos::geos_unnest(g, keep_multi = FALSE)
  
  if (length(g_lines) != nrow(x)) {
    warning("Multi-linestrings were unnested. Row counts may mismatch if input had MULTILINESTRINGs.")
    # For prototype, assume simple 1-1 mapping or robustify
  }
  
  # 2. Node Clustering
  starts <- geos::geos_point_start(g_lines)
  ends <- geos::geos_point_end(g_lines)
  all_points <- c(starts, ends)
  
  # Buffer points by dist/2 => Union => Connected Components
  node_buffs <- geos::geos_buffer(all_points, dist / 2)
  node_clusters_geom <- geos::geos_unary_union(geos::geos_make_collection(node_buffs))
  node_clusters_poly <- geos::geos_unnest(node_clusters_geom, keep_multi = FALSE)
  
  # Build Index for Node -> Cluster mapping
  cluster_tree <- geos::geos_strtree(node_clusters_poly)
  
  get_cluster_id <- function(pts) {
    q <- geos::geos_strtree_query(cluster_tree, pts)
    # Take first hitting polygon
    vapply(q, function(idx) if(length(idx) > 0) as.integer(idx[1]) else NA_integer_, integer(1))
  }
  
  # Map endpoints to clusters
  start_ids <- get_cluster_id(starts)
  end_ids <- get_cluster_id(ends)
  
  # 3. Graph Analysis
  # Create a lightweight dataframe for logic
  # We use base R / dplyr for speed
  df <- data.frame(
    edge_id = seq_along(g_lines),
    u = start_ids,
    v = end_ids
  )
  
  # Identify Loops (Internal to cluster)
  df$is_loop <- df$u == df$v
  
  # Identify Parallel Candidates
  # Sort u,v to handle direction
  pmin_uv <- pmin(df$u, df$v)
  pmax_uv <- pmax(df$u, df$v)
  df$pair_id <- paste(pmin_uv, pmax_uv, sep = "-")
  
  # Group info
  # We can't use complex dplyr pipes easily without importing everything, but basic is fine
  # Or use base R for 'ave'
  df$group_count <- ave(df$edge_id, df$pair_id, FUN = length)
  
  # Initialize classifications
  # types: "simple", "loop", "parallel"
  result_types <- rep("simple", nrow(df))
  result_types[df$is_loop] <- "loop"
  
  # Candidates for parallel: count > 1 AND not loop
  is_candidate <- (df$group_count > 1) & (!df$is_loop)
  
  parallel_group_id <- rep(NA_integer_, nrow(df))
  
  if (any(is_candidate)) {
    # Verify Frechet
    # Get unique pairs that are candidates
    candidate_pairs <- unique(df$pair_id[is_candidate])
    
    # We assign a unique group ID for verified parallel sets
    # Start IDs after max cluster ID to avoid confusion? Or just 1..N
    pg_counter <- 1
    
    for (pid in candidate_pairs) {
      idx <- which(df$pair_id == pid & !df$is_loop)
      if (length(idx) < 2) next # Should be covered by is_candidate logic
      
      geoms <- g_lines[idx]
      lens <- geos::geos_length(geoms)
      ref_idx <- which.min(lens)
      
      # Frechet check
      # densify = 0.5 might be slow for long lines. 
      # Adaptive densify? Or just rely on segment vertices?
      # Using NULL (default) relies on vertices.
      # User linked frechet docs.
      # Let's use 0.5 as in prototype if it wasn't too slow.
      # Prototype took < 1s for 1000 features. Safe.
      fd <- geos::geos_distance_frechet(geoms, geoms[ref_idx], densify = 0.5)
      
      ok <- fd <= frechet_threshold
      
      # If we have > 1 confirmed parallel
      if (sum(ok) > 1) {
        # Mark these as parallel
        result_types[idx[ok]] = "parallel"
        parallel_group_id[idx[ok]] = pg_counter
        pg_counter = pg_counter + 1
        
        # Note: those that failed check revert to 'simple' (default)
        # effectively splitting the group
      }
    }
  }
  
  # ============================================================
  # 4. SPATIAL PARALLEL DETECTION (Second Pass)
  # This catches parallel edges that don't share node clusters
  # (e.g., long dual carriageways with different junction patterns)
  # ============================================================
  
  # Only check "simple" edges that weren't already classified
  simple_idx = which(result_types == "simple")
  
  if (length(simple_idx) > 1) {
    # Build spatial index on simple edges
    simple_geoms = g_lines[simple_idx]
    simple_tree = geos::geos_strtree(simple_geoms)
    
    # Compute bearings for direction similarity
    get_bearing = function(geom) {
      p0 = geos::geos_point_start(geom)
      p1 = geos::geos_point_end(geom)
      dx = geos::geos_x(p1) - geos::geos_x(p0)
      dy = geos::geos_y(p1) - geos::geos_y(p0)
      ang = atan2(dy, dx) * 180 / pi
      ang %% 180  # Direction-invariant (0-180)
    }
    
    simple_bearings = vapply(seq_along(simple_geoms), function(i) get_bearing(simple_geoms[i]), numeric(1))
    
    # Query neighbors within dist
    simple_bboxes = geos::geos_envelope(simple_geoms)
    expanded = geos::geos_create_rectangle(
      geos::geos_xmin(simple_bboxes) - dist,
      geos::geos_ymin(simple_bboxes) - dist,
      geos::geos_xmax(simple_bboxes) + dist,
      geos::geos_ymax(simple_bboxes) + dist
    )
    candidates = geos::geos_strtree_query(simple_tree, expanded)
    
    # Track spatial parallel groups
    spatial_parallel_adj = vector("list", length(simple_geoms))
    bearing_tol = 20  # degrees
    
    for (i in seq_along(simple_geoms)) {
      nbrs = candidates[[i]]
      nbrs = nbrs[nbrs > i]  # Only forward to avoid duplicates
      if (length(nbrs) == 0) next
      
      # Filter by bearing similarity
      bi = simple_bearings[i]
      bj = simple_bearings[nbrs]
      bearing_diff = abs(((bj - bi + 90) %% 180) - 90)
      similar_bearing = bearing_diff <= bearing_tol
      
      if (!any(similar_bearing)) next
      nbrs = nbrs[similar_bearing]
      
      # Filter: must not intersect (parallel lines don't cross)
      does_intersect = geos::geos_intersects(simple_geoms[i], simple_geoms[nbrs])
      nbrs = nbrs[!does_intersect]
      
      if (length(nbrs) == 0) next
      
      # Buffer overlap check: if buffers overlap, they run alongside each other
      # This works even when segments have different lengths (e.g., one side has extra junction)
      buff_i = geos::geos_buffer(simple_geoms[i], frechet_threshold / 2)
      buffs_j = geos::geos_buffer(simple_geoms[nbrs], frechet_threshold / 2)
      overlaps = geos::geos_intersects(buff_i, buffs_j)
      
      # Additional check: actual minimum distance should be reasonable (< frechet_threshold)
      min_dists = geos::geos_distance(simple_geoms[i], simple_geoms[nbrs])
      is_close = min_dists <= frechet_threshold
      
      is_par = overlaps & is_close
      
      if (any(is_par)) {
        spatial_parallel_adj[[i]] = c(spatial_parallel_adj[[i]], nbrs[is_par])
        for (j in nbrs[is_par]) {
          spatial_parallel_adj[[j]] = c(spatial_parallel_adj[[j]], i)
        }
      }
    }
    
    # Build connected components from adjacency
    has_spatial_parallel = vapply(spatial_parallel_adj, function(x) length(x) > 0, logical(1))
    
    if (any(has_spatial_parallel)) {
      # Simple union-find for connected components
      parent = seq_along(simple_geoms)
      find_root = function(a) {
        while (parent[a] != a) {
          parent[a] <<- parent[parent[a]]
          a = parent[a]
        }
        a
      }
      unite = function(a, b) {
        ra = find_root(a)
        rb = find_root(b)
        if (ra != rb) parent[rb] <<- ra
      }
      
      for (i in which(has_spatial_parallel)) {
        for (j in spatial_parallel_adj[[i]]) {
          unite(i, j)
        }
      }
      
      roots = vapply(seq_along(simple_geoms), find_root, integer(1))
      unique_roots = unique(roots[has_spatial_parallel])
      
      for (r in unique_roots) {
        members = which(roots == r)
        if (length(members) > 1) {
          orig_idx = simple_idx[members]
          result_types[orig_idx] = "parallel"
          parallel_group_id[orig_idx] = pg_counter
          pg_counter = pg_counter + 1
        }
      }
    }
  }
  
  # Attach results to x
  x$net_type = result_types
  x$cluster_u = df$u
  x$cluster_v = df$v
  x$parallel_group = parallel_group_id
  
  x
}

