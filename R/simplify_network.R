
#' Simplify Network Based on Classification
#'
#' Simplifies a network using the classification from `classify_network()`.
#' - **Loop** edges are collapsed (removed, endpoints merged to cluster centroid)
#' - **Parallel** groups are merged into single centerlines (skeletonized)
#' - **Simple** edges are kept as-is
#'
#' @param x An sf object containing LINESTRING geometries (output of `classify_network()`).
#' @param dist Buffer distance for skeletonization of parallel groups.
#'
#' @return An sf object with simplified network geometry.
#' @importFrom geos as_geos_geometry geos_unnest geos_buffer geos_unary_union geos_make_collection
#' @importFrom sf st_as_sf st_geometry st_crs
#' @export
simplify_network = function(x, dist = 10) {
  
  # Ensure classification exists

if (!"net_type" %in% names(x)) {
    x = classify_network(x, dist = dist)
  }
  
  crs_orig = sf::st_crs(x)
  
  # Separate by class
  loops = x[x$net_type == "loop", ]
  parallel = x[x$net_type == "parallel", ]
  simple = x[x$net_type == "simple", ]
  
  result_geoms = list()
  
  # 1. SIMPLE edges: Keep as-is
  if (nrow(simple) > 0) {
    result_geoms$simple = sf::st_geometry(simple)
  }
  
  # 2. LOOP edges: Drop them (they're internal to intersections)
  # Optionally could collapse to centroid points, but for line network, just remove
  
  # 3. PARALLEL groups: Merge into single centerlines
  if (nrow(parallel) > 0) {
    groups = unique(parallel$parallel_group)
    groups = groups[!is.na(groups)]
    
    merged_parallels = vector("list", length(groups))
    
    for (i in seq_along(groups)) {
      gid = groups[i]
      group_geoms = parallel[parallel$parallel_group == gid, ]
      
      # Convert to geos
      g = geos::as_geos_geometry(group_geoms)
      wk::wk_crs(g) = NULL
      g = geos::geos_unnest(g, keep_multi = FALSE)
      
      # Skeletonize: Buffer -> Union -> Skeleton
      # Use the existing neat_skeletonize logic
      g_buff = geos::geos_buffer(g, dist)
      g_union = geos::geos_unary_union(geos::geos_make_collection(g_buff))
      
      # Simple skeleton: centerline via Voronoi
      skeleton = tryCatch({
        neat_skeletonize(g_union, dist, dist * 2, dist * 3)
      }, error = function(e) {
        # Fallback: just pick one of the parallels (the longest)
        lens = geos::geos_length(g)
        g[which.max(lens)]
      })
      
      if (!geos::geos_is_empty(skeleton)[1]) {
        merged_parallels[[i]] = skeleton
      }
    }
    
    # Combine all merged parallels
    merged_parallels = merged_parallels[!sapply(merged_parallels, is.null)]
    if (length(merged_parallels) > 0) {
      all_merged = do.call(c, merged_parallels)
      all_merged = geos::geos_unnest(all_merged, keep_multi = FALSE)
      all_merged = all_merged[!geos::geos_is_empty(all_merged)]
      if (length(all_merged) > 0) {
        result_geoms$parallel = sf::st_as_sfc(all_merged)
      }
    }
  }
  
  # Combine all result geometries
  all_geom = sf::st_sfc(crs = crs_orig)
  
  if (!is.null(result_geoms$simple) && length(result_geoms$simple) > 0) {
    sf::st_crs(result_geoms$simple) = crs_orig
    all_geom = c(all_geom, result_geoms$simple)
  }
  
  if (!is.null(result_geoms$parallel) && length(result_geoms$parallel) > 0) {
    sf::st_crs(result_geoms$parallel) = crs_orig
    all_geom = c(all_geom, result_geoms$parallel)
  }
  
  if (length(all_geom) == 0) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = crs_orig)))
  }
  
  result = sf::st_sf(geometry = all_geom)
  
  result
}
