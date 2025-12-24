
#' Simplify Network Using sfnetworks (Intersection-First Approach)
#'
#' Uses sfnetworks and DBSCAN clustering to identify and contract
#' complex intersection nodes, then merges parallel edges to centerlines.
#'
#' @param x An sf object containing LINESTRING geometries.
#' @param eps DBSCAN epsilon parameter - nodes within this distance are clustered.
#' @param min_pts DBSCAN minPts parameter (default 1 = assign all nodes to clusters).
#' @param merge_parallels Logical. If TRUE, merges parallel edges to centerlines.
#' @param remove_dangles Logical. If TRUE, removes dangling edges.
#' @param dangle_length Maximum length of dangling edges to remove (set to Inf to remove all).
#'
#' @return An sf object with simplified network geometry.
#' @importFrom sfnetworks as_sfnetwork activate
#' @importFrom tidygraph convert group_components
#' @importFrom dplyr mutate filter
#' @importFrom dbscan dbscan
#' @importFrom sf st_as_sf st_geometry st_coordinates st_crs st_length
#' @importFrom geos as_geos_geometry geos_centroid geos_make_collection
#' @export
simplify_network_sfn = function(x, eps = 20, min_pts = 1, 
                                 merge_parallels = TRUE,
                                 remove_dangles = TRUE,
                                 dangle_length = 50) {
  
  # Check required packages
  if (!requireNamespace("sfnetworks", quietly = TRUE)) {
    stop("Package 'sfnetworks' is required. Install with: install.packages('sfnetworks')")
  }
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required. Install with: install.packages('dbscan')")
  }
  if (!requireNamespace("tidygraph", quietly = TRUE)) {
    stop("Package 'tidygraph' is required. Install with: install.packages('tidygraph')")
  }
  
  crs_orig = sf::st_crs(x)
  
  # 1. Convert to sfnetwork
  net = sfnetworks::as_sfnetwork(x, directed = FALSE)
  
  # 2. Subdivide edges at intersections to ensure proper network topology
  net = tidygraph::convert(net, sfnetworks::to_spatial_subdivision)
  
  # 3. Get node coordinates for DBSCAN clustering
  nodes_sf = sf::st_as_sf(sfnetworks::activate(net, "nodes"))
  node_coords = sf::st_coordinates(nodes_sf)
  
  # 4. Cluster nodes with DBSCAN
  clusters = dbscan::dbscan(node_coords, eps = eps, minPts = min_pts)$cluster
  
  # 5. Add cluster and component info to network
  net = net |>
    sfnetworks::activate("nodes") |>
    dplyr::mutate(
      cls = clusters,
      cmp = tidygraph::group_components()
    )
  
  # 6. Contract nodes in same cluster+component
  # simplify = FALSE to keep parallels for centerline computation
  contracted = tidygraph::convert(
    net,
    sfnetworks::to_spatial_contracted,
    cls, cmp,
    simplify = FALSE  # Keep parallels for now
  )
  
  # 7. Remove loops (self-edges created by contraction)
  edge_data = contracted |> 
    sfnetworks::activate("edges") |> 
    tibble::as_tibble()
  
  is_loop = edge_data$from == edge_data$to
  if (any(is_loop)) {
    contracted = contracted |>
      sfnetworks::activate("edges") |>
      dplyr::slice(which(!is_loop))
  }
  
  # 8. Merge parallel edges to centerlines
  if (merge_parallels) {
    edge_sf = contracted |> sfnetworks::activate("edges") |> sf::st_as_sf()
    edge_data = contracted |> sfnetworks::activate("edges") |> tibble::as_tibble()
    
    # Assign pair_id (order-independent for undirected)
    edge_data$pair_id = paste(
      pmin(edge_data$from, edge_data$to),
      pmax(edge_data$from, edge_data$to),
      sep = "-"
    )
    edge_data$edge_idx = seq_len(nrow(edge_data))
    
    # Get node geometries for snapping
    node_sf = contracted |> sfnetworks::activate("nodes") |> sf::st_as_sf()
    node_points = sf::st_geometry(node_sf)
    
    # Process each unique pair
    unique_pairs = unique(edge_data$pair_id)
    new_edges = vector("list", length(unique_pairs))
    new_from = integer(length(unique_pairs))
    new_to = integer(length(unique_pairs))
    
    for (i in seq_along(unique_pairs)) {
      pid = unique_pairs[i]
      idx = which(edge_data$pair_id == pid)
      
      if (length(idx) == 1) {
        # Single edge - keep as is
        geom = sf::st_geometry(edge_sf[idx, ])
        sf::st_crs(geom) = NA_crs_
        new_edges[[i]] = geom
        new_from[i] = edge_data$from[idx]
        new_to[i] = edge_data$to[idx]
      } else {
        # Multiple parallel edges - compute centerline
        parallel_geoms = sf::st_geometry(edge_sf[idx, ])
        
        # Get start and end nodes
        from_node = edge_data$from[idx[1]]
        to_node = edge_data$to[idx[1]]
        
        # Compute centerline by averaging coordinates
        # Densify lines first, then average
        centerline = compute_centerline(parallel_geoms, 
                                         node_points[from_node], 
                                         node_points[to_node])
        
        new_edges[[i]] = centerline
        new_from[i] = from_node
        new_to[i] = to_node
      }
    }
    
    # Create new edge sf - explicitly strip and reset CRS
    all_geoms = lapply(new_edges, function(g) {
      sf::st_crs(g) = NA_crs_
      g
    })
    new_edge_geom = do.call(c, all_geoms)
    sf::st_crs(new_edge_geom) = crs_orig
    
    # Rebuild network with merged edges
    new_edge_sf = sf::st_sf(
      from = new_from,
      to = new_to,
      geometry = new_edge_geom
    )
    
    # Ensure node_sf has correct CRS
    sf::st_crs(node_sf) = crs_orig
    
    # Create sfnetwork from nodes and edges
    contracted = sfnetworks::sfnetwork(
      nodes = node_sf,
      edges = new_edge_sf,
      directed = FALSE,
      edges_as_lines = TRUE
    )
  }
  
  # 9. Smooth pseudo-nodes (degree-2 nodes)
  smoothed = tidygraph::convert(contracted, sfnetworks::to_spatial_smooth)
  
  # 10. Remove dangling edges
  if (remove_dangles) {
    for (iter in 1:5) {  # Iterate to remove chains
      smoothed = smoothed |>
        sfnetworks::activate("nodes") |>
        dplyr::mutate(degree = tidygraph::centrality_degree())
      
      edge_sf = smoothed |> sfnetworks::activate("edges") |> sf::st_as_sf()
      edge_lens = as.numeric(sf::st_length(edge_sf))
      
      node_degrees = smoothed |> 
        sfnetworks::activate("nodes") |>
        tibble::as_tibble() |>
        dplyr::pull(degree)
      
      edge_data = smoothed |>
        sfnetworks::activate("edges") |>
        tibble::as_tibble()
      
      # Dangle = one endpoint has degree 1
      is_dangle = (node_degrees[edge_data$from] == 1 | node_degrees[edge_data$to] == 1)
      
      # Only remove short dangles OR dangles that are parallel to other edges
      is_short = edge_lens < dangle_length
      
      # Check if dangle is parallel to nearby edges (redundant)
      is_redundant = rep(FALSE, nrow(edge_data))
      if (any(is_dangle)) {
        dangle_idx = which(is_dangle)
        non_dangle_idx = setdiff(seq_len(nrow(edge_data)), dangle_idx)
        
        if (length(non_dangle_idx) > 0) {
          non_dangle_geoms = sf::st_geometry(edge_sf[non_dangle_idx, ])
          for (d in dangle_idx) {
            dangle_geom = sf::st_geometry(edge_sf[d, ])
            # Check if dangle is within eps of any non-dangle edge
            dists = as.numeric(sf::st_distance(dangle_geom, non_dangle_geoms))
            if (min(dists) < eps) {
              is_redundant[d] = TRUE
            }
          }
        }
      }
      
      remove_edge = is_dangle & (is_short | is_redundant)
      
      if (!any(remove_edge)) break
      
      smoothed = smoothed |>
        sfnetworks::activate("edges") |>
        dplyr::slice(which(!remove_edge))
      
      # Re-smooth
      smoothed = tryCatch(
        tidygraph::convert(smoothed, sfnetworks::to_spatial_smooth),
        error = function(e) smoothed
      )
    }
  }
  
  # 11. Extract edges as sf
  edges_sf = sf::st_as_sf(sfnetworks::activate(smoothed, "edges"))
  sf::st_crs(edges_sf) = crs_orig
  
  edges_sf
}


#' Compute Centerline of Parallel Lines
#'
#' @param geoms sfc of parallel linestrings
#' @param start_pt Start point (node geometry)
#' @param end_pt End point (node geometry)
#' @return sfc with single centerline
#' @keywords internal
compute_centerline = function(geoms, start_pt, end_pt) {
  # Convert to geos for efficient processing
  g = geos::as_geos_geometry(geoms)
  
  # Densify lines to have similar vertex counts
  # Use small segment length
  g_dense = geos::geos_densify(g, tolerance = 5)
  
  # Extract all coordinates
  all_coords = lapply(seq_along(g_dense), function(i) {
    coords = wk::wk_coords(g_dense[i])
    data.frame(x = coords$x, y = coords$y)
  })
  
  # Normalize to same number of points along the curve
  n_pts = 50
  interp_coords = lapply(all_coords, function(df) {
    if (nrow(df) < 2) return(df)
    # Cumulative distance
    dx = diff(df$x)
    dy = diff(df$y)
    dists = c(0, cumsum(sqrt(dx^2 + dy^2)))
    total_len = max(dists)
    if (total_len == 0) return(df[1, ])
    # Interpolate at regular intervals
    target_dists = seq(0, total_len, length.out = n_pts)
    x_interp = approx(dists, df$x, xout = target_dists)$y
    y_interp = approx(dists, df$y, xout = target_dists)$y
    data.frame(x = x_interp, y = y_interp)
  })
  
  # Average coordinates
  n_lines = length(interp_coords)
  avg_x = rowMeans(do.call(cbind, lapply(interp_coords, `[[`, "x")))
  avg_y = rowMeans(do.call(cbind, lapply(interp_coords, `[[`, "y")))
  
  # Force start and end to match node points
  start_coord = sf::st_coordinates(start_pt)
  end_coord = sf::st_coordinates(end_pt)
  avg_x[1] = start_coord[1]
  avg_y[1] = start_coord[2]
  avg_x[n_pts] = end_coord[1]
  avg_y[n_pts] = end_coord[2]
  
  # Create linestring
  coords_mat = cbind(avg_x, avg_y)
  centerline = sf::st_sfc(sf::st_linestring(coords_mat))
  
  centerline
}
