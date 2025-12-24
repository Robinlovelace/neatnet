
#' Simplify Network Using sfnetworks (Intersection-First Approach)
#'
#' Uses sfnetworks and DBSCAN clustering to identify and contract
#' complex intersection nodes, then removes parallel edges and dangles.
#'
#' @param x An sf object containing LINESTRING geometries.
#' @param eps DBSCAN epsilon parameter - nodes within this distance are clustered.
#' @param min_pts DBSCAN minPts parameter (default 1 = assign all nodes to clusters).
#' @param remove_parallels Logical. If TRUE, removes parallel edges (multiple edges
#'   between the same node pair), keeping the shortest.
#' @param remove_dangles Logical. If TRUE, removes dangling edges (dead ends).
#' @param dangle_length Maximum length of dangling edges to remove (in map units).
#'
#' @return An sf object with simplified network geometry.
#' @importFrom sfnetworks as_sfnetwork activate
#' @importFrom tidygraph convert group_components
#' @importFrom dplyr mutate filter
#' @importFrom dbscan dbscan
#' @importFrom sf st_as_sf st_geometry st_coordinates st_crs st_length
#' @export
simplify_network_sfn = function(x, eps = 20, min_pts = 1, 
                                 remove_parallels = TRUE,
                                 remove_dangles = TRUE,
                                 dangle_length = 30) {
  
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
  # simplify = TRUE removes loops and multiple edges created by contraction
  contracted = tidygraph::convert(
    net,
    sfnetworks::to_spatial_contracted,
    cls, cmp,
    simplify = TRUE
  )
  
  # 7. Remove parallel edges (multiple edges between same node pair)
  # After contraction, parallels become multiple edges - keep shortest
  if (remove_parallels) {
    contracted = contracted |>
      sfnetworks::activate("edges") |>
      dplyr::mutate(
        edge_len = sf::st_length(sfnetworks::activate(contracted, "edges") |> sf::st_as_sf())
      )
    
    # Get edge data
    edges_data = contracted |> 
      sfnetworks::activate("edges") |> 
      tibble::as_tibble()
    
    # Identify parallels: same from-to pair (order-independent for undirected)
    edges_data$pair_id = paste(
      pmin(edges_data$from, edges_data$to),
      pmax(edges_data$from, edges_data$to),
      sep = "-"
    )
    
    # Keep only shortest edge per pair
    keep_idx = edges_data |>
      dplyr::mutate(row_id = dplyr::row_number()) |>
      dplyr::group_by(pair_id) |>
      dplyr::slice_min(edge_len, n = 1, with_ties = FALSE) |>
      dplyr::pull(row_id)
    
    # Filter to keep only selected edges
    contracted = contracted |>
      sfnetworks::activate("edges") |>
      dplyr::slice(keep_idx)
  }
  
  # 8. Smooth pseudo-nodes (degree-2 nodes that just connect two edges)
  smoothed = tidygraph::convert(contracted, sfnetworks::to_spatial_smooth)
  
  # 9. Remove dangling edges (short dead-ends)
  if (remove_dangles) {
    # Iterate a few times to handle chains of dangles
    for (i in 1:3) {
      smoothed = smoothed |>
        sfnetworks::activate("nodes") |>
        dplyr::mutate(degree = tidygraph::centrality_degree())
      
      # Get edge lengths
      edge_sf = smoothed |> sfnetworks::activate("edges") |> sf::st_as_sf()
      edge_lens = as.numeric(sf::st_length(edge_sf))
      
      # Get node degrees
      node_degrees = smoothed |> 
        sfnetworks::activate("nodes") |>
        tibble::as_tibble() |>
        dplyr::pull(degree)
      
      # Get from/to for each edge
      edge_data = smoothed |>
        sfnetworks::activate("edges") |>
        tibble::as_tibble()
      
      # Edge is dangle if one endpoint has degree 1 and edge is short
      is_dangle = (node_degrees[edge_data$from] == 1 | node_degrees[edge_data$to] == 1) &
                  edge_lens < dangle_length
      
      if (!any(is_dangle)) break
      
      # Remove dangles
      smoothed = smoothed |>
        sfnetworks::activate("edges") |>
        dplyr::slice(which(!is_dangle))
      
      # Re-smooth after removing dangles
      smoothed = tidygraph::convert(smoothed, sfnetworks::to_spatial_smooth)
    }
  }
  
  # 10. Extract edges as sf
  edges_sf = sf::st_as_sf(sfnetworks::activate(smoothed, "edges"))
  
  # Restore CRS
  sf::st_crs(edges_sf) = crs_orig
  
  edges_sf
}
