
#' Simplify Network Using sfnetworks (Intersection-First Approach)
#'
#' Uses sfnetworks and DBSCAN clustering to identify and contract
#' complex intersection nodes, then handles parallel edges.
#'
#' @param x An sf object containing LINESTRING geometries.
#' @param eps DBSCAN epsilon parameter - nodes within this distance are clustered.
#' @param min_pts DBSCAN minPts parameter (default 1 = assign all nodes to clusters).
#'
#' @return An sf object with simplified network geometry.
#' @importFrom sfnetworks as_sfnetwork activate st_network_paths
#' @importFrom tidygraph convert group_components
#' @importFrom dplyr mutate
#' @importFrom dbscan dbscan
#' @importFrom sf st_as_sf st_geometry st_coordinates st_crs
#' @export
simplify_network_sfn = function(x, eps = 10, min_pts = 1) {
  
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
  contracted = tidygraph::convert(
    net,
    sfnetworks::to_spatial_contracted,
    cls, cmp,
    simplify = TRUE
  )
  
  # 7. Smooth pseudo-nodes (degree-2 nodes that just connect two edges)
  smoothed = tidygraph::convert(contracted, sfnetworks::to_spatial_smooth)
  
  # 8. Extract edges as sf
  edges_sf = sf::st_as_sf(sfnetworks::activate(smoothed, "edges"))
  
  # Restore CRS
  sf::st_crs(edges_sf) = crs_orig
  
  edges_sf
}
