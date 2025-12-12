#' Network summary metrics
#'
#' Compute basic network summary metrics for a line network.
#'
#' @param x An sf object (or sfc) containing LINESTRING/MULTILINESTRING geometries.
#' @param grid_size Precision grid size in map units used to snap coordinates before
#'   topology calculations. Set to `NULL` to skip precision snapping.
#' @param node Logical; if `TRUE`, node the linework before computing connected
#'   components.
#'
#' @return A named list with elements `n_features`, `n_vertices`, `total_length`, and `n_components`.
#'   Also returns `n_features_topo` and `n_vertices_topo`, which reflect the
#'   (optional) precision-snapped + noded linework used for topology.
#' @export
net_summary <- function(x, grid_size = 0.1, node = TRUE) {
  if (inherits(x, "sf")) {
    geom <- sf::st_geometry(x)
  } else if (inherits(x, "sfc")) {
    geom <- x
  } else {
    stop("`x` must be an sf or sfc object.")
  }

  g <- geos::as_geos_geometry(geom)
  g <- g[!geos::geos_is_empty(g)]

  if (length(g) == 0) {
    return(list(n_features = 0L, n_vertices = 0L, total_length = 0, n_components = 0L))
  }

  # Clean up to lines
  g <- geos::geos_unnest(g, keep_multi = FALSE)
  g <- g[!geos::geos_is_empty(g)]

  if (length(g) == 0) {
    return(list(n_features = 0L, n_vertices = 0L, total_length = 0, n_components = 0L))
  }

  total_length <- sum(as.numeric(geos::geos_length(g)), na.rm = TRUE)
  n_features <- as.integer(length(g))
  n_vertices <- sum(as.integer(geos::geos_num_coordinates(g)), na.rm = TRUE)

  g_topo <- g

  # Snap to precision grid via unary union with precision
  if (!is.null(grid_size) && is.finite(grid_size) && grid_size > 0) {
    g_topo <- geos::geos_unary_union_prec(geos::geos_make_collection(g_topo), grid_size = grid_size)
    g_topo <- geos::geos_unnest(g_topo, keep_multi = FALSE)
  }

  if (node) {
    g_topo <- geos::geos_node(geos::geos_make_collection(g_topo))
    g_topo <- geos::geos_unnest(g_topo, keep_multi = FALSE)
  }

  g_topo <- g_topo[!geos::geos_is_empty(g_topo)]
  if (length(g_topo) == 0) {
    return(list(
      n_features = n_features,
      n_vertices = as.integer(n_vertices),
      total_length = total_length,
      n_components = 0L,
      n_features_topo = 0L,
      n_vertices_topo = 0L
    ))
  }

  n_features_topo <- as.integer(length(g_topo))
  n_vertices_topo <- sum(as.integer(geos::geos_num_coordinates(g_topo)), na.rm = TRUE)

  # Compute connected components using union-find on endpoint IDs.
  # Using HEX ensures stable point equality after precision snapping.
  starts_hex <- geos::geos_write_hex(geos::geos_point_start(g_topo))
  ends_hex <- geos::geos_write_hex(geos::geos_point_end(g_topo))

  nodes <- unique(c(starts_hex, ends_hex))
  idx <- seq_along(nodes)
  names(idx) <- nodes

  u <- as.integer(idx[starts_hex])
  v <- as.integer(idx[ends_hex])

  parent <- seq_along(nodes)
  find <- function(a) {
    while (parent[a] != a) {
      parent[a] <<- parent[parent[a]]
      a <- parent[a]
    }
    a
  }
  unite <- function(a, b) {
    ra <- find(a)
    rb <- find(b)
    if (ra != rb) parent[rb] <<- ra
  }

  for (i in seq_along(u)) {
    unite(u[[i]], v[[i]])
  }

  reps <- vapply(seq_along(nodes), find, integer(1))
  n_components <- length(unique(reps))

  list(
    n_features = n_features,
    n_vertices = as.integer(n_vertices),
    total_length = total_length,
    n_components = as.integer(n_components),
    n_features_topo = n_features_topo,
    n_vertices_topo = as.integer(n_vertices_topo)
  )
}
