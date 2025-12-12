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
#' @return A named list with elements `n_features`, `total_length`, and `n_components`.
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
    return(list(n_features = 0L, total_length = 0, n_components = 0L))
  }

  # Clean up to lines
  g <- geos::geos_unnest(g, keep_multi = FALSE)
  g <- g[!geos::geos_is_empty(g)]

  if (length(g) == 0) {
    return(list(n_features = 0L, total_length = 0, n_components = 0L))
  }

  total_length <- sum(as.numeric(geos::geos_length(g)), na.rm = TRUE)

  # Snap to precision grid via unary union with precision
  if (!is.null(grid_size) && is.finite(grid_size) && grid_size > 0) {
    g <- geos::geos_unary_union_prec(geos::geos_make_collection(g), grid_size = grid_size)
    g <- geos::geos_unnest(g, keep_multi = FALSE)
  }

  if (node) {
    g <- geos::geos_node(geos::geos_make_collection(g))
    g <- geos::geos_unnest(g, keep_multi = FALSE)
  }

  g <- g[!geos::geos_is_empty(g)]
  if (length(g) == 0) {
    return(list(n_features = 0L, total_length = total_length, n_components = 0L))
  }

  # Compute connected components using union-find on endpoint IDs.
  # Using HEX ensures stable point equality after precision snapping.
  starts_hex <- geos::geos_write_hex(geos::geos_point_start(g))
  ends_hex <- geos::geos_write_hex(geos::geos_point_end(g))

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
    n_features = as.integer(length(g)),
    total_length = total_length,
    n_components = as.integer(n_components)
  )
}
