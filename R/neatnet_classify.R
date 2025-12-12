#' Classify line features into groups (prototype)
#'
#' Prototype preprocessing that clusters input lines into groups and classifies
#' each group as `parallel`, `complex`, both, or `keep`.
#'
#' Grouping is based on overlap of buffers of width `dist`.
#' Parallelism is detected by low orientation dispersion within a group.
#' Complexity is detected by the presence of high-degree junctions after noding.
#'
#' @param x An sf object (or sfc) containing LINESTRING/MULTILINESTRING geometries.
#' @param dist Buffer distance (map units). Used both for grouping and for
#'   tolerance values in the topology check.
#' @param grid_size Precision grid size used before topology checks.
#'   Defaults to `dist/5`.
#' @param parallel_tol_deg Maximum angular deviation (degrees) within a group to
#'   be considered parallel. Default 20.
#' @param min_parallel_lines Minimum number of features in a group for it to be
#'   considered parallel. Default 2.
#' @param complex_min_neighbors Minimum number of nearby neighbours (within
#'   `dist`) required for a feature to be considered in a locally complex area.
#'   Default 8.
#' @param complex_min_disp_deg Minimum angular dispersion (degrees, 0..90) among
#'   nearby neighbours required to flag complexity. Default 60.
#'
#' @return An sf object with the original geometries and added columns:
#'   `group_id`, `n_in_group`, `bearing_deg`, `feature_parallel`, `feature_complex`,
#'   `group_parallel`, `group_complex`, and `group_class`.
#' @export
neatnet_classify_groups <- function(
  x,
  dist,
  grid_size = dist / 5,
  parallel_tol_deg = 20,
  min_parallel_lines = 2,
  complex_min_neighbors = 8,
  complex_min_disp_deg = 60
) {
  old_s2 <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
  suppressMessages(sf::sf_use_s2(FALSE))

  if (inherits(x, "sf")) {
    sf_x <- x
  } else if (inherits(x, "sfc")) {
    sf_x <- sf::st_sf(geometry = x)
  } else {
    stop("`x` must be an sf or sfc object.")
  }

  geom <- sf::st_geometry(sf_x)
  if (length(geom) == 0) {
    sf_x$group_id <- integer(0)
    sf_x$n_in_group <- integer(0)
    sf_x$bearing_deg <- numeric(0)
    sf_x$feature_parallel <- logical(0)
    sf_x$feature_complex <- logical(0)
    sf_x$group_parallel <- logical(0)
    sf_x$group_complex <- logical(0)
    sf_x$group_class <- character(0)
    return(sf_x)
  }

  # 1) Neighbourhood graph (within distance)
  # Using within-distance avoids transitive closure from the *road network*
  # connectivity itself (which otherwise merges almost the entire city into one
  # group). We'll build groups based on near-parallel adjacency instead.
  hits <- sf::st_is_within_distance(geom, geom, dist = dist, sparse = TRUE)
  # Precompute intersection/touch adjacency to avoid classifying connected
  # (collinear) street segments as "parallel".
  intersects <- sf::st_intersects(geom, geom, sparse = TRUE)

  # 2) Feature-level bearing (0..180, direction-invariant)
  sf_x$bearing_deg <- vapply(geom, .neatnet_bearing_deg_180, numeric(1))

  # 3) Feature-level parallel flag: has at least one near-parallel neighbour
  # within the neighbourhood.
  feature_parallel <- logical(length(geom))
  parallel_adj <- vector("list", length(geom))
  if (length(geom) >= 2) {
    for (i in seq_along(hits)) {
      nbrs <- hits[[i]]
      # exclude self and any geometries that intersect/touch this feature
      nbrs <- setdiff(nbrs, c(i, intersects[[i]]))
      if (length(nbrs) == 0) next
      bi <- sf_x$bearing_deg[[i]]
      if (!is.finite(bi)) next
      bj <- sf_x$bearing_deg[nbrs]
      bj <- bj[is.finite(bj)]
      if (length(bj) == 0) next
      diffs <- abs(((bj - bi + 90) %% 180) - 90)
      keep <- diffs <= parallel_tol_deg
      if (any(keep, na.rm = TRUE)) {
        feature_parallel[[i]] <- TRUE
        parallel_adj[[i]] <- nbrs[keep]
      }
    }
  }
  sf_x$feature_parallel <- feature_parallel

  # Build parallel-based groups: connected components of the parallel adjacency.
  # Features with no parallel neighbours become singleton groups.
  parallel_adj2 <- lapply(seq_along(geom), function(i) {
    unique(c(i, parallel_adj[[i]]))
  })
  sf_x$group_id <- .neatnet_components_from_adjlist(parallel_adj2)

  # 4) Feature-level complex flag (conservative): lots of nearby neighbours and
  # high angular dispersion.
  # This aims to capture roundabouts / complex junction areas without marking
  # most ordinary street segments as complex.
  feature_complex <- logical(length(geom))
  if (length(geom) >= 2) {
    for (i in seq_along(hits)) {
      nbrs <- hits[[i]]
      nbrs <- nbrs[nbrs != i]
      if (length(nbrs) < complex_min_neighbors) next

      b <- sf_x$bearing_deg[nbrs]
      b <- b[is.finite(b)]
      if (length(b) < complex_min_neighbors) next

      disp <- .neatnet_max_angular_dev_180(b)
      feature_complex[[i]] <- is.finite(disp) && (disp >= complex_min_disp_deg)
    }
  }
  sf_x$feature_complex <- feature_complex

  # 5) Group-level summaries
  group_sizes <- table(sf_x$group_id)
  sf_x$n_in_group <- as.integer(group_sizes[as.character(sf_x$group_id)])

  groups <- sort(unique(sf_x$group_id))
  group_parallel <- logical(length(groups))
  group_complex <- logical(length(groups))
  names(group_parallel) <- as.character(groups)
  names(group_complex) <- as.character(groups)

  for (gid in groups) {
    idx <- which(sf_x$group_id == gid)
    group_parallel[as.character(gid)] <- (length(idx) >= min_parallel_lines) && any(sf_x$feature_parallel[idx])

    group_complex[as.character(gid)] <- any(sf_x$feature_complex[idx])
  }

  sf_x$group_parallel <- unname(group_parallel[as.character(sf_x$group_id)])
  sf_x$group_complex <- unname(group_complex[as.character(sf_x$group_id)])

  # 4) Class label
  sf_x$group_class <- ifelse(
    sf_x$group_parallel & sf_x$group_complex,
    "parallel+complex",
    ifelse(
      sf_x$group_parallel,
      "parallel",
      ifelse(sf_x$group_complex, "complex", "keep")
    )
  )

  sf_x
}

.neatnet_components_from_adjlist <- function(adj) {
  n <- length(adj)
  parent <- seq_len(n)
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

  for (i in seq_len(n)) {
    nbrs <- adj[[i]]
    if (length(nbrs) == 0) next
    for (j in nbrs) unite(i, j)
  }

  reps <- vapply(seq_len(n), find, integer(1))
  # Make consecutive IDs starting at 1
  as.integer(factor(reps, levels = unique(reps)))
}

.neatnet_bearing_deg_180 <- function(geom) {
  # Handles LINESTRING/MULTILINESTRING; uses first and last coordinate.
  # Returns degrees in [0, 180), so opposite directions are considered parallel.
  coords <- sf::st_coordinates(geom)
  if (is.null(coords) || nrow(coords) < 2) return(NA_real_)
  first <- coords[1, 1:2]
  last <- coords[nrow(coords), 1:2]
  dx <- last[1] - first[1]
  dy <- last[2] - first[2]
  if (!is.finite(dx) || !is.finite(dy) || (dx == 0 && dy == 0)) return(NA_real_)
  ang <- atan2(dy, dx) * 180 / pi
  ang <- ang %% 180
  if (ang < 0) ang <- ang + 180
  ang
}

.neatnet_max_angular_dev_180 <- function(bearings_deg) {
  # Compute max deviation from the circular mean on a 180-degree circle.
  # Use doubled angles trick to account for 180 periodicity.
  if (length(bearings_deg) == 0) return(NA_real_)
  th <- bearings_deg * pi / 180
  th2 <- 2 * th
  cbar <- mean(cos(th2))
  sbar <- mean(sin(th2))
  if (!is.finite(cbar) || !is.finite(sbar)) return(NA_real_)
  mean_th2 <- atan2(sbar, cbar)
  mean_th <- (mean_th2 / 2)
  dev <- abs(((bearings_deg - (mean_th * 180 / pi) + 90) %% 180) - 90)
  max(dev, na.rm = TRUE)
}

.neatnet_junction_points <- function(geom, grid_size) {
  g <- geos::as_geos_geometry(geom)
  g <- g[!geos::geos_is_empty(g)]
  if (length(g) == 0) return(sf::st_sfc(crs = sf::st_crs(geom)))

  g <- geos::geos_unnest(g, keep_multi = FALSE)
  g <- g[!geos::geos_is_empty(g)]
  if (length(g) == 0) return(sf::st_sfc(crs = sf::st_crs(geom)))

  if (!is.null(grid_size) && is.finite(grid_size) && grid_size > 0) {
    g <- geos::geos_unary_union_prec(geos::geos_make_collection(g), grid_size = grid_size)
    g <- geos::geos_unnest(g, keep_multi = FALSE)
  }

  g <- geos::geos_node(geos::geos_make_collection(g))
  g <- geos::geos_unnest(g, keep_multi = FALSE)
  g <- g[!geos::geos_is_empty(g)]
  if (length(g) == 0) return(sf::st_sfc(crs = sf::st_crs(geom)))

  starts <- geos::geos_point_start(g)
  ends <- geos::geos_point_end(g)
  starts_hex <- geos::geos_write_hex(starts)
  ends_hex <- geos::geos_write_hex(ends)
  point_counts <- table(c(starts_hex, ends_hex))
  junction_hex <- names(point_counts[point_counts >= 3])
  if (length(junction_hex) == 0) return(sf::st_sfc(crs = sf::st_crs(geom)))

  junction_points <- c(starts[starts_hex %in% junction_hex], ends[ends_hex %in% junction_hex])
  junction_points <- junction_points[!geos::geos_is_empty(junction_points)]
  if (length(junction_points) == 0) return(sf::st_sfc(crs = sf::st_crs(geom)))

  junction_points <- geos::geos_unary_union(geos::geos_make_collection(junction_points))
  junction_points <- geos::geos_unnest(junction_points, keep_multi = FALSE)
  sf::st_set_crs(sf::st_as_sfc(junction_points), sf::st_crs(geom))
}
