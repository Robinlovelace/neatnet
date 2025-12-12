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
#'   tolerance values in the topology check. Default 20.
#' @param grid_size Precision grid size used before topology checks.
#'   Defaults to `dist/5`.
#' @param parallel_tol_deg Maximum angular deviation (degrees) within a group to
#'   be considered parallel. Default 20.
#' @param min_parallel_lines Minimum number of features in a group for it to be
#'   considered parallel. Default 2.
#' @param complex_dist Distance (map units) used for the local-neighbourhood
#'   check that flags complex areas. Defaults to `min(dist, 8)` so increasing
#'   `dist` for parallel detection does not automatically make everything
#'   "complex".
#' @param complex_min_neighbors Minimum number of nearby neighbours (within
#'   `complex_dist`) required for a feature to be considered in a locally complex area.
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
  dist = 20,
  grid_size = dist / 5,
  parallel_tol_deg = 20,
  min_parallel_lines = 2,
  complex_dist = min(dist, 8),
  complex_min_neighbors = 8,
  complex_min_disp_deg = 60
) {
  if (inherits(x, "sf")) {
    sf_x <- x
    geom_sfc <- sf::st_geometry(sf_x)
  } else if (inherits(x, "sfc")) {
    geom_sfc <- x
    sf_x <- sf::st_sf(geometry = geom_sfc)
  } else if (inherits(x, "geos_geometry")) {
    geom_sfc <- sf::st_as_sfc(x)
    sf_x <- sf::st_sf(geometry = geom_sfc)
  } else {
    stop("`x` must be an sf, sfc, or geos_geometry object.")
  }

  geom <- geom_sfc
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

  # 1) Neighbourhood graph (within distance) using geos STRtree.
  # We first query candidates using expanded bboxes (fast), then refine using
  # exact GEOS distances.
  g <- geos::as_geos_geometry(geom)
  g <- geos::geos_unnest(g, keep_multi = FALSE)
  if (length(g) != length(geom)) {
    stop("Input must contain only LINESTRING/MULTILINESTRING and no empty geometries.")
  }
  if (any(geos::geos_is_empty(g))) {
    stop("Input contains empty geometries; please drop empties before classifying.")
  }

  tree <- geos::geos_strtree(g)
  srid <- geos::geos_srid(g)
  srid <- if (length(srid) > 0) srid[[1]] else NA_integer_
  crs <- wk::wk_crs(g)

  .neatnet_hits_within <- function(distance) {
    rects <- geos::geos_create_rectangle(
      geos::geos_xmin(g) - distance,
      geos::geos_ymin(g) - distance,
      geos::geos_xmax(g) + distance,
      geos::geos_ymax(g) + distance
    )
    if (!is.null(crs)) {
      rects <- wk::wk_set_crs(rects, crs)
    }
    if (!is.na(srid)) {
      rects <- geos::geos_set_srid(rects, srid)
    }
    cand <- geos::geos_strtree_query(tree, rects)
    hits <- vector("list", length(g))
    for (i in seq_along(cand)) {
      idx <- cand[[i]]
      if (length(idx) == 0) {
        hits[[i]] <- integer(0)
        next
      }
      d <- geos::geos_distance(g[idx], g[i])
      keep <- as.numeric(d) <= distance
      keep[is.na(keep)] <- FALSE
      hits[[i]] <- idx[keep]
    }
    hits
  }

  hits_parallel <- .neatnet_hits_within(dist)
  hits_complex <- .neatnet_hits_within(complex_dist)

  # 2) Feature-level bearing (0..180, direction-invariant)
  sf_x$bearing_deg <- vapply(seq_along(g), function(i) .neatnet_bearing_deg_180(g[i]), numeric(1))

  # 3) Feature-level parallel flag: has at least one near-parallel neighbour
  # within the neighbourhood.
  feature_parallel <- logical(length(geom))
  parallel_adj <- vector("list", length(geom))
  if (length(geom) >= 2) {
    for (i in seq_along(hits_parallel)) {
      nbrs <- hits_parallel[[i]]
      nbrs <- setdiff(nbrs, i)
      if (length(nbrs) == 0) next

      # Exclude any geometries that intersect/touch this feature so we don't
      # label connected/collinear segments as "parallel".
      it <- geos::geos_intersects(g[i], g[nbrs])
      it <- as.logical(it)
      it[is.na(it)] <- FALSE
      nbrs <- nbrs[!it]
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
    for (i in seq_along(hits_complex)) {
      nbrs <- hits_complex[[i]]
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
  # Handles LINESTRING/MULTILINESTRING.
  # Uses an overall principal-axis bearing (PCA of sampled points) so curvy
  # carriageways still get a stable direction.
  # Returns degrees in [0, 180), so opposite directions are considered parallel.
  if (!inherits(geom, "geos_geometry")) {
    geom <- geos::as_geos_geometry(geom)
  }
  if (length(geom) != 1 || geos::geos_is_empty(geom)) return(NA_real_)

  L <- as.numeric(geos::geos_length(geom))
  if (!is.finite(L) || L <= 0) return(NA_real_)

  # Sample points along the line and compute the principal axis.
  # Keep this small for speed; we only need a stable orientation.
  n <- 7L
  dists <- seq(0, L, length.out = n)
  pts <- vector("list", length(dists))
  for (i in seq_along(dists)) {
    pts[[i]] <- geos::geos_interpolate(geom, dists[[i]])
  }
  pts <- do.call(c, pts)
  pts <- pts[!geos::geos_is_empty(pts)]
  if (length(pts) < 2) {
    # Fallback to start/end.
    p0 <- geos::geos_point_start(geom)
    p1 <- geos::geos_point_end(geom)
    if (geos::geos_is_empty(p0) || geos::geos_is_empty(p1)) return(NA_real_)
    x0 <- as.numeric(geos::geos_x(p0)); y0 <- as.numeric(geos::geos_y(p0))
    x1 <- as.numeric(geos::geos_x(p1)); y1 <- as.numeric(geos::geos_y(p1))
    dx <- x1 - x0; dy <- y1 - y0
    if (!is.finite(dx) || !is.finite(dy) || (dx == 0 && dy == 0)) return(NA_real_)
    ang <- atan2(dy, dx) * 180 / pi
    ang <- ang %% 180
    if (ang < 0) ang <- ang + 180
    return(ang)
  }

  x <- as.numeric(geos::geos_x(pts))
  y <- as.numeric(geos::geos_y(pts))
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NA_real_)

  x <- x - mean(x)
  y <- y - mean(y)
  if (all(x == 0) && all(y == 0)) return(NA_real_)

  Sxx <- mean(x * x)
  Syy <- mean(y * y)
  Sxy <- mean(x * y)
  cov <- matrix(c(Sxx, Sxy, Sxy, Syy), nrow = 2)
  ev <- eigen(cov, symmetric = TRUE)
  v <- ev$vectors[, which.max(ev$values)]
  dx <- v[[1]]
  dy <- v[[2]]
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
