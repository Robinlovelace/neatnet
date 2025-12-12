
library(geos)
library(sf)

# Recreate the scenario
l1 <- matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)
l2 <- matrix(c(0, 5, 100, 5), ncol = 2, byrow = TRUE)
lines_sf <- st_sf(id = 1:2, geometry = st_sfc(st_linestring(l1), st_linestring(l2)), crs = 27700)
dist <- 5
max_segment_length <- 5

geometry <- as_geos_geometry(lines_sf)
buffered <- geos_buffer(geometry, dist, params = geos_buffer_params(quad_segs = 16))
collection <- geos_make_collection(buffered)
geometry_union <- geos_unary_union(collection)

boundary <- geos_boundary(geometry_union)
boundary_densified <- geos_densify(boundary, tolerance = max_segment_length)
voronoi_edges <- geos_voronoi_edges(boundary_densified)
skeleton_raw <- geos_intersection(voronoi_edges, geometry_union)

edges <- geos_unnest(skeleton_raw, keep_multi = FALSE)
is_spoke <- geos_intersects(edges, boundary)
spine_edges <- edges[!is_spoke]

dists <- geos_distance(spine_edges, boundary)
print(paste("Max distance of remaining edges to boundary:", max(dists)))
print(paste("Min distance of remaining edges to boundary:", min(dists)))
print(paste("Number of edges with dist < 1e-6:", sum(dists < 1e-6)))
