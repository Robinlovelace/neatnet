
library(sf)
library(geos)
source("R/neatnet.R")

# Parallel lines example
l1 <- matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)
l2 <- matrix(c(0, 5, 100, 5), ncol = 2, byrow = TRUE)
lines_sf <- st_sf(id = 1:2, geometry = st_sfc(st_linestring(l1), st_linestring(l2)), crs = 27700)

print("Running neatnet...")
simplified <- neatnet(lines_sf, dist = 5, max_segment_length = 5)

print(paste("Number of features returned:", nrow(simplified)))

# Check structure of intermediate steps (simulated)
# Copying logic from neat_skeletonize to debug
geometry <- as_geos_geometry(lines_sf)
buffered <- geos_buffer(geometry, 5, params = geos_buffer_params(quad_segs = 16))
geometry_union <- geos_unary_union(geos_make_collection(buffered))
boundary <- geos_boundary(geometry_union)
boundary_densified <- geos_densify(boundary, tolerance = 5)
voronoi_edges <- geos_voronoi_edges(boundary_densified)
skeleton_raw <- geos_intersection(voronoi_edges, geometry_union)
edges <- geos_unnest(skeleton_raw, keep_multi = FALSE)
dist_to_boundary <- geos_distance(edges, boundary)
is_spoke <- dist_to_boundary < 0.1
spine_edges <- edges[!is_spoke]

print(paste("Spine edges vector length:", length(spine_edges)))

# The suspected bug:
wrong_union <- geos_unary_union(spine_edges)
print(paste("Result of geos_unary_union(spine_edges) length:", length(wrong_union)))

# The fix:
collection <- geos_make_collection(spine_edges)
right_union <- geos_unary_union(collection)
print(paste("Result of geos_unary_union(collection) length:", length(right_union)))
print(paste("Result of geos_unary_union(collection) type:", geos_type(right_union)))

merged <- geos_line_merge(right_union)
print(paste("Merged type:", geos_type(merged)))
print(paste("Merged num geometries:", geos_num_geometries(merged)))
