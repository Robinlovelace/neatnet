
library(sf)
library(geos)
library(dplyr)
library(igraph)

# Load Data
f <- system.file("extdata", "rnet_princes_street.geojson", package = "neatnet")
if (f == "") f <- "/home/robin/github/robinlovelace/neatnet/inst/extdata/rnet_princes_street.geojson"
princes_st <- st_read(f, quiet = TRUE) %>% st_transform(27700)

cat("Features:", nrow(princes_st), "\n")

# Parameters
dist <- 10        # Buffer radius (approx 20m width junction tolerance)
frechet_thresh <- 15 

# 1. Geometry Prep
g <- as_geos_geometry(princes_st)
wk::wk_crs(g) <- NULL
g_lines <- geos_unnest(g, keep_multi = FALSE)

# 2. Node Clustering (Identify Complex Intersections)
# Extract all endpoints
starts <- geos_point_start(g_lines)
ends <- geos_point_end(g_lines)
all_points <- c(starts, ends)

# Cluster nodes: Buffer, Union, connected components
# Since we want to merge nodes that are close, we buffer by dist/2
node_buffs <- geos_buffer(all_points, dist / 2)
node_clusters_geom <- geos_unary_union(geos_make_collection(node_buffs))
node_clusters_poly <- geos_unnest(node_clusters_geom, keep_multi = FALSE) 
# Note: geos_unnest on MPOLY describes connected components

# Create Cluster Index
# We need to map each original Point to a Cluster ID
# Quick way: Spatial Join (Intersects). Or Centroids.
# Build spatial index on clusters
cluster_tree <- geos_strtree(node_clusters_poly)

# Map Start/End points to Cluster IDs
get_cluster_id <- function(pts, tree, poly) {
  # Query tree
  # pts is vector of points
  # strictly, each point should fall in exactly one cluster polygon
  # query might return multiple if overlaps (shouldn't happen after union)
  # But union might produce touching polys?
  # Let's assume unique containment.
  
  # geos_strtree_query returns candidates. We iterate or verify.
  # Optimization: Use Centroid distance?
  # Or simply: intersection check.
  
  # Since we buffered points, the point IS inside the polygon.
  # geos_strtree_query is fast.
  q <- geos_strtree_query(tree, pts)
  vapply(q, function(x) if(length(x)>0) as.integer(x[1]) else NA_integer_, integer(1))
}

start_cluster <- get_cluster_id(starts, cluster_tree, node_clusters_poly)
end_cluster <- get_cluster_id(ends, cluster_tree, node_clusters_poly)

# 3. Graph Construction
# Edges now connect (StartCluster, EndCluster)
edges_df <- data.frame(
  edge_id = seq_along(g_lines),
  u = start_cluster,
  v = end_cluster
) %>%
  rowwise() %>%
  mutate(
    pair_id = paste(sort(c(u, v)), collapse = "-"),
    is_loop = u == v
  ) %>%
  ungroup()

# 4. Analyze Edges
# Loops: Edges entirely within a cluster -> Collapse to point
loops <- edges_df %>% filter(is_loop)

# Parallel Candidates: Multiple edges between same (u, v) (where u != v)
parallel_groups <- edges_df %>%
  filter(!is_loop) %>%
  group_by(pair_id) %>%
  filter(n() > 1) %>% # At least 2 edges
  mutate(group_id = cur_group_id()) %>%
  ungroup()

cat("\n--- Graph Analysis ---\n")
cat("Total Clusters (Nodes):", length(node_clusters_poly), "\n")
cat("Internal Edges (Loops to collapse):", nrow(loops), "\n")
cat("Parallel Candidate Groups:", length(unique(parallel_groups$group_id)), "\n")
cat("Edges in Parallel Groups:", nrow(parallel_groups), "\n")

# 5. Verify Parallelism (Frechet) within Candidate Groups
# Only check if shapes are similar.
# Since they share endpoints (Cluster-wise), they form a "Parallel Link".
# We can refine this: if one path is straight and another is very long (detour), maybe keep both?
# Frechet distance check against the "shortest" edge in the group.

verified_parallel <- 0

# Check Frechet
# Iterate groups
groups <- unique(parallel_groups$group_id)
final_parallel_edges <- integer(0)

if (length(groups) > 0) {
  cat("\nVerifying Frechet Distances...\n")
  pb <- txtProgressBar(min = 0, max = length(groups), style = 3)
  
  for (i in seq_along(groups)) {
    gid <- groups[i]
    ids <- parallel_groups$edge_id[parallel_groups$group_id == gid]
    
    # Extract geometries
    geoms <- g_lines[ids]
    
    # Pick reference (shortest length)
    lens <- geos_length(geoms)
    ref_idx <- which.min(lens)
    ref_geom <- geoms[ref_idx]
    
    # Calc Frechet from reference
    fd <- geos_distance_frechet(geoms, ref_geom, densify = 0.5)
    
    # Keep those within threshold relative to reference
    # Note: ref itself has dist 0
    is_par <- fd < frechet_thresh
    
    # If multiple are parallel, we mark them for merging
    if (sum(is_par) > 1) {
      final_parallel_edges <- c(final_parallel_edges, ids[is_par])
      verified_parallel <- verified_parallel + sum(is_par)
    }
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
}

cat("\nVerified Parallel Edges:", length(final_parallel_edges), "\n")

# Summary
non_parallel <- nrow(princes_st) - length(final_parallel_edges) - nrow(loops)
cat("Single Edges (Keep):", non_parallel, "\n")

