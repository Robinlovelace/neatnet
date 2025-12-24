
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
buffer_dist <- 20 # Search radius for candidates
frechet_thresh <- 15 # If Frechet dist < this, they are "the same road"

# 1. Convert to Geos
g <- as_geos_geometry(princes_st)
wk::wk_crs(g) <- NULL # Strip CRS to avoid mismatches
g_lines <- geos_unnest(g, keep_multi = FALSE) # Ensure linestrings

# 2. Candidate Search
tree <- geos_strtree(g_lines)
# Query candidates within buffer_dist
# We use a slightly larger box query then precise distance check if needed, 
# but strtree_query is just bbox.
bboxes <- geos_envelope(g_lines)
expanded_bboxes <- geos_create_rectangle(
  geos_xmin(bboxes) - buffer_dist,
  geos_ymin(bboxes) - buffer_dist,
  geos_xmax(bboxes) + buffer_dist,
  geos_ymax(bboxes) + buffer_dist
)
candidates <- geos_strtree_query(tree, expanded_bboxes)

# 3. Parallelism Check (Frechet)
# This can be N^2 in worst case local clusters, but efficient with index.
adj_list <- vector("list", length(g_lines))
count_parallel <- 0

cat("Classifying...\n")
pb <- txtProgressBar(min = 0, max = length(g_lines), style = 3)

for (i in seq_along(g_lines)) {
  nbrs <- candidates[[i]]
  nbrs <- nbrs[nbrs > i] # Only check forward to avoid double counting and self
  
  if (length(nbrs) == 0) {
    setTxtProgressBar(pb, i)
    next
  }
  
  # Filter by Fréchet Distance
  # geos_distance_frechet is expensive, so maybe pre-filter by bearing?
  # But let's try raw Frechet first on this smallish dataset (1000 features).
  
  # Optimization: Only check if bbox actually intersects (strtree is rough)
  # And maybe simple distance first?
  # But user specifically asked for Frechet.
  
  # Note: geos_distance_frechet computes similarity of shape.
  # If two lines are parallel but offset by 10m, Frechet dist ~ 10m.
  dists <- geos_distance_frechet(g_lines[i], g_lines[nbrs], densify = 0.5) 
  # densify=0.5 ensures we measure shape accurately even if vertices are sparse
  
  is_parallel <- dists < frechet_thresh
  
  if (any(is_parallel)) {
    matched <- nbrs[is_parallel]
    adj_list[[i]] <- matched
    # Also add back-links for graph
    # Actually simpler to build edge list for igraph
    count_parallel <- count_parallel + length(matched)
  }
  
  setTxtProgressBar(pb, i)
}
close(pb)

cat("\nFound", count_parallel, "parallel pairs.\n")

# 4. Grouping using Graph
# Build edge list
edges_mat <- do.call(rbind, lapply(seq_along(adj_list), function(i) {
  if (is.null(adj_list[[i]])) return(NULL)
  cbind(i, adj_list[[i]])
}))

if (!is.null(edges_mat)) {
  g_graph <- graph_from_edgelist(edges_mat, directed = FALSE)
  # Add all vertices
  g_graph <- add_vertices(g_graph, length(g_lines) - vcount(g_graph)) # Ensure all IDs exist
  cl <- components(g_graph)
  
  groups <- cl$membership
} else {
  groups <- seq_along(g_lines)
}

princes_st$group_id <- groups
princes_st$group_size <- ave(groups, groups, FUN = length)

cat("Group Sizes:\n")
print(table(princes_st$group_size))

# Visualize Results
plot(st_geometry(princes_st), col = ifelse(princes_st$group_size > 1, "red", "gray"), lwd = ifelse(princes_st$group_size > 1, 2, 1))
title("Parallel Groups (Red) vs Single (Gray)")

