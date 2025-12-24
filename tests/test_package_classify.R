
library(sf)
library(geos)
library(dplyr)

# Source the new function (mock package loading)
source("R/classify_network.R")

# Load Data
f <- system.file("extdata", "rnet_princes_street.geojson", package = "neatnet")
if (f == "") f <- "/home/robin/github/robinlovelace/neatnet/inst/extdata/rnet_princes_street.geojson"
princes_st <- st_read(f, quiet = TRUE) %>% st_transform(27700)

cat("Original Features:", nrow(princes_st), "\n")

# Run Classification
start_t <- Sys.time()
classified <- classify_network(princes_st, dist = 10, frechet_threshold = 15)
end_t <- Sys.time()

cat("\n--- Classification Results ---\n")
cat("Time elapsed:", as.numeric(end_t - start_t, units="secs"), "s\n")
print(table(classified$net_type))

cat("\nParallel Groups:", length(unique(na.omit(classified$parallel_group))), "\n")
cat("Clusters (Nodes):", length(unique(c(classified$cluster_u, classified$cluster_v))), "\n")
