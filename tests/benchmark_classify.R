# Benchmark: neatnet_classify_groups

suppressPackageStartupMessages({
  library(sf)
  library(neatnet)
})

f <- system.file("extdata", "rnet_princes_street.geojson", package = "neatnet")
princes_st <- st_transform(st_read(f, quiet = TRUE), 27700)

cat("n features:", nrow(princes_st), "\n")

# Warm up
invisible(neatnet_classify_groups(princes_st, dist = 8))

bm <- system.time({
  groups <- neatnet_classify_groups(princes_st, dist = 8)
})

print(bm)
print(table(groups$group_class))

elapsed <- unname(bm[["elapsed"]])
if (is.finite(elapsed) && elapsed > 1) {
  warning(sprintf("Classification took %.2fs (>1s)", elapsed))
}
