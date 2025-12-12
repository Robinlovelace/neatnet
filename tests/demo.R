
library(sf)
library(geos)
library(wk)

source("R/neatnet.R")

# Create a simple test case: two parallel lines close together
l1 <- matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)
l2 <- matrix(c(0, 5, 100, 5), ncol = 2, byrow = TRUE)

lines_sf <- st_sf(
  id = 1:2,
  geometry = st_sfc(st_linestring(l1), st_linestring(l2)),
  crs = 27700 # Projected CRS (British National Grid)
)

print("Original Network:")
print(lines_sf)

# Run neatnet
# dist = 5 (so buffer is 5m radius). The lines are 5m apart.
# They should merge into a single buffer of approx width 15m.
simplified <- neatnet(lines_sf, dist = 5, max_segment_length = 5)

print("Simplified Network:")
print(simplified)

# Check if we got geometry back
if (nrow(simplified) > 0) {
  print("Success: Geometry returned.")
  # Optional: check length or type
} else {
  stop("Failure: No geometry returned.")
}
