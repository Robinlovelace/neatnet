library(sf)
library(geos)
library(neatnet)

# Create a simple test case: two parallel lines
l1 <- matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)
l2 <- matrix(c(0, 5, 100, 5), ncol = 2, byrow = TRUE)

lines_sf <- st_sf(
  id = 1:2,
  geometry = st_sfc(st_linestring(l1), st_linestring(l2)),
  crs = 27700 # Projected CRS
)

print("Original Network:")
print(lines_sf)

# Run neatnet with default parameters
simplified <- neatnet(lines_sf, dist = 5) # max_segment_length and final_min_length are derived from dist

print("Simplified Network:")
print(simplified)

# Check if we got geometry back
if (nrow(simplified) > 0) {
  print("Success: Geometry returned.")
  # Optional: check length or type
} else {
  stop("Failure: No geometry returned.")
}