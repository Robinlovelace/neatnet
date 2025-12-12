
library(sf)
library(geos)
library(wk)

# 1. Define Input
l1 <- matrix(c(0, 0, 100, 0), ncol = 2, byrow = TRUE)
l2 <- matrix(c(0, 5, 100, 5), ncol = 2, byrow = TRUE)

lines_sf <- st_sf(
  id = 1:2,
  geometry = st_sfc(st_linestring(l1), st_linestring(l2)),
  crs = 27700
)

# 2. Run R neatnet
source("R/neatnet.R")
out_r_path <- "tests/output_neatnet_r.geojson"
tryCatch({
  simplified_r <- neatnet(lines_sf, dist = 5, max_segment_length = 5)
  st_write(simplified_r, out_r_path, delete_dsn = TRUE, quiet = TRUE)
  print(paste("R output saved to", out_r_path))
}, error = function(e) {
  print(paste("R neatnet failed:", e$message))
})

# 3. Run Python Comparison (neatnet & parenx)
# We assume 'tests/compare_python.py' exists and venv is set up
cmd <- "bash -c 'source venv/bin/activate && python3 tests/compare_python.py'"
print("Running Python comparison script...")
system(cmd, intern = FALSE)

# 4. Compare Results
print("--- Comparison Results ---")

read_safe <- function(path) {
  if (file.exists(path)) {
    return(st_read(path, quiet = TRUE))
  } else {
    return(NULL)
  }
}

sf_r <- read_safe(out_r_path)
sf_py_neat <- read_safe("tests/output_neatnet_python.geojson")
sf_py_par <- read_safe("tests/output_parenx.geojson")

compare_geoms <- function(name1, sf1, name2, sf2) {
  if (is.null(sf1) || is.null(sf2)) {
    print(paste("Cannot compare", name1, "and", name2, "- one is missing."))
    return()
  }
  
  # Simple check: Bounding box overlap and length similarity
  bb1 <- st_bbox(sf1)
  bb2 <- st_bbox(sf2)
  
  len1 <- as.numeric(sum(st_length(sf1)))
  len2 <- as.numeric(sum(st_length(sf2)))
  
  print(paste("Comparing", name1, "vs", name2))
  print(paste("  Total Length:", round(len1, 2), "vs", round(len2, 2)))
  
  # Hausdorff distance (max distance between them)
  # We unite them to single multilinestrings for comparison
  g1 <- st_union(sf1)
  g2 <- st_union(sf2)
  
  hd <- st_distance(g1, g2, which = "Hausdorff")
  print(paste("  Hausdorff Distance:", round(as.numeric(hd), 3)))
  
  # Visual check description
  if (as.numeric(hd) < 2.0) {
    print("  -> Result: VERY SIMILAR")
  } else if (as.numeric(hd) < 10.0) {
    print("  -> Result: SOMEWHAT SIMILAR")
  } else {
    print("  -> Result: DIFFERENT")
  }
}

compare_geoms("R (neatnet)", sf_r, "Python (neatnet)", sf_py_neat)
compare_geoms("R (neatnet)", sf_r, "Python (parenx)", sf_py_par)
compare_geoms("Python (neatnet)", sf_py_neat, "Python (parenx)", sf_py_par)

# 5. Check for "sameness" or "best of" validity
# Ideally, for parallel lines 5m apart with 5m buffer (radius), 
# they should merge into a single centerline approx at y=2.5.
# Length should be approx 100.

check_quality <- function(name, sf) {
  if (is.null(sf)) return()
  len <- as.numeric(sum(st_length(sf)))
  # Expected length around 100 (merged) or 200 (not merged)
  # If merged, it should be near 100.
  
  print(paste("Quality Check for", name))
  if (len < 120 && len > 80) {
    print("  -> Merged successfully (Length approx 100)")
  } else if (len >= 180) {
    print("  -> Did NOT merge (Length approx 200)")
  } else {
    print(paste("  -> Partial or weird merge (Length", round(len, 2), ")"))
  }
}

check_quality("R (neatnet)", sf_r)
check_quality("Python (neatnet)", sf_py_neat)
check_quality("Python (parenx)", sf_py_par)
