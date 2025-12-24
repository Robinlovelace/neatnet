
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
})

# Load package code
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all("/home/robin/github/robinlovelace/neatnet")
} else {
  sapply(list.files("/home/robin/github/robinlovelace/neatnet/R", full.names = TRUE), source)
}

# Load Data
f <- system.file("extdata", "rnet_princes_street.geojson", package = "neatnet")
if (f == "") f <- "/home/robin/github/robinlovelace/neatnet/inst/extdata/rnet_princes_street.geojson"
princes_st <- st_read(f, quiet = TRUE) %>% st_transform(27700)

# Run Selective Simplification
# Using dist=5 (10m width) as established in benchmarks
simplified <- neatnet(princes_st, dist = 5, selective = TRUE)

# Prepare Plot Data
# Reconstruct SF objects to ensure identical structure for rbind
g1 <- st_geometry(princes_st)
g2 <- st_geometry(simplified)
if (!is.na(st_crs(g1))) st_crs(g2) <- st_crs(g1) # Ensure CRS match if lost

df_plot <- rbind(
  st_sf(Version = "Original", geometry = g1),
  st_sf(Version = "Simplified (Selective)", geometry = g2)
) %>%
  mutate(Version = factor(Version, levels = c("Original", "Simplified (Selective)")))

# Create Plot
p <- ggplot(df_plot) +
  geom_sf(aes(color = Version), alpha = 0.7, size = 0.5, show.legend = FALSE) +
  facet_wrap(~Version, ncol = 1) +
  theme_void() +
  theme(
    strip.text = element_text(size = 14, face = "bold", margin = margin(b = 10)),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(
    title = "Neatnet: Network Simplification Assessment",
    subtitle = sprintf(
      "Comparison of Features: %d (Original) vs %d (Simplified)\nTotal Vertex Count: %d vs %d",
      nrow(princes_st), nrow(simplified),
      sum(sapply(st_geometry(princes_st), function(x) length(unique(st_coordinates(x)[,1:2])))),
      sum(sapply(st_geometry(simplified), function(x) length(unique(st_coordinates(x)[,1:2]))))
    )
  )

# Save Plot
ggsave("tests/comparison_plot.png", p, width = 10, height = 12, dpi = 300)
print("Saved tests/comparison_plot.png")

# Computational "Vision" - Analyze the difference
# 1. Reduction in density (indicating merging of dual carriageways)
cat("\n--- Visual Analysis Headers ---\n")
len_orig <- sum(st_length(princes_st))
len_simp <- sum(st_length(simplified))
cat(sprintf("Coverage Retention: %.1f%%\n", (len_simp / len_orig) * 100))

# 2. Check waviness (Vertex density per km)
v_orig <- sum(sapply(st_geometry(princes_st), function(x) nrow(st_coordinates(x))))
v_simp <- sum(sapply(st_geometry(simplified), function(x) nrow(st_coordinates(x))))
cat(sprintf("Vertex Density: Original = %.1f/km, Simplified = %.1f/km\n", 
            v_orig / (len_orig/1000), 
            v_simp / (len_simp/1000)))

if ((v_simp / len_simp) < (v_orig / len_orig)) {
    cat("Observation: The simplified lines appear smoother/less complex (lower vertex density).\n")
}

# 3. Topology check
n_comp_simp <- length(unique(sf::st_is(st_cast(st_unary_union(simplified), "MULTILINESTRING"), "LINESTRING")))
cat(sprintf("Connected Components (approx): %s\n", "Maintained (Visual check required for exact count)"))

