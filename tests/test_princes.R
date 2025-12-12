
library(sf)
library(geos)
source("R/neatnet.R")

f <- system.file("extdata", "rnet_princes_street.geojson", package = "neatnet")
if (f == "") f <- "inst/extdata/rnet_princes_street.geojson"
princes_st <- st_read(f, quiet = TRUE)
princes_st <- st_transform(princes_st, 27700)

print(paste("Original features:", nrow(princes_st)))

# Run neatnet with default parameters
simplified <- neatnet(princes_st, dist = 8) 
print(paste("Simplified features:", nrow(simplified)))
