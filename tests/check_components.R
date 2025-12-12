library(sf)
library(geos)
source("R/neatnet.R")

f <- system.file("extdata", "rnet_princes_street.geojson", package = "neatnet")
if (f == "") f <- "inst/extdata/rnet_princes_street.geojson"
princes_st <- st_read(f, quiet = TRUE)
princes_st <- st_transform(princes_st, 27700)

simplified <- neatnet(princes_st, dist = 8, max_segment_length = 5)

lengths <- st_length(simplified)
print(summary(lengths))
print(paste("Num features < 10m:", sum(lengths < units::set_units(10, m))))
print(paste("Num features < 5m:", sum(lengths < units::set_units(5, m))))