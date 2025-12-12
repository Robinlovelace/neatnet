library(sf)
library(geos)
library(neatnet)

f <- system.file("extdata", "rnet_princes_street.geojson", package = "neatnet")
if (f == "") f <- "inst/extdata/rnet_princes_street.geojson"
princes_st <- st_read(f, quiet = TRUE)
princes_st <- st_transform(princes_st, 27700)

simplified <- neatnet(princes_st, dist = 8)

lengths_m <- as.numeric(st_length(simplified))
print(summary(lengths_m))
print(paste("Num features < 10m:", sum(lengths_m < 10)))
print(paste("Num features < 5m:", sum(lengths_m < 5)))