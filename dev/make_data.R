# Builds the example data in data/ from spdep / sf (reference-only packages in C:/dev/reflibs). Run in the full local R.
.libPaths(c("C:/dev/reflibs", .libPaths()))
suppressMessages({library(spdep); library(sf)})
out <- "data"
# Columbus, Ohio: 49 neighbourhoods (Anselin, 1988)
cs <- st_read(system.file("shapes/columbus.gpkg", package = "spData"), quiet = TRUE)
d <- st_drop_geometry(cs)
col <- data.frame(id = d$POLYID, CRIME = d$CRIME, INC = d$INC, HOVAL = d$HOVAL, DISCBD = d$DISCBD, OPEN = d$OPEN, PLUMB = d$PLUMB, X = d$X, Y = d$Y,
                  NSA = d$NSA, EW = d$EW, CP = d$CP)
write.csv(col, file.path(out, "columbus.csv"), row.names = FALSE, na = "")
cs$id <- cs$POLYID
st_write(cs[, "id"], file.path(out, "columbus.geojson"), delete_dsn = TRUE, quiet = TRUE)
# the textbook weights of the spdep columbus data (col.gal.nb); same row order as the polygons
ed <- data.frame(from = rep(cs$POLYID, card(col.gal.nb)), to = cs$POLYID[unlist(col.gal.nb)])
write.csv(ed, file.path(out, "columbus_edges.csv"), row.names = FALSE)
# North Carolina SIDS (Cressie, 1993): 100 counties
nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
cen <- suppressWarnings(st_coordinates(st_centroid(st_geometry(nc))))
ncd <- data.frame(fips = nc$FIPSNO, name = nc$NAME, lon = cen[, 1], lat = cen[, 2], BIR74 = nc$BIR74, SID74 = nc$SID74, NWBIR74 = nc$NWBIR74,
                  BIR79 = nc$BIR79, SID79 = nc$SID79, NWBIR79 = nc$NWBIR79)
ncd$sids_rate74 <- 1000 * ncd$SID74 / ncd$BIR74
ncd$nonwhite_share74 <- 100 * ncd$NWBIR74 / ncd$BIR74
write.csv(ncd, file.path(out, "nc_sids.csv"), row.names = FALSE)
nc$fips <- nc$FIPSNO
st_write(nc[, "fips"], file.path(out, "nc.geojson"), delete_dsn = TRUE, quiet = TRUE)
st_write(nc[, "fips"], file.path(out, "nc.shp"), delete_dsn = TRUE, quiet = TRUE)
cat(list.files(out), "\n")
