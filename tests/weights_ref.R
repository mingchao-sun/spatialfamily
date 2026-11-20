# Weights engine against spdep: contiguity from shapefile / GeoJSON, kNN, distance bands, orders, edge-list readers.
lib <- Sys.getenv("SP_REFLIB", "C:/dev/reflibs"); if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
suppressMessages(library(Matrix))
for (f in list.files("R", "^engine_.*[.]R$", full.names = TRUE)) source(f)
if (!requireNamespace("spdep", quietly = TRUE) || !requireNamespace("sf", quietly = TRUE)) { cat("spdep / sf not available: skipped\n"); quit(status = 0) }
suppressMessages({library(spdep); library(sf)})
bad <- 0
chk <- function(label, ok) { cat(sprintf("%-52s %s\n", label, if (isTRUE(ok)) "ok" else "DIFF")); if (!isTRUE(ok)) bad <<- bad + 1 }
nbmat <- function(nb) { m <- nb2mat(nb, style = "B", zero.policy = TRUE); m }
nc <- st_read("data/nc.shp", quiet = TRUE)
for (fmt in c("nc.shp", "nc.geojson")) {
    b <- sp_read_boundary(file.path("data", fmt))
    chk(paste(fmt, ": units"), length(b$units) == nrow(nc))
    chk(paste(fmt, ": id field"), identical(as.character(b$attrs$fips), as.character(nc$fips)))
    for (ty in c("queen", "rook")) {
        A <- sp_contig(b$units, ty)
        ref <- nbmat(poly2nb(nc, queen = ty == "queen"))
        chk(sprintf("%s %s contiguity identical to spdep", fmt, ty), all(as.matrix(A) == ref))
    }
}
b <- sp_read_boundary("data/nc.shp")
cen <- sp_centroids(b$units); ref <- suppressWarnings(st_coordinates(st_centroid(st_geometry(nc))))
chk("polygon centroids close to sf (planar)", max(abs(cen - ref)) < 0.05)
# order 2 (cumulative and pure)
A <- sp_contig(b$units, "queen"); nbq <- poly2nb(nc)
A2 <- sp_order(A, 2, FALSE); ref2 <- nbmat(nblag(nbq, 2)[[2]])
chk("order 2 (pure) identical to nblag", all(as.matrix(A2) == ref2))
A2c <- sp_order(A, 2, TRUE); ref2c <- nbmat(nblag_cumul(nblag(nbq, 2)))
chk("order 2 (cumulative) identical to nblag_cumul", all(as.matrix(A2c) == ref2c))
# coordinates
set.seed(3); xy <- cbind(runif(150), runif(150))
L <- sp_links_coords(xy, "knn", k = 6)$B
chk("kNN identical to knearneigh", all(as.matrix(L) == nbmat(knn2nb(knearneigh(xy, k = 6)))))
d <- sp_links_coords(xy, "dist", dband = 0.15)$B
chk("distance band identical to dnearneigh", all(as.matrix(d) == nbmat(dnearneigh(xy, 0, 0.15))))
da <- sp_links_coords(xy, "dist", dband = 0)
thr <- max(unlist(nbdists(knn2nb(knearneigh(xy, k = 1)), xy)))
chk("automatic threshold = largest nearest-neighbour distance", all(as.matrix(da$B) == nbmat(dnearneigh(xy, 0, thr * (1 + 1e-9)))))
ll <- cbind(runif(60, -100, -80), runif(60, 30, 45))
dk <- sp_links_coords(ll, "knn", k = 4, lonlat = TRUE)$B
chk("kNN with great-circle distance identical to spdep", all(as.matrix(dk) == nbmat(knn2nb(knearneigh(ll, k = 4, longlat = TRUE)))))
idw <- sp_links_coords(xy, "idw", power = 2, cutoff = 0.2)$B
D <- as.matrix(dist(xy)); Rm <- ifelse(D > 0 & D <= 0.2, D^-2, 0)
chk("inverse distance weights", isTRUE(all.equal(as.matrix(idw), Rm, check.attributes = FALSE)))
# full build with ids / alignment
d <- read.csv("data/nc_sids.csv"); d <- d[sample(nrow(d)), ]
s <- list(source = "boundary", file = "data/nc.shp", idfield = "fips", method = "queen", k = 5, dband = 0, power = 1, cutoff = 0, order = 1, cumul = FALSE, symm = FALSE, lonlat = FALSE, idvar = "fips", xvar = NULL, yvar = NULL)
w <- sp_build_weights(d, s)
chk("shuffled rows are aligned through the id", {
    ref <- nbmat(poly2nb(nc)); dimnames(ref) <- list(as.character(nc$fips), as.character(nc$fips))
    all(as.matrix(w$B) == ref[as.character(d$fips), as.character(d$fips)])
})
s$source <- "adjacency"; s$file <- "data/columbus_edges.csv"; cd <- read.csv("data/columbus.csv"); s$idvar <- "id"
w2 <- sp_build_weights(cd, s)
cs <- st_read("data/columbus.geojson", quiet = TRUE)
data(columbus, package = "spdep")
chk("edge list reproduces the textbook Columbus weights (col.gal.nb)", all(as.matrix(w2$B) == nbmat(col.gal.nb)))
s$source <- "boundary"; s$file <- "data/columbus.geojson"; s$idfield <- "id"; s$method <- "queen"
w3 <- sp_build_weights(cd, s); chk("GeoJSON queen contiguity of Columbus", all(as.matrix(w3$B) == nbmat(poly2nb(cs))))
# errors
e <- tryCatch(sp_build_weights(cd[-1, ], modifyList(s, list(idfield = "id"))), error = function(e) "err"); chk("extra units in the file are dropped", !identical(e, "err"))
e <- tryCatch(sp_build_weights(transform(cd, id = id + 1000), s), error = function(e) conditionMessage(e)); chk("unmatched ids give an error", grepl("no match", e))
# GAL / GWT round trip
gal <- tempfile(fileext = ".gal"); nbq <- poly2nb(cs); write.nb.gal(nbq, gal)
s$source <- "adjacency"; s$file <- gal; s$idvar <- "id"
w4 <- sp_build_weights(data.frame(id = seq_len(49)), s)
chk("GAL file read (ids 1..n)", all(as.matrix(w4$B) == nbmat(nbq)))
cat(if (bad == 0) "\nAll weights checks passed.\n" else sprintf("\n%d weights check(s) failed.\n", bad))
if (bad > 0) quit(status = 1)
