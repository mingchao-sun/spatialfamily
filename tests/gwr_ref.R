# GWR engine against GWmodel (reference-only package in C:/dev/reflibs; run in the full local R).
lib <- Sys.getenv("SP_REFLIB", "C:/dev/reflibs"); if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
suppressMessages(library(Matrix)); for (f in list.files("R", "^engine_.*[.]R$", full.names = TRUE)) source(f)
if (!requireNamespace("GWmodel", quietly = TRUE)) { cat("GWmodel not available: skipped\n"); quit(status = 0) }
suppressMessages({library(GWmodel); library(sp)})
bad <- 0
chk <- function(label, a, b, tol = 1e-6) { ok <- isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))
    cat(sprintf("%-50s %s\n", label, if (ok) "ok" else sprintf("DIFF mine=%s ref=%s", paste(signif(as.numeric(a), 7)[1:min(3, length(a))], collapse = ","), paste(signif(as.numeric(b), 7)[1:min(3, length(b))], collapse = ",")))); if (!ok) bad <<- bad + 1 }
col <- read.csv("data/columbus.csv"); xy <- cbind(col$X, col$Y)
sdf <- SpatialPointsDataFrame(xy, col)
y <- col$CRIME; X <- cbind(1, col$INC, col$HOVAL)
for (kern in c("bisquare", "gaussian", "exponential", "tricube", "boxcar")) {
    for (ad in c(TRUE, FALSE)) {
        bw <- if (ad) 25 else 6
        ref <- suppressWarnings(gwr.basic(CRIME ~ INC + HOVAL, data = sdf, bw = bw, kernel = kern, adaptive = ad))
        fit <- sp_gwr(y, X, xy, FALSE, kern, ad, bw)
        lab <- sprintf("%s %s", kern, if (ad) "adaptive" else "fixed")
        chk(paste(lab, "coefficients"), fit$beta, as.matrix(ref$SDF@data[, c("Intercept", "INC", "HOVAL")]), 1e-6)
        chk(paste(lab, "standard errors"), fit$se, as.matrix(ref$SDF@data[, c("Intercept_SE", "INC_SE", "HOVAL_SE")]), 1e-6)
        chk(paste(lab, "AICc"), fit$aicc, ref$GW.diagnostic$AICc, 1e-8)
        chk(paste(lab, "R2"), fit$r2, ref$GW.diagnostic$gw.R2, 1e-8)
    }
}
ref <- suppressWarnings(gwr.basic(CRIME ~ INC + HOVAL, data = sdf, bw = 25, adaptive = TRUE))
fit <- sp_gwr(y, X, xy, FALSE, "bisquare", TRUE, 25)
chk("effective number of parameters", fit$enp, ref$GW.diagnostic$enp, 1e-8)
# bandwidth selection: AICc at the selected bandwidth must be no worse than GWmodel's
bwr <- suppressWarnings(bw.gwr(CRIME ~ INC + HOVAL, data = sdf, approach = "AICc", kernel = "bisquare", adaptive = TRUE))
sel <- sp_gwr(y, X, xy, FALSE, "bisquare", TRUE, NULL)
ra <- suppressWarnings(gwr.basic(CRIME ~ INC + HOVAL, data = sdf, bw = bwr, kernel = "bisquare", adaptive = TRUE))$GW.diagnostic$AICc
cat(sprintf("   bandwidth: mine %g (AICc %.4f), GWmodel %g (AICc %.4f)\n", sel$bw, sel$aicc, bwr, ra))
chk("selected bandwidth has AICc <= GWmodel's", sel$aicc <= ra + 1e-8, TRUE)
# lon/lat distances
nc <- read.csv("data/nc_sids.csv"); ll <- cbind(nc$lon, nc$lat)
sdl <- SpatialPointsDataFrame(ll, nc, proj4string = CRS("+proj=longlat +datum=WGS84"))
ref <- suppressWarnings(gwr.basic(sids_rate74 ~ nonwhite_share74, data = sdl, bw = 30, kernel = "bisquare", adaptive = TRUE, longlat = TRUE))
fit <- sp_gwr(nc$sids_rate74, cbind(1, nc$nonwhite_share74), ll, TRUE, "bisquare", TRUE, 30)
chk("great-circle distances: coefficients", fit$beta, as.matrix(ref$SDF@data[, c("Intercept", "nonwhite_share74")]), 5e-3)   # GWmodel computes great-circle distances slightly differently
cat(if (bad == 0) "\nAll GWR checks passed.\n" else sprintf("\n%d GWR check(s) failed.\n", bad)); if (bad > 0) quit(status = 1)
