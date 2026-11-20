# Spatial panel engine against splm (reference-only packages in C:/dev/reflibs; run in the full local R).
lib <- Sys.getenv("SP_REFLIB", "C:/dev/reflibs"); if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
suppressMessages(library(Matrix)); for (f in list.files("R", "^engine_.*[.]R$", full.names = TRUE)) source(f)
if (!requireNamespace("splm", quietly = TRUE)) { cat("splm not available: skipped\n"); quit(status = 0) }
suppressMessages({library(spdep); library(spatialreg); library(splm); library(plm)})
bad <- 0
chk <- function(label, a, b, tol = 1e-5) { ok <- isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))
    cat(sprintf("%-46s %s\n", label, if (ok) "ok" else sprintf("DIFF mine=%s ref=%s", paste(signif(as.numeric(a), 6)[1:min(4, length(a))], collapse = ","), paste(signif(as.numeric(b), 6)[1:min(4, length(b))], collapse = ",")))); if (!ok) bad <<- bad + 1 }
data(columbus, package = "spdep"); N <- 49; T <- 6
lw <- nb2listw(col.gal.nb); W <- as(as_dgRMatrix_listw(lw), "CsparseMatrix")
set.seed(5)
x1 <- matrix(rnorm(N * T), N, T); x2 <- matrix(rnorm(N * T), N, T) + matrix(rnorm(N), N, T)
mu <- rnorm(N, 0, 2); lam <- rnorm(T)
eps <- matrix(rnorm(N * T), N, T)
Y <- matrix(0, N, T)
A <- solve(diag(N) - 0.4 * as.matrix(W))
for (t in 1:T) Y[, t] <- as.numeric(A %*% (1 + 1.5 * x1[, t] - 0.8 * x2[, t] + mu + lam[t] + eps[, t]))
d <- data.frame(id = rep(1:N, T), time = rep(1:T, each = N), y = as.numeric(Y), x1 = as.numeric(x1), x2 = as.numeric(x2))
write.csv(d, "tests/panel_sim.csv", row.names = FALSE)
y <- d$y; X <- cbind(x1 = d$x1, x2 = d$x2)
f <- y ~ x1 + x2
for (eff in c("individual", "twoway")) {
    for (m in c("sar", "sem")) {
        fit <- sp_panel_fit(y, X, W, N, T, m, eff)
        ref <- spml(f, d, index = c("id", "time"), listw = lw, model = "within", effect = eff, lag = m == "sar", spatial.error = if (m == "sem") "b" else "none")
        sp <- if (m == "sar") ref$coefficients[1] else ref$spat.coef
        mine <- if (m == "sar") fit$rho else fit$lambda
        chk(sprintf("%s %s: spatial parameter", eff, m), mine, sp, 1e-2)
        bs <- coef(ref)[c("x1", "x2")]
        chk(sprintf("%s %s: beta", eff, m), fit$beta, bs, 1e-3)
        cat(sprintf("   ll mine %.4f  splm %.4f\n", fit$ll, if (is.null(ref$logLik)) NA else as.numeric(ref$logLik)))
    }
}
# pooled: identical to a cross-section SAR on the stacked data with the block-diagonal weights
fit <- sp_panel_fit(y, X, W, N, T, "sar", "none")
Wb <- kronecker(Diagonal(T), W); lwb <- mat2listw(as.matrix(Wb), style = "M")
ref <- lagsarlm(y ~ x1 + x2, d, lwb, method = "eigen")
chk("pooled sar: rho vs lagsarlm", fit$rho, ref$rho, 1e-4); chk("pooled sar: beta", fit$beta, coef(ref)[c(2, 3, 4)], 1e-4)
chk("pooled sar: logLik", fit$ll, as.numeric(logLik(ref)), 1e-6)
# SDM with individual effects: lag of the covariates, compared with splm's spml (Durbin)
fit <- sp_panel_fit(y, X, W, N, T, "sdm", "individual", lagcols = 1:2)
d$Wx1 <- as.numeric(as.matrix(W) %*% x1); d$Wx2 <- as.numeric(as.matrix(W) %*% x2)      # explicit WX regressors in a plain SAR
ref <- spml(y ~ x1 + x2 + Wx1 + Wx2, d, index = c("id", "time"), listw = lw, model = "within", effect = "individual", lag = TRUE, spatial.error = "none")
chk("individual sdm: rho", fit$rho, ref$coefficients[1], 1e-3); chk("individual sdm: beta", fit$beta, coef(ref)[c("x1", "x2", "Wx1", "Wx2")], 1e-3)
cat(if (bad == 0) "\nAll panel checks passed.\n" else sprintf("\n%d panel check(s) failed.\n", bad)); if (bad > 0) quit(status = 1)
