# Spatial DID engine: TWFE against fixest, rings by direct computation, Conley meat by brute force.
lib <- Sys.getenv("SP_REFLIB", "C:/dev/reflibs"); if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
suppressMessages(library(Matrix)); for (f in list.files("R", "^engine_.*[.]R$", full.names = TRUE)) source(f)
bad <- 0
chk <- function(label, a, b, tol = 1e-6) { ok <- isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))
    cat(sprintf("%-54s %s\n", label, if (ok) "ok" else sprintf("DIFF mine=%s ref=%s", paste(signif(as.numeric(a), 7)[1:min(4, length(a))], collapse = ","), paste(signif(as.numeric(b), 7)[1:min(4, length(b))], collapse = ",")))); if (!ok) bad <<- bad + 1 }
set.seed(8); N <- 60; T <- 8
xy <- cbind(runif(N, 0, 10), runif(N, 0, 10)); D <- as.matrix(dist(xy))
treated <- sample(N, 12); adopt <- rep(Inf, N); adopt[treated] <- sample(4:6, 12, TRUE)
d <- expand.grid(i = 1:N, t = 1:T); d$treat <- as.integer(d$t >= adopt[d$i])
TR <- matrix(d$treat, N, T)
dmin <- sp_did_nearest(D, TR)
# rings by direct computation
ring <- sp_did_rings(dmin, TR, c(1.5, 3))
chk_ring <- TRUE
for (t in 1:T) for (i in 1:N) { tr <- which(TR[, t] == 1); r <- 0
    if (TR[i, t] == 0 && length(tr)) { m <- min(D[i, tr]); r <- if (m < 1.5) 1 else if (m < 3) 2 else 0 }
    if (ring[i, t] != r) chk_ring <- FALSE }
chk("rings equal the direct computation", chk_ring, TRUE)
d$r1 <- as.numeric(as.vector(ring) == 1); d$r2 <- as.numeric(as.vector(ring) == 2)
d$x1 <- rnorm(nrow(d)); ai <- rnorm(N); gt <- rnorm(T)
d$y <- 1 + ai[d$i] + gt[d$t] + 2 * d$treat + 1.2 * d$r1 + 0.5 * d$r2 + 0.7 * d$x1 + rnorm(nrow(d))
X <- cbind(treat = d$treat, ring1 = d$r1, ring2 = d$r2, x1 = d$x1)
if (requireNamespace("fixest", quietly = TRUE)) {
    suppressMessages(library(fixest))
    fit <- sp_did_fit(d$y, X, d$i, d$t, "cluster")
    ref <- feols(y ~ treat + r1 + r2 + x1 | i + t, d, cluster = ~i)
    chk("coefficients", fit$coef$est, coef(ref), 1e-8)
    chk("clustered s.e. (fixest default small-sample correction)", fit$coef$se, se(ref), 1e-6)
    fit2 <- sp_did_fit(d$y, X, d$i, d$t, "robust")
    ref2 <- feols(y ~ treat + r1 + r2 + x1 | i + t, d, vcov = "hetero")
    cat(sprintf("   robust se: mine %s | fixest %s\n", paste(signif(fit2$coef$se, 5), collapse = ","), paste(signif(se(ref2), 5), collapse = ",")))
    fit3 <- sp_did_fit(d$y, X, d$i, d$t, "cluster", effects = "unit")
    ref3 <- feols(y ~ treat + r1 + r2 + x1 | i, d, cluster = ~i)
    chk("unit effects only: coefficients", fit3$coef$est, coef(ref3), 1e-6); chk("unit effects only: s.e.", fit3$coef$se, se(ref3), 1e-6)
} else cat("fixest not available\n")
# Conley meat against explicit double loops (lag 1)
ds <- d; Xs <- X
fit4 <- sp_did_fit(ds$y, Xs, ds$i, ds$t, "conley", D = D, bw = 3, lag = 1)
Z <- sp_absorb(cbind(ds$y, Xs), rep(1, nrow(ds)), list(ds$i, ds$t)); Xd <- Z[, -1]; e <- fit4$resid; sc <- Xd * e
M <- matrix(0, 4, 4)
for (a in seq_len(nrow(ds))) for (b in seq_len(nrow(ds))) {
    lg <- abs(ds$t[a] - ds$t[b]); if (lg > 1) next
    k <- max(0, 1 - D[ds$i[a], ds$i[b]] / 3) * (1 - lg / 2)
    M <- M + k * tcrossprod(sc[a, ], sc[b, ])
}
V <- solve(crossprod(Xd)) %*% ((M + t(M)) / 2) %*% solve(crossprod(Xd))
chk("Conley s.e. equal the double-loop definition", fit4$coef$se, sqrt(diag(V)), 1e-8)
cat(if (bad == 0) "\nAll spatial DID checks passed.\n" else sprintf("\n%d check(s) failed.\n", bad)); if (bad > 0) quit(status = 1)
