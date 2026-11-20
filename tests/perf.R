# Timing and sanity checks on larger samples (simulated spatial data on random locations).
suppressMessages({library(Matrix)}); for (f in list.files("R", "^engine_.*[.]R$", full.names = TRUE)) source(f)
set.seed(11)
for (n in c(1000, 2500, 5000)) {
    xy <- cbind(runif(n), runif(n))
    B <- sp_links_coords(xy, "knn", k = 6)$B; B <- sp_symmetrize(B); W <- sp_style(B, "W")
    X <- cbind(1, rnorm(n), rnorm(n))
    rho <- 0.5; A <- Diagonal(n) - rho * W
    y <- as.numeric(solve(A, X %*% c(1, 2, -1) + rnorm(n)))
    t0 <- proc.time()[3]; m <- sp_moran(y, W, "randomisation"); t1 <- proc.time()[3]
    li <- sp_local_moran(y, W, nperm = 199); t2 <- proc.time()[3]
    fit <- sp_sreg(y, X, W, "sar"); t3 <- proc.time()[3]
    ef <- sp_effects(fit, data.frame(b = 2:3, t = NA), W, nsim = 100); t4 <- proc.time()[3]
    cat(sprintf("n=%d  moran %.1fs (I=%.3f)  local moran 199 perm %.1fs  SAR %.1fs (rho=%.3f, b=%.3f,%.3f; %s)  effects(100 sims) %.1fs direct=%.3f\n",
                n, t1 - t0, m$stat, t2 - t1, t3 - t2, fit$rho, fit$beta[2], fit$beta[3], fit$ldmethod, t4 - t3, ef$est[1, 1]))
}
# traces: Hutchinson series versus exact eigenvalues
n <- 1200; xy <- cbind(runif(n), runif(n)); W <- sp_style(sp_symmetrize(sp_links_coords(xy, "knn", k = 6)$B), "W")
le <- sp_logdet_engine(W, TRUE); tr <- sp_traces(W, NULL, K = 60, probes = 40)
v <- data.frame(b = 2, t = 3); beta <- c(1, 2, 0.7)
for (rho in c(0.3, 0.6, 0.85)) {
    ex <- sp_effects_point(rho, beta, v, W, le); ap <- sp_effects_point(rho, beta, v, W, list(ev = NULL, rs = TRUE), tr)
    cat(sprintf("rho=%.2f exact direct=%.4f approx=%.4f (diff %.2f%%)\n", rho, ex[1, 1], ap[1, 1], 100 * (ap[1, 1] / ex[1, 1] - 1)))
}
