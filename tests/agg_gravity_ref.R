# Agglomeration indices (hand-computed cases) and PPML gravity (against fixest::fepois).
lib <- Sys.getenv("SP_REFLIB", "C:/dev/reflibs"); if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
for (f in list.files("R", "^engine_.*[.]R$", full.names = TRUE)) source(f)
bad <- 0
chk <- function(label, a, b, tol = 1e-8) { ok <- isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))
    cat(sprintf("%-52s %s\n", label, if (ok) "ok" else sprintf("DIFF mine=%s ref=%s", paste(signif(as.numeric(a), 7)[1:min(3, length(a))], collapse = ","), paste(signif(as.numeric(b), 7)[1:min(3, length(b))], collapse = ",")))); if (!ok) bad <<- bad + 1 }
# --- agglomeration: 4 regions, 3 industries
E <- matrix(c(10, 20, 30,
              10, 20, 30,
              10, 20, 30,
              70, 0, 10), 4, 3, byrow = TRUE, dimnames = list(paste0("r", 1:4), c("A", "B", "C")))
ind <- sp_agg_industry(E)
tot <- sum(E); x <- rowSums(E) / tot
for (j in 1:3) {
    s <- E[, j] / sum(E[, j]); G <- sum((s - x)^2)
    chk(paste("raw concentration G, industry", colnames(E)[j]), ind$G[j], G)
    chk(paste("EG gamma with H = 0, industry", colnames(E)[j]), ind$gamma[j], G / (1 - sum(x^2)))
    chk(paste("HHI of regional shares, industry", colnames(E)[j]), ind$hhi[j], sum(s^2))
}
# an industry spread exactly like total employment has no concentration
E2 <- cbind(E, D = rowSums(E)); i2 <- sp_agg_industry(E2)
chk("proportional industry: G = 0, Gini = 0, gamma = 0", c(i2$G[4], i2$gini[4], i2$gamma[4]), c(0, 0, 0), 1e-12)
# an industry in a single region: Gini of the Lorenz curve by direct integration
E3 <- cbind(E, single = c(100, 0, 0, 0)); i3 <- sp_agg_industry(E3)
x3 <- rowSums(E3) / sum(E3); lq <- ifelse(x3 > 0, c(1, 0, 0, 0) / x3, NA); o <- order(lq)
X <- c(0, cumsum(x3[o])); Y <- c(0, cumsum(c(1, 0, 0, 0)[o]))
chk("Gini of a single-region industry (trapezoid rule)", i3$gini[4], 1 - sum(diff(X) * (Y[-1] + Y[-length(Y)])))
# plant Herfindahl: one plant per cell -> H = sum s^2 and gamma = (G - (1 - sum x^2) H) / ((1 - sum x^2)(1 - H))
Nf <- E * 0 + 1; Nf[E == 0] <- 0; ind2 <- sp_agg_industry(E, Nf)
s <- E[, 1] / sum(E[, 1]); chk("plant Herfindahl with one plant per cell", ind2$H[1], sum(s^2))
chk("EG gamma with the plant Herfindahl", ind2$gamma[1], (ind2$G[1] - (1 - sum(x^2)) * sum(s^2)) / ((1 - sum(x^2)) * (1 - sum(s^2))))
reg <- sp_agg_region(E)
s4 <- E[4, ] / sum(E[4, ]); rest <- colSums(E[1:3, ]) / sum(E[1:3, ])
chk("Krugman dissimilarity index of region 4", reg$krugman[4], sum(abs(s4 - rest)))
chk("regional HHI", reg$hhi[4], sum(s4^2))
lq <- sp_agg_lq(E); chk("location quotient (region 4, industry A)", lq[4, 1], (E[4, 1] / sum(E[4, ])) / (sum(E[, 1]) / tot))
# --- PPML gravity against fixest
if (requireNamespace("fixest", quietly = TRUE)) {
    suppressMessages(library(fixest))
    set.seed(21); no <- 25; nd <- 25
    g <- expand.grid(o = 1:no, d = 1:nd); g <- g[g$o != g$d, ]
    g$dist <- exp(rnorm(nrow(g), 5, 0.6)); g$contig <- rbinom(nrow(g), 1, 0.1); g$lgdp <- rnorm(nrow(g))
    mo <- rnorm(no); md <- rnorm(nd)
    g$flow <- rpois(nrow(g), exp(2 + mo[g$o] + md[g$d] - 0.9 * log(g$dist) + 0.5 * g$contig))
    g$flow[sample(nrow(g), 40)] <- 0
    g$ldist <- log(g$dist)
    X <- cbind(ldist = g$ldist, contig = g$contig)
    fit <- sp_ppml(g$flow, X, fe = list(g$o, g$d), cluster = NULL)
    ref <- fepois(flow ~ ldist + contig | o + d, g, vcov = "hetero")
    chk("PPML origin + destination FE: coefficients", fit$coef$est, coef(ref), 1e-6)
    ref2 <- fepois(flow ~ ldist + contig | o + d, g, vcov = ~o + d, ssc = ssc(adj = FALSE, cluster.adj = FALSE))
    fit2 <- sp_ppml(g$flow, X, fe = list(g$o, g$d), cluster = list(g$o, g$d))
    chk("PPML two-way clustered s.e. (no small-sample adjustment)", fit2$coef$se, se(ref2), 5e-2)
    g$pair <- as.integer(interaction(g$o, g$d)); pair <- g$pair
    ref3 <- fepois(flow ~ ldist + contig | o + d, g, vcov = ~pair, ssc = ssc(adj = FALSE, cluster.adj = FALSE))
    fit3 <- sp_ppml(g$flow, X, fe = list(g$o, g$d), cluster = list(pair))
    cat(sprintf("   pair-clustered se: mine %s | fixest %s\n", paste(signif(fit3$coef$se, 5), collapse = ","), paste(signif(se(ref3), 5), collapse = ",")))
    ref4 <- fepois(flow ~ ldist + contig, g, vcov = "hetero")
    fit4 <- sp_ppml(g$flow, X, fe = list())
    chk("PPML without fixed effects: coefficients", fit4$coef$est, coef(ref4), 1e-6)
    ref5 <- fepois(flow ~ ldist + contig | o, g, vcov = "hetero"); fit5 <- sp_ppml(g$flow, X, fe = list(g$o))
    chk("PPML origin FE only: coefficients", fit5$coef$est, coef(ref5), 1e-6)
    chk("PPML pseudo log-likelihood", fit$ll, as.numeric(logLik(fepois(flow ~ ldist + contig | o + d, g))), 1e-6)
    # a destination with only zero flows is dropped
    g$flow[g$d == 3] <- 0; fit6 <- sp_ppml(g$flow, X, fe = list(g$o, g$d))
    chk("separated group is dropped (at least its observations)", fit6$dropped >= sum(g$d == 3), TRUE)
} else cat("fixest not available: PPML checks skipped\n")
cat(if (bad == 0) "\nAll agglomeration / gravity checks passed.\n" else sprintf("\n%d check(s) failed.\n", bad)); if (bad > 0) quit(status = 1)
