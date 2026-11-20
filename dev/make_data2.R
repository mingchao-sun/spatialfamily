# Simulated example data for the agglomeration and gravity analyses (base R only). Run from the module folder.
set.seed(2026)
# ---- region x industry employment (25 regions, 6 industries)
R <- 25; regions <- sprintf("R%02d", 1:R)
size <- exp(rnorm(R, 0, 0.6))                                    # overall size of the regions
mk <- function(boost, base = 1) { w <- rep(base, R); w[boost] <- w[boost] * 12; w * size }
w <- list(Mining = mk(2, 0.2), Textiles = mk(c(5, 6, 7), 0.5), Software = mk(c(11, 12)), Autos = mk(c(15, 16, 17, 18), 0.8), Food = size, Retail = size)
tot <- c(Mining = 4000, Textiles = 30000, Software = 20000, Autos = 50000, Food = 60000, Retail = 120000)
rows <- do.call(rbind, lapply(names(w), function(ind) {
    e <- round(tot[[ind]] * w[[ind]] / sum(w[[ind]]) * exp(rnorm(R, 0, 0.1)))
    plant <- switch(ind, Mining = 400, Textiles = 90, Software = 60, Autos = 250, Food = 40, Retail = 12)
    data.frame(region = regions, industry = ind, employment = e, establishments = pmax(1, round(e / (plant * exp(rnorm(R, 0, 0.2))))))
}))
rows <- rows[rows$employment > 0, ]
write.csv(rows, "data/region_industry.csv", row.names = FALSE)
# ---- gravity: 20 countries, bilateral flows with zeros
N <- 20; ctry <- sprintf("C%02d", 1:N)
xy <- cbind(runif(N, 0, 10), runif(N, 0, 10)); g <- expand.grid(o = 1:N, d = 1:N); g <- g[g$o != g$d, ]
g$dist <- round(sqrt(rowSums((xy[g$o, ] - xy[g$d, ])^2)) * 400, 1) + 50
reg <- sample(1:4, N, TRUE); g$sameRegion <- as.integer(reg[g$o] == reg[g$d])
lang <- sample(1:5, N, TRUE); g$language <- as.integer(lang[g$o] == lang[g$d])
mo <- rnorm(N, 0, 0.7); md <- rnorm(N, 0, 0.7)
mu <- exp(3 + mo[g$o] + md[g$d] - 0.9 * log(g$dist) + 0.6 * g$sameRegion + 0.4 * g$language)
g$trade <- rpois(nrow(g), mu * 40)
g$origin <- ctry[g$o]; g$destination <- ctry[g$d]; g$ldist <- round(log(g$dist), 4)
write.csv(g[, c("origin", "destination", "trade", "dist", "ldist", "sameRegion", "language")], "data/gravity_sim.csv", row.names = FALSE)
cat(sprintf("%d region-industry rows, %d gravity rows (%d zero flows)\n", nrow(rows), nrow(g), sum(g$trade == 0)))
# ---- spatial DID: 100 units on a 100 x 100 km square, 10 years, 20 neighbouring units treated from year 6
set.seed(77); N <- 100; T <- 10
ux <- runif(N, 0, 100); uy <- runif(N, 0, 100)
centre <- c(30, 30); dc <- sqrt((ux - centre[1])^2 + (uy - centre[2])^2); treated_units <- order(dc)[1:20]
dd <- as.matrix(dist(cbind(ux, uy)))
sd <- expand.grid(unit = 1:N, year = 2001:2010)
sd$cx <- round(ux[sd$unit], 3); sd$cy <- round(uy[sd$unit], 3)
sd$treat <- as.integer(sd$unit %in% treated_units & sd$year >= 2006)
dmin <- vapply(1:N, function(i) min(dd[i, treated_units][treated_units != i]), 0)
ring1 <- as.integer(!(sd$unit %in% treated_units) & dmin[sd$unit] < 10 & sd$year >= 2006)
ring2 <- as.integer(!(sd$unit %in% treated_units) & dmin[sd$unit] >= 10 & dmin[sd$unit] < 20 & sd$year >= 2006)
sd$x1 <- round(rnorm(nrow(sd)), 3)
ai <- rnorm(N, 0, 2); gt <- cumsum(rnorm(10, 0.3, 0.4))
sd$outcome <- round(5 + ai[sd$unit] + gt[sd$year - 2000] + 2 * sd$treat + 1 * ring1 + 0.4 * ring2 + 0.5 * sd$x1 + rnorm(nrow(sd), 0, 1), 3)
write.csv(sd[order(sd$unit, sd$year), c("unit", "year", "outcome", "treat", "x1", "cx", "cy")], "data/spillover_did.csv", row.names = FALSE)
cat(sprintf("spillover_did: %d rows, %d ring-1 and %d ring-2 unit-periods\n", nrow(sd), sum(ring1), sum(ring2)))
