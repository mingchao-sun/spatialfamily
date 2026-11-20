# Spatial autocorrelation statistics (global and local). Base R and Matrix only.
# W is a sparse n x n weights matrix (already styled); x a numeric vector.

sp_S <- function(W) {
    W <- as(W, "CsparseMatrix")
    S0 <- sum(W)
    Wt <- Matrix::t(W)
    S1 <- 0.5 * sum((W + Wt)^2)
    S2 <- sum((Matrix::rowSums(W) + Matrix::colSums(W))^2)
    list(S0 = S0, S1 = S1, S2 = S2)
}

sp_test_p <- function(z, alt) {
    switch(alt, two.sided = 2 * stats::pnorm(-abs(z)), greater = stats::pnorm(z, lower.tail = FALSE), less = stats::pnorm(z))
}

# Permutation p-value from the observed statistic and the simulated ones (pseudo p-value, with the observed value counted).
sp_perm_p <- function(obs, sims, alt) {
    ns <- length(sims)
    ge <- sum(sims >= obs) + 1; le <- sum(sims <= obs) + 1
    switch(alt, greater = ge / (ns + 1), less = le / (ns + 1), two.sided = min(1, 2 * min(ge, le) / (ns + 1)))
}

# ------------------------------------------------------------------ Moran's I --

sp_moran <- function(x, W, inference = "randomisation", alt = "two.sided", nperm = 999, S = sp_S(W)) {
    n <- length(x); S0 <- S$S0; S1 <- S$S1; S2 <- S$S2
    z <- x - mean(x); m2 <- sum(z^2)
    if (m2 == 0) stop("The variable has no variation.", call. = FALSE)
    Wz <- as.numeric(W %*% z)
    I <- (n / S0) * sum(z * Wz) / m2
    EI <- -1 / (n - 1)
    out <- list(stat = I, expect = EI, S = S)
    if (inference %in% c("randomisation", "normality")) {
        if (inference == "normality") {
            V <- (n^2 * S1 - n * S2 + 3 * S0^2) / ((n^2 - 1) * S0^2) - EI^2
        } else {
            b2 <- n * sum(z^4) / m2^2
            V <- (n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * S0^2) - b2 * ((n^2 - n) * S1 - 2 * n * S2 + 6 * S0^2)) /
                ((n - 1) * (n - 2) * (n - 3) * S0^2) - EI^2
        }
        out$var <- V; out$z <- (I - EI) / sqrt(V); out$p <- sp_test_p(out$z, alt)
    } else {
        sims <- vapply(seq_len(nperm), function(b) {
            zz <- z[sample.int(n)]
            (n / S0) * sum(zz * as.numeric(W %*% zz)) / m2
        }, 0)
        out$sims <- sims; out$var <- stats::var(sims); out$z <- (I - mean(sims)) / stats::sd(sims); out$p <- sp_perm_p(I, sims, alt)
    }
    out
}

# ------------------------------------------------------------------ Geary's C --

sp_geary <- function(x, W, inference = "randomisation", alt = "two.sided", nperm = 999, S = sp_S(W)) {
    n <- length(x); S0 <- S$S0; S1 <- S$S1; S2 <- S$S2
    z <- x - mean(x); m2 <- sum(z^2)
    if (m2 == 0) stop("The variable has no variation.", call. = FALSE)
    cfun <- function(v) {
        zz <- v - mean(v); Wz <- as.numeric(W %*% v)
        # sum_ij w_ij (v_i - v_j)^2 = sum_i v_i^2 r_i + sum_j v_j^2 c_j - 2 sum_ij w_ij v_i v_j
        num <- sum(v^2 * Matrix::rowSums(W)) + sum(v^2 * Matrix::colSums(W)) - 2 * sum(v * Wz)
        (n - 1) * num / (2 * S0 * sum(zz^2))
    }
    C <- cfun(x)
    out <- list(stat = C, expect = 1)
    if (inference %in% c("randomisation", "normality")) {
        if (inference == "normality") {
            V <- ((2 * S1 + S2) * (n - 1) - 4 * S0^2) / (2 * (n + 1) * S0^2)
        } else {
            k <- n * sum(z^4) / m2^2
            V <- ((n - 1) * S1 * (n^2 - 3 * n + 3 - (n - 1) * k) - 0.25 * (n - 1) * S2 * (n^2 + 3 * n - 6 - (n^2 - n + 2) * k) +
                  S0^2 * (n^2 - 3 - (n - 1)^2 * k)) / (n * (n - 2) * (n - 3) * S0^2)
        }
        out$var <- V; out$z <- (1 - C) / sqrt(V)          # sign: positive z = positive autocorrelation, as for I
        out$p <- sp_test_p(out$z, if (alt == "greater") "greater" else if (alt == "less") "less" else "two.sided")
    } else {
        sims <- vapply(seq_len(nperm), function(b) cfun(x[sample.int(n)]), 0)
        out$var <- stats::var(sims); out$z <- (mean(sims) - C) / stats::sd(sims)
        out$p <- sp_perm_p(C, sims, switch(alt, greater = "less", less = "greater", alt))
    }
    out
}

# ---------------------------------------------------- Getis-Ord general G ------

sp_globalG <- function(x, W, inference = "normality", alt = "two.sided", nperm = 999, S = sp_S(W)) {
    n <- length(x); S0 <- S$S0; S1 <- S$S1; S2 <- S$S2
    if (any(x < 0)) stop("Getis-Ord G needs a variable without negative values.", call. = FALSE)
    Gf <- function(v) { den <- sum(v)^2 - sum(v^2); if (den <= 0) NA_real_ else sum(v * as.numeric(W %*% v)) / den - 0 }
    # numerator sum_{i != j} w_ij x_i x_j; the diagonal of W is zero, so the quadratic form is exactly that
    G <- Gf(x); EG <- S0 / (n * (n - 1))
    out <- list(stat = G, expect = EG)
    if (inference != "permutation") {
        n1 <- n - 1
        B0 <- (n^2 - 3 * n + 3) * S1 - n * S2 + 3 * S0^2
        B1 <- -((n^2 - n) * S1 - 2 * n * S2 + 6 * S0^2)
        B2 <- -(2 * n * S1 - (n + 3) * S2 + 6 * S0^2)
        B3 <- 4 * (n - 1) * S1 - 2 * (n + 1) * S2 + 8 * S0^2
        B4 <- S1 - S2 + S0^2
        s1 <- sum(x); s2 <- sum(x^2); s3 <- sum(x^3); s4 <- sum(x^4)
        EG2 <- (B0 * s2^2 + B1 * s4 + B2 * s1^2 * s2 + B3 * s1 * s3 + B4 * s1^4) / ((s1^2 - s2)^2 * n * (n - 1) * (n - 2) * (n - 3))
        V <- EG2 - EG^2
        out$var <- V; out$z <- (G - EG) / sqrt(V); out$p <- sp_test_p(out$z, alt)
    } else {
        sims <- vapply(seq_len(nperm), function(b) Gf(x[sample.int(n)]), 0)
        out$var <- stats::var(sims); out$z <- (G - mean(sims)) / stats::sd(sims); out$p <- sp_perm_p(G, sims, alt)
    }
    out
}

# ------------------------------------------------------------- local statistics --

# Local Moran's I_i = z_i * sum_j w_ij z_j / m2 with m2 = sum(z^2)/n, analytical moments (randomisation, Anselin 1995)
# and conditional permutation p-values (z_i fixed, the other values shuffled).
sp_local_moran <- function(x, W, nperm = 999, alt = "two.sided") {
    n <- length(x)
    z <- x - mean(x); m2 <- sum(z^2) / n
    if (m2 == 0) stop("The variable has no variation.", call. = FALSE)
    Wz <- as.numeric(W %*% z)
    Ii <- z * Wz / m2
    wi <- Matrix::rowSums(W); wi2 <- Matrix::rowSums(W^2)
    b2 <- sum(z^4) / n / m2^2
    E <- -wi / (n - 1)
    V <- wi2 * (n - b2) / (n - 1) + (wi^2 - wi2) * (2 * b2 - n) / ((n - 1) * (n - 2)) - wi^2 / (n - 1)^2
    zscore <- (Ii - E) / sqrt(V)
    # conditional permutation through the neighbour lists
    Wc <- as(W, "CsparseMatrix")
    nb <- sp_neighbours(Wc)
    wts <- lapply(seq_len(n), function(i) { j <- nb[[i]]; if (length(j)) as.numeric(Wc[i, j]) else numeric(0) })
    pp <- numeric(n); Wzs <- vector("list", n)
    for (i in seq_len(n)) {
        j <- nb[[i]]; k <- length(j)
        if (k == 0) { pp[i] <- NA; next }
        pool <- z[-i]
        sims <- numeric(nperm)
        for (b in seq_len(nperm)) sims[b] <- sum(wts[[i]] * pool[sample.int(n - 1, k)])
        simI <- z[i] * sims / m2
        pp[i] <- sp_perm_p(Ii[i], simI, alt)
    }
    list(Ii = Ii, E = E, V = V, z = zscore, p = pp, Wz = Wz, zs = z, m2 = m2)
}

# Getis-Ord Gi / Gi*: z-score of the local sum, analytical; permutation p-values conditional on x_i.
sp_local_g <- function(x, B, style = "W", star = TRUE, nperm = 999, alt = "two.sided") {
    n <- length(x)
    Bs <- if (star) B + Matrix::Diagonal(n, ifelse(Matrix::diag(B) == 0, if (style == "B") 1 else 1, 0)) else B
    Ws <- sp_style(Bs, style)
    xbar <- mean(x); s <- sqrt(sum(x^2) / n - xbar^2)
    if (s == 0) stop("The variable has no variation.", call. = FALSE)
    lag <- as.numeric(Ws %*% x)
    wi <- Matrix::rowSums(Ws); si <- Matrix::rowSums(Ws^2)
    if (star) {
        Ei <- xbar * wi
        Vi <- s^2 * (n * si - wi^2) / (n - 1)
    } else {
        # Gi without the own value: moments of the other n-1 values
        xb <- (sum(x) - x) / (n - 1)
        s2 <- (sum(x^2) - x^2) / (n - 1) - xb^2
        Ei <- xb * wi
        Vi <- s2 * (((n - 1) * si - wi^2) / (n - 2))
        lag <- lag
    }
    zv <- (lag - Ei) / sqrt(Vi)
    Wc <- as(Ws, "CsparseMatrix"); nb <- sp_neighbours(Wc)
    pp <- numeric(n)
    for (i in seq_len(n)) {
        j <- nb[[i]]; k <- length(j)
        if (k == 0) { pp[i] <- NA; next }
        w <- as.numeric(Wc[i, j])
        if (star) {
            # the self term is part of the neighbour list; keep x_i in place, shuffle the rest
            jself <- which(j == i); jo <- j[j != i]; wo <- w[j != i]; wself <- if (length(jself)) w[jself] else 0
            ko <- length(jo); pool <- x[-i]
            sims <- vapply(seq_len(nperm), function(b) wself * x[i] + sum(wo * pool[sample.int(n - 1, ko)]), 0)
        } else {
            pool <- x[-i]
            sims <- vapply(seq_len(nperm), function(b) sum(w * pool[sample.int(n - 1, k)]), 0)
        }
        pp[i] <- sp_perm_p(lag[i], sims, alt)
    }
    list(G = zv, lag = lag, p = pp, E = Ei)
}

# Quadrant / cluster class of a local statistic. cluster = "HH","LL","HL","LH" among significant units, else "Not significant".
sp_lisa_class <- function(z, Wz, p, alpha) {
    q <- ifelse(z >= 0 & Wz >= 0, "High-High", ifelse(z < 0 & Wz < 0, "Low-Low", ifelse(z >= 0, "High-Low", "Low-High")))
    q[is.na(p) | p > alpha] <- "Not significant"
    q
}

sp_adjust_p <- function(p, how) {
    if (how == "none") return(p)
    ok <- !is.na(p); out <- p
    out[ok] <- stats::p.adjust(p[ok], method = if (how == "fdr") "BH" else "bonferroni")
    out
}
