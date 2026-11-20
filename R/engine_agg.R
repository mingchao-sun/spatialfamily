# Industrial agglomeration and regional specialisation indices; PPML gravity estimation. Base R only.
#
# Input for the indices: a long table region x industry with an employment (or output) value and optionally the number of
# establishments. Missing region-industry combinations count as zero.

# Matrix E[region, industry] from long vectors (duplicates are summed).
sp_agg_matrix <- function(region, industry, value, firms = NULL) {
    r <- factor(as.character(region)); i <- factor(as.character(industry))
    E <- matrix(0, nlevels(r), nlevels(i), dimnames = list(levels(r), levels(i)))
    idx <- cbind(as.integer(r), as.integer(i))
    for (k in seq_along(value)) E[idx[k, 1], idx[k, 2]] <- E[idx[k, 1], idx[k, 2]] + value[k]
    Nf <- NULL
    if (!is.null(firms)) {
        Nf <- E * 0
        for (k in seq_along(value)) Nf[idx[k, 1], idx[k, 2]] <- Nf[idx[k, 1], idx[k, 2]] + firms[k]
    }
    list(E = E, N = Nf, r = r, i = i)
}

# Location quotient of every cell: (E_ri / E_r.) / (E_.i / E_..).
sp_agg_lq <- function(E) {
    tot <- sum(E); rs <- rowSums(E); cs <- colSums(E)
    LQ <- E / rs %o% rep(1, ncol(E)) / (rep(1, nrow(E)) %o% (cs / tot))
    LQ[rs == 0, ] <- NA; LQ[, cs == 0] <- NA
    LQ
}

# Industry-level indices. Returns a data frame with one row per industry.
sp_agg_industry <- function(E, N = NULL) {
    tot <- sum(E); x <- rowSums(E) / tot            # region shares of total employment
    sx2 <- sum(x^2)
    out <- lapply(seq_len(ncol(E)), function(j) {
        e <- E[, j]; ej <- sum(e)
        if (ej <= 0) return(NULL)
        s <- e / ej                                   # region shares of the industry
        G <- sum((s - x)^2)                           # raw concentration (Ellison and Glaeser, 1997)
        H <- if (is.null(N)) NA_real_ else { n <- N[, j]; ok <- n > 0; sum(e[ok]^2 / n[ok]) / ej^2 }   # plant Herfindahl, equal plants within a cell
        h <- if (is.na(H)) 0 else H
        gam <- (G - (1 - sx2) * h) / ((1 - sx2) * (1 - h))
        # locational Gini (Krugman, 1991): Lorenz curve of the industry share against the region share, regions ordered by LQ
        lq <- ifelse(x > 0, s / x, NA)
        o <- order(lq); xs <- cumsum(x[o]); ss <- cumsum(s[o])
        xs0 <- c(0, xs); ss0 <- c(0, ss)
        gini <- 1 - sum((xs0[-1] - xs0[-length(xs0)]) * (ss0[-1] + ss0[-length(ss0)]))
        c(total = ej, share = ej / tot, hhi = sum(s^2), G = G, H = H, gamma = gam, gini = gini, nreg = sum(e > 0))
    })
    nm <- colnames(E)
    keep <- !vapply(out, is.null, TRUE)
    d <- as.data.frame(do.call(rbind, out[keep])); d$industry <- nm[keep]
    d[, c("industry", "total", "share", "nreg", "hhi", "G", "H", "gamma", "gini")]
}

# Region-level indices: Herfindahl of the industry shares and the Krugman dissimilarity index against the rest of the economy.
sp_agg_region <- function(E) {
    tot <- sum(E); cs <- colSums(E)
    res <- lapply(seq_len(nrow(E)), function(r) {
        e <- E[r, ]; er <- sum(e)
        if (er <= 0) return(NULL)
        s <- e / er; rest <- (cs - e) / (tot - er)
        c(total = er, share = er / tot, hhi = sum(s^2), entropy = -sum(ifelse(s > 0, s * log(s), 0)), krugman = sum(abs(s - rest)))
    })
    keep <- !vapply(res, is.null, TRUE)
    d <- as.data.frame(do.call(rbind, res[keep])); d$region <- rownames(E)[keep]
    d[, c("region", "total", "share", "hhi", "entropy", "krugman")]
}

# ---------------------------------------------------------------------- PPML ----

# Weighted demeaning by one or more fixed effects (alternating projections). `fe` is a list of integer group vectors.
sp_absorb <- function(M, w, fe, tol = 1e-10, maxit = 1000) {
    M <- as.matrix(M)
    if (length(fe) == 0) return(M)
    for (j in seq_len(ncol(M))) {
        v <- M[, j]
        for (it in seq_len(maxit)) {
            vold <- v
            for (g in fe) {
                m <- rowsum(w * v, g, reorder = TRUE) / rowsum(w, g, reorder = TRUE)
                v <- v - m[g]
            }
            if (length(fe) == 1 || max(abs(v - vold)) < tol * (1 + max(abs(vold)))) break
        }
        M[, j] <- v
    }
    M
}

# Poisson pseudo maximum likelihood with high-dimensional fixed effects (Santos Silva and Tenreyro, 2006; Correia et al., 2020).
# X: n x p covariates (no constant when fe is used; a constant is added otherwise). cluster: list of integer vectors (multiway) or NULL.
sp_ppml <- function(y, X, fe = list(), cluster = NULL, tol = 1e-9, maxit = 100) {
    n <- length(y)
    if (any(y < 0)) stop("The dependent variable (flow) must not be negative.", call. = FALSE)
    keep <- rep(TRUE, n); dropped <- 0
    repeat {                                         # drop observations whose fixed-effect group has only zero flows (separation)
        ch <- FALSE
        for (g in fe) {
            gs <- tapply(y[keep], g[keep], sum)
            bad <- names(gs)[gs <= 0]
            if (length(bad)) { rm <- keep & (as.character(g) %in% bad); if (any(rm)) { keep <- keep & !rm; ch <- TRUE } }
        }
        if (!ch) break
    }
    dropped <- sum(!keep)
    if (sum(keep) < ncol(X) + 5) stop("Too few observations are left after removing groups with only zero flows.", call. = FALSE)
    y <- y[keep]; X <- X[keep, , drop = FALSE]; fe <- lapply(fe, function(g) match(g[keep], unique(g[keep])))
    if (!is.null(cluster)) cluster <- lapply(cluster, function(g) g[keep])
    n <- length(y)
    if (length(fe) == 0) X <- cbind(`(Intercept)` = 1, X)
    p <- ncol(X)
    mu <- (y + mean(y)) / 2; eta <- log(mu); dev_old <- Inf; conv <- FALSE
    dev <- function(mu) 2 * sum(ifelse(y > 0, y * log(y / mu), 0) - (y - mu))
    for (it in seq_len(maxit)) {
        z <- eta + (y - mu) / mu; w <- mu
        Zd <- sp_absorb(cbind(z, X), w, fe)
        zd <- Zd[, 1]; Xd <- Zd[, -1, drop = FALSE]
        zero <- colnames(X)[apply(Xd, 2, function(v) max(abs(v)) < 1e-8 * (1 + max(abs(X))))]
        if (length(zero)) stop(sprintf("Covariate(s) absorbed by the fixed effects (they do not vary within the fixed-effect groups): %s.", paste(zero, collapse = ", ")), call. = FALSE)
        XtWX <- crossprod(Xd, Xd * w)
        beta <- tryCatch(solve(XtWX, crossprod(Xd, w * zd)), error = function(e) NULL)
        if (is.null(beta)) stop("The covariates are collinear (or absorbed by the fixed effects).", call. = FALSE)
        e <- zd - Xd %*% beta
        eta <- as.numeric(z - e); eta <- pmin(eta, 700); mu <- exp(eta)
        d <- dev(mu)
        if (abs(d - dev_old) / (0.1 + abs(d)) < tol) { conv <- TRUE; break }
        dev_old <- d
    }
    # final weights and scores
    z <- eta + (y - mu) / mu; w <- mu
    Zd <- sp_absorb(cbind(z, X), w, fe); Xd <- Zd[, -1, drop = FALSE]
    bread <- solve(crossprod(Xd, Xd * w))
    score <- Xd * (y - mu)
    meat1 <- function(cl) {
        if (is.null(cl)) return(crossprod(score))
        S <- rowsum(score, cl)
        crossprod(S) * (nlevels(factor(cl)) / (nlevels(factor(cl)) - 1))
    }
    if (is.null(cluster)) V <- bread %*% (crossprod(score) * n / (n - 1)) %*% bread
    else {
        # multiway clustering (Cameron, Gelbach and Miller, 2011): add the one-way matrices, subtract the intersections
        k <- length(cluster); M <- 0
        for (m in seq_len(2^k - 1)) {
            idx <- which(bitwAnd(m, 2^(seq_len(k) - 1)) > 0)
            cl <- if (length(idx) == 1) cluster[[idx]] else as.integer(interaction(cluster[idx], drop = TRUE))
            M <- M + (-1)^(length(idx) + 1) * meat1(cl)
        }
        V <- bread %*% M %*% bread
    }
    se <- sqrt(pmax(diag(V), 0))
    ll <- sum(ifelse(y > 0, y * log(mu), 0) - mu) - sum(lgamma(y + 1))
    list(coef = data.frame(term = colnames(X), est = as.numeric(beta), se = se, stringsAsFactors = FALSE), vcov = V,
         mu = mu, y = y, n = n, dropped = dropped, keep = keep, converged = conv, iterations = it, deviance = dev(mu), ll = ll,
         pr2 = stats::cor(y, mu)^2, nfe = vapply(fe, function(g) length(unique(g)), 1L), nclusters = if (is.null(cluster)) NULL else vapply(cluster, function(g) length(unique(g)), 1L))
}
