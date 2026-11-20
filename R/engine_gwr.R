# Geographically weighted regression (Brunsdon, Fotheringham and Charlton, 1996; Fotheringham et al., 2002). Base R only.
# Conventions follow the GWmodel package: adaptive bandwidth = number of nearest neighbours (the unit itself counts as the first),
# AICc = 2 n ln(sqrt(RSS / n)) + n ln(2 pi) + n (n + tr(S)) / (n - 2 - tr(S)),
# standard errors use sigma^2 = RSS / (n - 2 tr(S) + tr(S'S)); the effective number of parameters is 2 tr(S) - tr(S'S).
# Local R2 = 1 - weighted RSS / weighted TSS around the weighted mean (GWmodel's Local_R2 is defined differently).

sp_dist_matrix <- function(xy, lonlat) {
    n <- nrow(xy)
    if (n > 4000) stop("Geographically weighted regression is limited to 4000 units (it needs the full distance matrix).", call. = FALSE)
    if (!lonlat) return(as.matrix(stats::dist(xy)))
    D <- matrix(0, n, n)
    for (i in seq_len(n)) D[i, ] <- sp_dist_from(xy, i, TRUE)
    D
}

sp_gwr_weights <- function(d, b, kernel) {
    switch(kernel,
           gaussian = exp(-0.5 * (d / b)^2),
           exponential = exp(-d / b),
           bisquare = ifelse(d < b, (1 - (d / b)^2)^2, 0),
           tricube = ifelse(d < b, (1 - (d / b)^3)^3, 0),
           boxcar = ifelse(d <= b, 1, 0),
           stop("Unknown kernel.", call. = FALSE))
}

# Bandwidth of unit i: adaptive = distance to the N-th nearest unit (itself first), fixed = bw.
sp_gwr_bw_i <- function(d, bw, adaptive) {
    if (!adaptive) return(bw)
    k <- min(length(d), max(1, round(bw)))
    sort(d, partial = k)[k]
}

# One pass: coefficients, hat diagonal, fitted values. `full` also returns standard-error ingredients, tr(S'S) and local R2.
sp_gwr_pass <- function(y, X, D, bw, adaptive, kernel, full = FALSE) {
    n <- length(y); p <- ncol(X)
    beta <- matrix(NA_real_, n, p); sii <- numeric(n); yhat <- numeric(n)
    cc <- if (full) matrix(NA_real_, n, p) else NULL
    stS <- numeric(n); r2 <- numeric(n)
    for (i in seq_len(n)) {
        di <- D[i, ]
        b <- sp_gwr_bw_i(di, bw, adaptive)
        w <- sp_gwr_weights(di, b, kernel)
        if (sum(w > 0) <= p) return(NULL)
        XtW <- t(X * w)
        A <- XtW %*% X
        C <- tryCatch(solve(A, XtW), error = function(e) NULL)
        if (is.null(C)) return(NULL)
        bi <- as.numeric(C %*% y)
        beta[i, ] <- bi; yhat[i] <- sum(X[i, ] * bi); sii[i] <- sum(X[i, ] * C[, i])
        if (full) {
            Si <- as.numeric(X[i, ] %*% C)
            stS[i] <- sum(Si^2)
            cc[i, ] <- rowSums(C^2)
            fl <- as.numeric(X %*% bi); sw <- sum(w); ybar <- sum(w * y) / sw
            tss <- sum(w * (y - ybar)^2); r2[i] <- if (tss > 0) 1 - sum(w * (y - fl)^2) / tss else NA_real_
        }
    }
    list(beta = beta, sii = sii, yhat = yhat, cc = cc, stS = stS, r2 = r2)
}

sp_gwr_aicc <- function(y, yhat, trS) {
    n <- length(y); rss <- sum((y - yhat)^2)
    if (n - 2 - trS <= 0) return(Inf)
    2 * n * log(sqrt(rss / n)) + n * log(2 * pi) + n * (n + trS) / (n - 2 - trS)
}

sp_gwr_cv <- function(y, X, D, bw, adaptive, kernel) {
    n <- length(y); p <- ncol(X); cv <- 0
    for (i in seq_len(n)) {
        di <- D[i, ]; w <- sp_gwr_weights(di, sp_gwr_bw_i(di, bw, adaptive), kernel); w[i] <- 0
        if (sum(w > 0) <= p) return(Inf)
        XtW <- t(X * w); bi <- tryCatch(solve(XtW %*% X, XtW %*% y), error = function(e) NULL)
        if (is.null(bi)) return(Inf)
        cv <- cv + (y[i] - sum(X[i, ] * bi))^2
    }
    cv
}

# Bandwidth selection by golden-section search on the AICc (or the cross-validation score).
sp_gwr_select <- function(y, X, D, adaptive, kernel, criterion = "aicc") {
    n <- length(y); p <- ncol(X)
    score <- function(bw) {
        if (criterion == "cv") return(sp_gwr_cv(y, X, D, bw, adaptive, kernel))
        r <- sp_gwr_pass(y, X, D, bw, adaptive, kernel)
        if (is.null(r)) Inf else sp_gwr_aicc(y, r$yhat, sum(r$sii))
    }
    if (adaptive) { lo <- max(p + 2, 5); hi <- n } else {
        dmax <- max(D); dmin <- stats::quantile(D[D > 0], 0.01); lo <- as.numeric(dmin); hi <- dmax
    }
    o <- stats::optimize(function(b) { s <- score(b); if (is.finite(s)) s else 1e300 }, c(lo, hi), tol = if (adaptive) 0.5 else (hi - lo) * 1e-4)
    bw <- o$minimum
    if (adaptive) {
        cand <- unique(pmin(hi, pmax(lo, c(floor(bw), ceiling(bw)))))
        sc <- vapply(cand, score, 0); bw <- cand[which.min(sc)]
    }
    list(bw = bw, score = score(bw))
}

sp_gwr <- function(y, X, xy, lonlat = FALSE, kernel = "bisquare", adaptive = TRUE, bw = NULL, criterion = "aicc") {
    n <- length(y); p <- ncol(X)
    D <- sp_dist_matrix(xy, lonlat)
    if (is.null(bw)) {
        sel <- sp_gwr_select(y, X, D, adaptive, kernel, criterion)
        bw <- sel$bw
    } else sel <- NULL
    r <- sp_gwr_pass(y, X, D, bw, adaptive, kernel, full = TRUE)
    if (is.null(r)) stop("The bandwidth is too small: a local regression has too few observations with positive weight. Increase the bandwidth.", call. = FALSE)
    trS <- sum(r$sii); trStS <- sum(r$stS)
    rss <- sum((y - r$yhat)^2)
    sigma2 <- rss / (n - 2 * trS + trStS)
    se <- sqrt(sigma2 * r$cc)
    tss <- sum((y - mean(y))^2)
    ols <- stats::lm.fit(X, y); rss0 <- sum(ols$residuals^2)
    out <- list(beta = r$beta, se = se, tval = r$beta / se, yhat = r$yhat, resid = y - r$yhat, localR2 = r$r2, trS = trS, trStS = trStS, rss = rss,
                enp = 2 * trS - trStS, sigma2 = sigma2, bw = bw, adaptive = adaptive, kernel = kernel, n = n, p = p, hat = r$sii,
                aicc = sp_gwr_aicc(y, r$yhat, trS), aic = 2 * n * log(sqrt(rss / n)) + n * log(2 * pi) + n + trS,
                r2 = 1 - rss / tss, r2adj = 1 - (1 - (1 - rss / tss)) * (n - 1) / (n - trS - 1),
                ols_coef = as.numeric(ols$coefficients), ols_aicc = sp_gwr_aicc(y, y - ols$residuals, p), ols_r2 = 1 - rss0 / tss,
                ols_se = sqrt(diag(solve(crossprod(X))) * rss0 / (n - p)), selection = sel, criterion = criterion)
        out
}
