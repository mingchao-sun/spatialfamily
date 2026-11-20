# Spatial panel regression by maximum likelihood (Elhorst, 2003, 2014): fixed effects by the within transformation,
# the spatial lag of y and of the covariates computed before the transformation, log-determinant T * ln|I - rho W|.
# Balanced panel stacked time by time (units in the order of the rows and columns of W).
#   y: NT vector; X: NT x p matrix of covariates WITHOUT a constant; W: N x N sparse weights.
#   effects: "none" (pooled, with a constant), "individual", "time", "twoway".
#   model: "ols", "slx", "sar", "sem", "sdm", "sdem". `lagcols`: columns of X that get a spatial lag.

sp_panel_transform <- function(v, N, T, effects) {
    m <- matrix(v, N, T)
    if (effects == "individual") m <- m - rowMeans(m)
    else if (effects == "time") m <- sweep(m, 2, colMeans(m))
    else if (effects == "twoway") m <- m - rowMeans(m) - matrix(colMeans(m), N, T, byrow = TRUE) + mean(m)
    as.numeric(m)
}

sp_panel_lag <- function(v, W, N, T) as.numeric(as.matrix(W %*% matrix(v, N, T)))

sp_panel_fit <- function(y, X, W, N, T, model = "sar", effects = "individual", lagcols = NULL, rowstd = TRUE, eigen_max = 1500) {
    stopifnot(length(y) == N * T, nrow(X) == N * T)
    if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
    hasWX <- model %in% c("slx", "sdm", "sdem")
    Xa <- X
    if (hasWX) {
        if (length(lagcols) == 0) stop("Select at least one numeric covariate to lag (a spatially lagged covariate WX).", call. = FALSE)
        WX <- vapply(lagcols, function(j) sp_panel_lag(X[, j], W, N, T), numeric(N * T))
        WX <- matrix(WX, N * T); colnames(WX) <- paste0("W:", colnames(X)[lagcols])
        Xa <- cbind(X, WX)
    }
    if (effects == "none") { Xa <- cbind(`(Intercept)` = 1, Xa); lagshift <- 1 } else lagshift <- 0
    ys <- sp_panel_transform(y, N, T, effects)
    Xs <- apply(Xa, 2, sp_panel_transform, N = N, T = T, effects = effects)
    Xs <- matrix(Xs, N * T, ncol(Xa), dimnames = list(NULL, colnames(Xa)))
    keep <- apply(Xs, 2, function(v) any(abs(v) > 1e-10))
    if (!all(keep)) {
        bad <- colnames(Xs)[!keep]
        stop(sprintf("Covariate(s) without variation after removing the fixed effects: %s. Remove them (a time-invariant covariate cannot be estimated with unit fixed effects).", paste(bad, collapse = ", ")), call. = FALSE)
    }
    if (qr(Xs)$rank < ncol(Xs)) stop("The regressors are collinear after removing the fixed effects.", call. = FALSE)
    out <- list(N = N, T = T, effects = effects, model = model, pnames = colnames(Xs))
    SST <- sum((ys - mean(ys))^2)
    if (model %in% c("ols", "slx")) {
        fit <- stats::lm.fit(Xs, ys)
        n <- N * T; p <- ncol(Xs); e <- fit$residuals; s2 <- sum(e^2) / n
        # degrees of freedom lost to the fixed effects
        dfe <- switch(effects, none = 0, individual = N, time = T, twoway = N + T - 1)
        V <- solve(crossprod(Xs)) * sum(e^2) / (n - p - dfe)
        beta <- as.numeric(fit$coefficients)
        tab <- data.frame(term = colnames(Xs), est = beta, se = sqrt(diag(V)), stringsAsFactors = FALSE)
        tab$z <- tab$est / tab$se; tab$p <- 2 * stats::pt(-abs(tab$z), n - p - dfe)
        out <- c(out, list(coef = tab, beta = beta, vcov = V, rho = NA_real_, lambda = NA_real_, ll = -n / 2 * (log(2 * pi * s2) + 1), k = p + 1,
                           sigma = sqrt(s2), resid = e, fitted = ys - e, stat_name = "t", vcov_ok = TRUE, dfres = n - p - dfe))
        out$r2w <- 1 - sum(e^2) / SST
        return(out)
    }
    Wbig <- Matrix::kronecker(Matrix::Diagonal(T), W)
    le <- sp_logdet_engine(W, rowstd, eigen_max)
    lew <- list(ld = function(r) T * le$ld(r), dld = function(r) T * le$dld(r), lo = le$lo, hi = le$hi, method = le$method, ev = le$ev, rs = le$rs)
    core <- if (model %in% c("sar", "sdm")) "sar" else "sem"
    Wys <- if (core == "sar") sp_panel_transform(sp_panel_lag(y, W, N, T), N, T, effects) else NULL
    f <- sp_sreg(ys, Xs, Wbig, core, NULL, rowstd, ldeng = lew, Wy = Wys, diag = FALSE)
    f$N <- N; f$T <- T; f$effects <- effects; f$model <- model; f$W <- W; f$ldeng <- lew
    f$r2w <- 1 - sum(f$resid^2) / SST
    f$lagshift <- lagshift
    # degrees-of-freedom correction for the fixed effects removed by the within transformation (as in Stata's spxtregress)
    dfe <- switch(effects, none = 0, individual = N, time = T, twoway = N + T - 1)
    if (dfe > 0 && isTRUE(f$vcov_ok)) {
        fac <- N * T / (N * T - dfe)
        f$vcov <- f$vcov * fac
        se <- sqrt(diag(f$vcov)); names(se) <- rownames(f$vcov)
        f$coef$se <- se[f$coef$term]; f$coef$z <- f$coef$est / f$coef$se; f$coef$p <- 2 * stats::pnorm(-abs(f$coef$z))
        f$sigma <- f$sigma * sqrt(fac); f$dffac <- fac
    }
    f
}

# Effects for the panel models: same formulas as the cross-section on the N x N weights.
sp_panel_effects_vars <- function(fit, p0, lagcols, shift = 0) {
    cols <- (1 + shift):p0
    tpos <- match(cols - shift, lagcols)
    data.frame(b = cols, t = ifelse(is.na(tpos), NA, p0 + tpos))
}
