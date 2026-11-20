# Spatial regression by maximum likelihood (Ord, 1975; Anselin, 1988). Base R and Matrix only.
#
# General form:  y = rho W y + X b + u,  u = lambda W u + e,  e ~ N(0, s2 I)
#   ols  rho = 0, lambda = 0          slx   OLS on [X, WX]
#   sar  rho free                     sdm   rho free, X replaced by [X, WX]
#   sem  lambda free                  sdem  lambda free, X replaced by [X, WX]
#   sac  rho and lambda free (same W)
# The concentrated log-likelihood is maximised over (rho, lambda); standard errors come from the Hessian of the full
# log-likelihood in (rho, lambda, b, s2) using the analytic gradient.

# ------------------------------------------------------------ log-determinant --

# Returns list(ld = function(r), dld = function(r), lo, hi, method, ev)
sp_logdet_engine <- function(W, rowstd, eigen_max = 1500) {
    n <- nrow(W)
    if (n <= eigen_max) {
        ev <- eigen(as.matrix(W), only.values = TRUE)$values
        real <- Re(ev[abs(Im(ev)) < 1e-8])
        lo <- if (any(real < 0)) 1 / min(real) else -1
        hi <- if (any(real > 0)) 1 / max(real) else 1
        eps <- 1e-4
        list(ld = function(r) sum(log(Mod(1 - r * ev))),
             dld = function(r) sum(Re(-ev / (1 - r * ev))),
             lo = lo + eps * abs(lo), hi = hi - eps * abs(hi), method = "eigenvalues", ev = ev, rs = rowstd)
    } else {
        if (!rowstd) stop("With more than 1500 units the spatial regression needs row-standardised weights (style W).", call. = FALSE)
        lo <- -0.99; hi <- 0.99
        gr <- seq(lo, hi, length.out = 81)
        I <- Matrix::Diagonal(n)
        ld <- vapply(gr, function(r) as.numeric(Matrix::determinant(I - r * W, logarithm = TRUE)$modulus), 0)
        f <- stats::splinefun(gr, ld, method = "natural")
        list(ld = function(r) f(r), dld = function(r) f(r, deriv = 1), lo = lo, hi = hi, method = "sparse LU grid + spline", ev = NULL, rs = rowstd)
    }
}

# ------------------------------------------------------------- OLS diagnostics --

sp_lm_tests <- function(y, X, W) {
    n <- length(y)
    fit <- stats::lm.fit(X, y)
    u <- fit$residuals; b <- fit$coefficients
    s2 <- sum(u^2) / n
    Wy <- as.numeric(W %*% y); Wu <- as.numeric(W %*% u)
    Tr <- sum(W * Matrix::t(W)) + sum(W^2)
    WXb <- as.numeric(W %*% (X %*% b))
    # M WXb : residuals of WXb regressed on X
    MW <- stats::lm.fit(X, WXb)$residuals
    D <- sum(WXb * MW) / s2 + Tr              # = n * J
    a <- sum(u * Wu) / s2; l <- sum(u * Wy) / s2
    LMerr <- a^2 / Tr
    LMlag <- l^2 / D
    RLMerr <- (a - Tr / D * l)^2 / (Tr - Tr^2 / D)
    RLMlag <- (l - a)^2 / (D - Tr)
    SARMA <- LMerr + RLMlag
    data.frame(test = c("LM error", "LM lag", "Robust LM error", "Robust LM lag", "SARMA (LM)"),
               stat = c(LMerr, LMlag, RLMerr, RLMlag, SARMA), df = c(1, 1, 1, 1, 2),
               p = stats::pchisq(c(LMerr, LMlag, RLMerr, RLMlag, SARMA), c(1, 1, 1, 1, 2), lower.tail = FALSE),
               stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------- the fit ----

sp_sreg <- function(y, X, W, model = "sar", lagcols = NULL, rowstd = TRUE, eigen_max = 1500, ldeng = NULL, Wy = NULL, diag = TRUE) {
    n <- length(y)
    if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
    pnames <- colnames(X)
    diag_plain <- if (diag) sp_lm_tests(y, X, W) else NULL          # the LM tests always refer to the plain OLS regression
    hasW <- model %in% c("slx", "sdm", "sdem")
    if (hasW) {
        if (length(lagcols) == 0) stop("Select at least one numeric covariate to lag (a spatially lagged covariate WX).", call. = FALSE)
        WX <- as.matrix(W %*% X[, lagcols, drop = FALSE]); colnames(WX) <- paste0("W:", pnames[lagcols])
        X <- cbind(X, WX); pnames <- colnames(X)
    }
    if (qr(X)$rank < ncol(X)) stop("The regressors are collinear (including the spatial lags); remove a covariate or change which covariates are lagged.", call. = FALSE)
    p <- ncol(X)
    ols <- stats::lm.fit(X, y)
    s2_ols <- sum(ols$residuals^2) / n
    ll_ols <- -n / 2 * (log(2 * pi * s2_ols) + 1)
    k_ols <- p + 1
    out <- list(n = n, pnames = pnames, model = model, X = X, ols = ols, y = y)
    out$diag <- diag_plain
    if (model %in% c("ols", "slx")) {
        XtXi <- solve(crossprod(X)); beta <- as.numeric(ols$coefficients)
        V <- XtXi * s2_ols * n / (n - p)        # usual OLS covariance (df-corrected)
        out$coef <- data.frame(term = pnames, est = beta, se = sqrt(diag(V)), stringsAsFactors = FALSE)
        out$vcov <- V; out$beta <- beta; out$rho <- NA_real_; out$lambda <- NA_real_
        out$fitted <- as.numeric(X %*% beta); out$resid <- y - out$fitted
        out$ll <- ll_ols; out$k <- k_ols; out$s2 <- s2_ols; out$dfres <- n - p; out$ldm <- NULL
        out$sigma <- sqrt(sum(ols$residuals^2) / (n - p))
        out$coef$z <- out$coef$est / out$coef$se
        out$coef$p <- 2 * stats::pt(-abs(out$coef$z), n - p)
        out$stat_name <- "t"
        out$ll_ols <- ll_ols
        return(out)
    }
    # --- maximum likelihood
    le <- if (is.null(ldeng)) sp_logdet_engine(W, rowstd, eigen_max) else ldeng
    out$bounds <- c(le$lo, le$hi); out$ldmethod <- le$method; out$ldeng <- le
    if (is.null(Wy)) Wy <- as.numeric(W %*% y)          # panel models pass the transformed spatial lag of y
    WX <- as.matrix(W %*% X); WWy <- as.numeric(W %*% Wy)
    useRho <- model %in% c("sar", "sdm", "sac"); useLam <- model %in% c("sem", "sdem", "sac")
    conc <- function(rho, lam) {
        yt <- y - rho * Wy
        if (lam != 0) { yt <- yt - lam * (Wy - rho * WWy); Xt <- X - lam * WX } else Xt <- X
        q <- qr(Xt)
        beta <- qr.coef(q, yt); e <- qr.resid(q, yt)
        s2 <- sum(e^2) / n
        ll <- (if (useRho) le$ld(rho) else 0) + (if (useLam) le$ld(lam) else 0) - n / 2 * (log(2 * pi * s2) + 1)
        list(ll = ll, beta = as.numeric(beta), s2 = s2, e = as.numeric(e))
    }
    lo <- le$lo; hi <- le$hi
    if (model %in% c("sar", "sdm")) {
        o <- stats::optimize(function(r) conc(r, 0)$ll, c(lo, hi), maximum = TRUE, tol = 1e-9)
        rho <- o$maximum; lam <- 0
    } else if (model %in% c("sem", "sdem")) {
        o <- stats::optimize(function(l) conc(0, l)$ll, c(lo, hi), maximum = TRUE, tol = 1e-9)
        lam <- o$maximum; rho <- 0
    } else {
        st1 <- stats::optimize(function(r) conc(r, 0)$ll, c(lo, hi), maximum = TRUE)$maximum
        st2 <- stats::optimize(function(l) conc(0, l)$ll, c(lo, hi), maximum = TRUE)$maximum
        best <- NULL
        for (st in list(c(st1, 0), c(0, st2), c(0.1, 0.1))) {
            r <- tryCatch(stats::optim(st, function(v) -conc(v[1], v[2])$ll, method = "L-BFGS-B", lower = rep(lo, 2), upper = rep(hi, 2),
                                       control = list(factr = 1e3, pgtol = 0)), error = function(e) NULL)
            if (!is.null(r) && (is.null(best) || r$value < best$value)) best <- r
        }
        if (is.null(best)) stop("The spatial lag + error model did not converge.", call. = FALSE)
        rho <- best$par[1]; lam <- best$par[2]
    }
    cc <- conc(rho, lam)
    beta <- cc$beta; s2 <- cc$s2
    out$rho <- if (useRho) rho else NA_real_; out$lambda <- if (useLam) lam else NA_real_
    out$ll <- cc$ll; out$s2 <- s2; out$beta <- beta
    out$hit_bound <- (useRho && (rho - lo < 1e-3 * abs(lo) + 1e-6 || hi - rho < 1e-3 * abs(hi) + 1e-6)) ||
                     (useLam && (lam - lo < 1e-3 * abs(lo) + 1e-6 || hi - lam < 1e-3 * abs(hi) + 1e-6))
    # --- covariance from the Hessian of the full log-likelihood
    nm <- c(if (useRho) "rho", if (useLam) "lambda", pnames, "sigma2")
    ir <- if (useRho) 1 else NULL; il <- if (useLam) 1 + useRho else NULL
    nr <- useRho + useLam
    unpack <- function(ph) list(rho = if (useRho) ph[1] else 0, lam = if (useLam) ph[1 + useRho] else 0, b = ph[nr + seq_len(p)], s2 = ph[nr + p + 1])
    fullpieces <- function(ph) {
        u <- unpack(ph)
        Ay <- y - u$rho * Wy; r <- Ay - as.numeric(X %*% u$b)
        Br <- if (u$lam != 0) r - u$lam * as.numeric(W %*% r) else r
        list(u = u, r = r, e = Br)
    }
    negll <- function(ph) {
        f <- fullpieces(ph); if (f$u$s2 <= 0) return(1e300)
        -((if (useRho) le$ld(f$u$rho) else 0) + (if (useLam) le$ld(f$u$lam) else 0) - n / 2 * log(2 * pi * f$u$s2) - sum(f$e^2) / (2 * f$u$s2))
    }
    neggr <- function(ph) {
        f <- fullpieces(ph); u <- f$u; e <- f$e
        Xt <- if (u$lam != 0) X - u$lam * WX else X
        g <- numeric(length(ph))
        if (useRho) {
            BWy <- if (u$lam != 0) Wy - u$lam * as.numeric(W %*% Wy) else Wy
            g[1] <- le$dld(u$rho) + sum(e * BWy) / u$s2
        }
        if (useLam) g[1 + useRho] <- le$dld(u$lam) + sum(e * as.numeric(W %*% f$r)) / u$s2
        g[nr + seq_len(p)] <- as.numeric(crossprod(Xt, e)) / u$s2
        g[nr + p + 1] <- -n / (2 * u$s2) + sum(e^2) / (2 * u$s2^2)
        -g
    }
    ph0 <- c(if (useRho) rho, if (useLam) lam, beta, s2)
    H <- tryCatch(stats::optimHess(ph0, negll, neggr), error = function(e) NULL)
    V <- if (!is.null(H)) tryCatch(solve(H), error = function(e) NULL) else NULL
    if (is.null(V) || any(!is.finite(diag(V))) || any(diag(V) <= 0)) {
        V <- matrix(NA_real_, length(ph0), length(ph0)); out$vcov_ok <- FALSE
    } else out$vcov_ok <- TRUE
    dimnames(V) <- list(nm, nm)
    se <- sqrt(diag(V))
    tab <- data.frame(term = nm, est = ph0, se = se, stringsAsFactors = FALSE)
    out$vcov <- V
    out$coef <- tab[tab$term != "sigma2", , drop = FALSE]
    out$coef$z <- out$coef$est / out$coef$se
    out$coef$p <- 2 * stats::pnorm(-abs(out$coef$z))
    out$stat_name <- "z"
    out$sigma <- sqrt(s2)
    out$k <- nr + p + 1
    out$fitted <- as.numeric(X %*% beta) + (if (useRho) rho * Wy else 0)
    out$resid <- y - out$fitted
    out$ll_ols <- ll_ols; out$dfres <- n - p
    # likelihood ratio test of the spatial parameter(s) against OLS
    out$lr <- 2 * (out$ll - ll_ols); out$lr_df <- nr
    out$lr_p <- stats::pchisq(out$lr, nr, lower.tail = FALSE)
    out
}

# ----------------------------------------------------------- impact measures ---

# Traces of W^k, k = 0..K: exact from eigenvalues when available, Hutchinson estimate otherwise.
sp_traces <- function(W, ev = NULL, K = 60, probes = 40) {
    n <- nrow(W)
    if (!is.null(ev)) return(list(tr = vapply(0:K, function(k) Re(sum(ev^k)), 0), exact = TRUE))
    set.seed(1)
    V <- matrix(sample(c(-1, 1), n * probes, TRUE), n, probes)
    tr <- numeric(K + 1); P <- V
    tr[1] <- n
    for (k in seq_len(K)) { P <- as.matrix(W %*% P); tr[k + 1] <- mean(colSums(V * P)) }
    list(tr = tr, exact = FALSE)
}

# Average direct / indirect / total effects for the covariates (rows of `vars`: name, beta index, theta index or NA).
# `par` holds rho, beta (all columns of X).
sp_effects_point <- function(rho, beta, vars, W, ldeng, tr = NULL) {
    n <- nrow(W)
    I <- Matrix::Diagonal(n)
    A <- I - rho * W
    one <- rep(1, n)
    if (isTRUE(ldeng$rs)) { a <- b <- one } else {
        a <- as.numeric(Matrix::solve(A, one))
        b <- as.numeric(Matrix::solve(A, as.numeric(W %*% one)))
    }
    if (!is.null(ldeng$ev)) {
        ev <- ldeng$ev
        d0 <- Re(sum(1 / (1 - rho * ev))) / n
        d1 <- Re(sum(ev / (1 - rho * ev))) / n
    } else {
        K <- length(tr$tr) - 1
        pw <- rho^(0:K)
        d0 <- sum(pw * tr$tr) / n
        d1 <- sum(pw[1:K] * tr$tr[2:(K + 1)]) / n
    }
    m0 <- mean(a); m1 <- mean(b)
    if (isTRUE(ldeng$rs)) { m0 <- 1 / (1 - rho); m1 <- m0 }
    res <- lapply(seq_len(nrow(vars)), function(v) {
        bk <- beta[vars$b[v]]; th <- if (is.na(vars$t[v])) 0 else beta[vars$t[v]]
        dir <- bk * d0 + th * d1; tot <- bk * m0 + th * m1
        c(direct = dir, indirect = tot - dir, total = tot)
    })
    do.call(rbind, res)
}

# Simulated inference (draws from the asymptotic normal distribution of rho and beta).
sp_effects <- function(fit, vars, W, nsim = 500, seed = 1) {
    n <- nrow(W)
    tr <- if (is.null(fit$ldeng$ev)) sp_traces(W) else NULL
    est <- sp_effects_point(fit$rho, fit$beta, vars, W, fit$ldeng, tr)
    out <- list(est = est, se = NULL, sims = NULL)
    if (isTRUE(fit$vcov_ok) && nsim > 0) {
        nm <- c("rho", fit$pnames)
        mu <- c(fit$rho, fit$beta)
        V <- fit$vcov[nm, nm, drop = FALSE]
        set.seed(seed)
        ev <- eigen(V, symmetric = TRUE)
        L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), length(mu))
        draws <- 0; keep <- list()
        tries <- 0
        while (sum(vapply(keep, nrow, 1L)) < nsim && tries < 5) {
            Z <- matrix(stats::rnorm(length(mu) * nsim), nsim)
            D <- sweep(Z %*% t(L), 2, mu, "+")
            D <- D[D[, 1] > fit$bounds[1] & D[, 1] < fit$bounds[2], , drop = FALSE]
            keep[[length(keep) + 1]] <- D; tries <- tries + 1
        }
        D <- do.call(rbind, keep); D <- D[seq_len(min(nrow(D), nsim)), , drop = FALSE]
        sims <- lapply(seq_len(nrow(D)), function(s) sp_effects_point(D[s, 1], D[s, -1], vars, W, fit$ldeng, tr))
        arr <- simplify2array(sims)                      # vars x 3 x nsim
        out$se <- apply(arr, c(1, 2), stats::sd)
        out$nsim <- nrow(D)
    }
    out
}
