# Difference-in-differences with spatial spillovers: two-way fixed effects, distance rings around treated units, spatially lagged
# treatment, and standard errors clustered by unit, robust, or Conley spatial HAC (Conley, 1999; Hsiang, 2010). Base R and Matrix.

# Distance from every unit to the nearest treated unit in each period. `D` is the N x N distance matrix, `treat` an N x T 0/1 matrix
# (units in the rows of D, periods in the columns, NA = unit not observed). Returns an N x T matrix (Inf when nobody is treated).
sp_did_nearest <- function(D, treat) {
    N <- nrow(D); T <- ncol(treat)
    out <- matrix(Inf, N, T)
    for (t in seq_len(T)) {
        tr <- which(treat[, t] == 1)
        if (length(tr) == 0) next
        out[, t] <- if (length(tr) == 1) D[, tr] else apply(D[, tr, drop = FALSE], 1, min)
        # a treated unit is not its own spillover source; its own row is ignored later (ring dummies apply to untreated units only)
    }
    out
}

# Ring index of every unit-period (0 = no ring: treated, or beyond the last cut-off). cuts: increasing distances.
sp_did_rings <- function(dmin, treat, cuts) {
    ring <- matrix(0L, nrow(dmin), ncol(dmin))
    for (k in seq_along(cuts)) {
        lo <- if (k == 1) 0 else cuts[k - 1]
        ring[dmin >= lo & dmin < cuts[k] & !is.na(treat) & treat == 0] <- k
    }
    ring
}

# Conley spatial HAC meat: sum over periods t, s with |t - s| <= L of Psi_t' K Psi_s, Bartlett kernels in distance and time.
# psi: list over periods of (N_t x p) score matrices; unit: list of unit indices (rows of D) for each period.
sp_conley_meat <- function(psi, unit, D, bw, lag) {
    T <- length(psi); p <- ncol(psi[[1]])
    meat <- matrix(0, p, p)
    K <- 1 - D / bw; K[K < 0] <- 0
    for (t in seq_len(T)) {
        for (s in max(1, t - lag):min(T, t + lag)) {
            wt <- 1 - abs(t - s) / (lag + 1)
            Kts <- K[unit[[t]], unit[[s]], drop = FALSE]
            meat <- meat + wt * crossprod(psi[[t]], Kts %*% psi[[s]])
        }
    }
    (meat + t(meat)) / 2
}

# Two-way fixed effects regression. y: n; X: n x p (no constant); unit, time: integer ids; `D` N x N distances (only for Conley).
# se: "cluster" (by unit), "robust", "conley" (bandwidth `bw`, serial lag `lag`).
sp_did_fit <- function(y, X, unit, time, se = "cluster", D = NULL, bw = NULL, lag = 0, effects = "twoway") {
    n <- length(y); p <- ncol(X)
    fe <- switch(effects, twoway = list(unit, time), unit = list(unit), time = list(time))
    w <- rep(1, n)
    Z <- sp_absorb(cbind(y, X), w, fe)
    yd <- Z[, 1]; Xd <- Z[, -1, drop = FALSE]
    zero <- colnames(X)[apply(Xd, 2, function(v) max(abs(v)) < 1e-9 * (1 + max(abs(X))))]
    if (length(zero)) stop(sprintf("Regressor(s) without variation after removing the fixed effects: %s.", paste(zero, collapse = ", ")), call. = FALSE)
    XtX <- crossprod(Xd)
    beta <- tryCatch(solve(XtX, crossprod(Xd, yd)), error = function(e) NULL)
    if (is.null(beta)) stop("The regressors are collinear after removing the fixed effects.", call. = FALSE)
    e <- as.numeric(yd - Xd %*% beta)
    bread <- solve(XtX)
    nunit <- length(unique(unit)); ntime <- length(unique(time))
    Kfe <- switch(effects, twoway = nunit + ntime - 1, unit = nunit, time = ntime)
    score <- Xd * e
    if (se == "cluster") {
        S <- rowsum(score, unit); G <- nrow(S)
        # unit effects are nested in the clusters and are not counted in K (fixest and Stata convention)
        K <- p + switch(effects, unit = 1, twoway = ntime, time = ntime)
        adj <- G / (G - 1) * (n - 1) / (n - K)
        V <- bread %*% (crossprod(S) * adj) %*% bread
        dfres <- G - 1
    } else if (se == "robust") {
        adj <- n / (n - p - Kfe)
        V <- bread %*% (crossprod(score) * adj) %*% bread
        dfres <- n - p - Kfe
    } else {
        if (is.null(D) || is.null(bw)) stop("Conley standard errors need coordinates and a distance bandwidth.", call. = FALSE)
        tl <- sort(unique(time))
        psi <- lapply(tl, function(t) score[time == t, , drop = FALSE]); uu <- lapply(tl, function(t) unit[time == t])
        meat <- sp_conley_meat(psi, uu, D, bw, lag)
        V <- bread %*% meat %*% bread
        dfres <- nunit - 1
    }
    se_v <- sqrt(pmax(diag(V), 0))
    tss <- sum(yd^2)
    list(coef = data.frame(term = colnames(X), est = as.numeric(beta), se = se_v, stringsAsFactors = FALSE), vcov = V, resid = e, r2w = 1 - sum(e^2) / tss,
         n = n, nunit = nunit, ntime = ntime, dfres = dfres)
}
