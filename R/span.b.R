# jamovi analysis class: spatial panel regression (maximum likelihood, fixed effects by the within transformation).

# Balanced panel from the long data: complete cases, units observed in every period, rows ordered by period then unit.
sp_panel_data <- function(df, unit, time) {
    u <- as.character(df[[unit]]); tt <- df[[time]]
    tl <- if (is.factor(tt)) levels(droplevels(tt)) else as.character(sort(unique(tt)))
    tc <- as.character(tt)
    if (anyDuplicated(paste(u, tc))) stop("A unit appears more than once in the same period; the data need one row per unit and period.", call. = FALSE)
    T <- length(tl)
    if (T < 2) stop("The panel needs at least two time periods.", call. = FALSE)
    cnt <- table(u)
    keep_u <- names(cnt)[cnt == T]
    n_unit_drop <- length(cnt) - length(keep_u)
    if (length(keep_u) < 10) stop("Fewer than 10 units are observed in every period.", call. = FALSE)
    uu <- as.character(sort(if (is.numeric(df[[unit]])) as.numeric(keep_u) else keep_u))
    keep <- u %in% keep_u
    df <- df[keep, , drop = FALSE]; u <- u[keep]; tc <- tc[keep]
    ord <- order(match(tc, tl), match(u, uu))
    list(df = df[ord, , drop = FALSE], units = uu, times = tl, N = length(uu), T = T, n_unit_drop = n_unit_drop, n_obs_drop = sum(!keep))
}

spanClass <- R6::R6Class(
    "spanClass",
    inherit = spanBase,
    private = list(

        .run = function() {
            o <- self$options; res <- self$results
            sp_show_text(res, o, "span")
            if (is.null(o$dep) || length(c(o$covs, o$noLag)) == 0 || is.null(o$timeVar) || is.null(o$idVar) || !sp_ready(o)) {
                res$notes$setContent(sp_ready_msg(o, "Select a dependent variable, covariates, the unit id and the time period, and set up the spatial weights. "))
                return()
            }
            model <- o$model; eff <- o$panelEffects
            if (model %in% c("slx", "sdm", "sdem") && length(o$covs) == 0)
                jmvcore::reject("Put at least one numeric covariate in the first covariate box to lag it.", code = "")
            sp <- sp_guard(sp_sample(self$data, c(o$dep, o$covs, o$noLag, o$timeVar), o))
            pd <- sp_guard(sp_panel_data(sp$df, o$idVar, o$timeVar))
            df <- pd$df; N <- pd$N; T <- pd$T
            # weights for the N units (one row per unit, in the order of pd$units)
            first <- df[match(pd$units, as.character(df[[o$idVar]])), , drop = FALSE]
            wt <- sp_guard(sp_weights(first, o))
            W <- wt$W
            sm <- sp_summary(wt$w$B)
            dz <- sp_guard(sp_design(df, o$covs, o$noLag))
            X <- dz$X[, -1, drop = FALSE]; lagcols <- dz$lagcols - 1
            y <- sp_num(df[[o$dep]])
            if (!identical(o$wStyle, "W") && N > 1500) jmvcore::reject("With more than 1500 units the spatial regression needs row-standardised weights.", code = "")

            note <- sp_weights_note(wt, 0)
            note <- paste0(note, sprintf("<p>Balanced panel: %d units x %d periods = %d observations.%s</p>", N, T, N * T,
                if (pd$n_unit_drop > 0 || sp$n_drop > 0) sprintf(" %d unit(s) not observed in every period and %d row(s) with missing values were dropped.", pd$n_unit_drop, sp$n_drop) else ""))
            if (sm$islands > 0) note <- paste0(note, "<p><b>Warning:</b> units without neighbours are in the sample.</p>")
            fit <- sp_guard(sp_panel_fit(y, X, W, N, T, model, eff, lagcols, wt$rowstd))
            if (isFALSE(fit$vcov_ok)) note <- paste0(note, "<p><b>Warning:</b> the covariance matrix could not be computed; standard errors are missing.</p>")
            if (isTRUE(fit$hit_bound)) note <- paste0(note, "<p><b>Warning:</b> a spatial parameter is at the edge of its admissible range; check the weights.</p>")
            note <- paste0(note, sprintf("<p>%s. Maximum likelihood (Elhorst, 2014); the spatial lags are computed before the within transformation.</p>",
                switch(eff, individual = "Unit fixed effects", time = "Time fixed effects", twoway = "Unit and time fixed effects", none = "Pooled model with a constant")))
            res$notes$setContent(note)

            cf <- fit$coef
            q <- if (identical(fit$stat_name, "t")) stats::qt(1 - (1 - o$level / 100) / 2, fit$dfres) else stats::qnorm(1 - (1 - o$level / 100) / 2)
            ord <- c(which(!cf$term %in% c("rho", "lambda")), which(cf$term %in% c("rho", "lambda")))
            lab <- c(rho = "rho (spatial lag of y)", lambda = "lambda (spatial error)")
            tb <- res$coef; tb$deleteRows()
            for (i in ord) {
                t <- cf$term[i]
                tb$addRow(rowKey = t, values = list(term = if (t %in% names(lab)) lab[[t]] else t, est = cf$est[i], se = cf$se[i], z = cf$z[i], p = cf$p[i],
                                                    lo = cf$est[i] - q * cf$se[i], hi = cf$est[i] + q * cf$se[i]))
            }
            tb$setNote("n", sprintf("Dependent variable: %s. The statistic is %s.", o$dep, if (identical(fit$stat_name, "t")) "t (degrees of freedom net of the fixed effects)" else "z (asymptotic)"))

            nobs <- N * T; ll <- fit$ll; k <- fit$k
            ft <- res$fit; ft$deleteRows()
            addf <- function(key, label, v, p = NULL) ft$addRow(rowKey = key, values = list(item = label, value = v, p = if (is.null(p)) NA else p))
            addf("N", "Units", N); addf("T", "Periods", T); addf("n", "Observations", nobs)
            addf("ll", "Log-likelihood (within-transformed data)", ll); addf("aic", "AIC", -2 * ll + 2 * k); addf("bic", "BIC", -2 * ll + k * log(nobs))
            addf("sigma", "Residual standard deviation", fit$sigma); addf("r2", "R-squared of the transformed model", fit$r2w)
            if (!is.null(fit$lr)) addf("lr", sprintf("Likelihood ratio test against the non-spatial model (%d df)", fit$lr_df), fit$lr, fit$lr_p)

            if (isTRUE(o$showEffects)) {
                et <- res$effectsTable; et$deleteRows()
                ev <- private$.effects(fit, X, lagcols, W, o)
                if (is.null(ev)) et$setNote("n", "Effects are defined for models with a spatial lag of y (SAR, SDM) or spatially lagged covariates (SLX, SDEM).")
                else {
                    for (v in seq_len(nrow(ev$est))) for (j in 1:3) {
                        kind <- c("Direct", "Indirect", "Total")[j]
                        se <- if (is.null(ev$se)) NA_real_ else ev$se[v, j]; z <- ev$est[v, j] / se
                        et$addRow(rowKey = paste(v, kind), values = list(term = ev$names[v], kind = kind, est = ev$est[v, j], se = se, z = z, p = if (is.na(z)) NA_real_ else 2 * stats::pnorm(-abs(z))))
                    }
                    et$setNote("n", ev$note)
                }
            }
        },

        .effects = function(fit, X, lagcols, W, o) {
            m <- fit$model
            if (!m %in% c("slx", "sar", "sdm", "sdem")) return(NULL)
            shift <- if (fit$effects == "none") 1 else 0
            nX <- ncol(X); cols <- (1 + shift):(nX + shift)
            tpos <- match(cols - shift, lagcols)
            vars <- data.frame(b = cols, t = ifelse(is.na(tpos), NA, nX + shift + tpos))
            if (m %in% c("sar", "sdm")) {
                if (m == "sar") vars$t <- NA
                e <- sp_effects(fit, vars, W, nsim = o$nsim, seed = o$seed)
                return(list(est = e$est, se = e$se, names = fit$pnames[cols],
                            note = if (is.null(e$se)) "Point estimates only (no simulations)." else sprintf("Average effects (LeSage and Pace, 2009); standard errors from %d simulations of rho and the coefficients.", e$nsim)))
            }
            V <- fit$vcov; B <- fit$beta
            ix <- function(i) if (is.null(rownames(V))) i else fit$pnames[i]       # the ML covariance also holds lambda and sigma2
            est <- matrix(NA_real_, nrow(vars), 3); se <- est
            for (i in seq_len(nrow(vars))) {
                b <- vars$b[i]; t <- vars$t[i]; th <- if (is.na(t)) 0 else B[t]
                est[i, ] <- c(B[b], th, B[b] + th)
                vb <- V[ix(b), ix(b)]
                if (is.na(t)) se[i, ] <- c(sqrt(vb), NA, sqrt(vb))
                else se[i, ] <- c(sqrt(vb), sqrt(V[ix(t), ix(t)]), sqrt(vb + V[ix(t), ix(t)] + 2 * V[ix(b), ix(t)]))
            }
            list(est = est, se = se, names = fit$pnames[cols], note = "Direct = own coefficient, indirect = coefficient of the spatially lagged covariate, total = their sum.")
        }
    )
)
