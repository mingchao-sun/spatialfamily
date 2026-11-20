# jamovi analysis class: spatial regression (maximum likelihood).

# Regressor matrix: intercept, then the covariates (factors as dummies). `lag` marks the columns that get a spatial lag.
sp_design <- function(df, covs, nolag) {
    vars <- c(covs, nolag)
    if (length(vars) == 0) stop("Select at least one covariate.", call. = FALSE)
    both <- intersect(covs, nolag)
    if (length(both)) stop("A variable cannot be in both covariate boxes.", call. = FALSE)
    cols <- list(`(Intercept)` = rep(1, nrow(df))); lag <- logical(0); src <- character(0)
    for (v in vars) {
        x <- df[[v]]
        if (is.factor(x)) {
            x <- droplevels(x)
            if (nlevels(x) < 2) stop(sprintf("The covariate '%s' has fewer than two levels.", v), call. = FALSE)
            mm <- stats::model.matrix(~x)[, -1, drop = FALSE]
            for (j in seq_len(ncol(mm))) { nm <- paste0(v, " (", levels(x)[j + 1], ")"); cols[[nm]] <- mm[, j]; lag[nm] <- FALSE; src[nm] <- v }
        } else {
            x <- as.numeric(x)
            if (stats::sd(x) == 0) stop(sprintf("The covariate '%s' has no variation.", v), call. = FALSE)
            cols[[v]] <- x; lag[v] <- v %in% covs; src[v] <- v
        }
    }
    X <- do.call(cbind, cols); colnames(X) <- names(cols)
    list(X = X, lagcols = which(c(FALSE, lag)), src = c(NA, src))
}

sregClass <- R6::R6Class(
    "sregClass",
    inherit = sregBase,
    private = list(

        .run = function() {
            o <- self$options; res <- self$results
            sp_show_text(res, o, "sreg")
            if (is.null(o$dep) || length(c(o$covs, o$noLag)) == 0 || !sp_ready(o)) {
                res$notes$setContent(sp_ready_msg(o, "Select a dependent variable, at least one covariate, and set up the spatial weights. "))
                return()
            }
            model <- o$model
            sp <- sp_guard(sp_sample(self$data, c(o$dep, o$covs, o$noLag), o))
            df <- sp$df; n <- nrow(df)
            wt <- sp_guard(sp_weights(df, o))
            y <- sp_num(df[[o$dep]])
            if (stats::sd(y) == 0) jmvcore::reject("The dependent variable has no variation.", code = "")
            dz <- sp_guard(sp_design(df, o$covs, o$noLag))
            X <- dz$X
            if (qr(X)$rank < ncol(X)) jmvcore::reject("The covariates are collinear; remove a redundant variable.", code = "")
            W <- wt$W
            sm <- sp_summary(wt$w$B)
            note <- sp_weights_note(wt, sp$n_drop)
            if (sm$islands > 0 && model %in% c("sar", "sdm", "sac", "sem", "sdem"))
                note <- paste0(note, "<p><b>Warning:</b> units without neighbours are in the sample; the spatial terms are not identified for them.</p>")
            if (!identical(o$wStyle, "W") && n > 1500)
                jmvcore::reject("With more than 1500 units the spatial regression needs row-standardised weights (set the style to Row-standardised).", code = "")

            ldeng <- NULL
            if (model %in% c("sar", "sem", "sdm", "sdem", "sac") || isTRUE(o$compare))
                ldeng <- sp_guard(sp_logdet_engine(W, wt$rowstd))
            fit <- sp_guard(sp_sreg(y, X, W, model, dz$lagcols, wt$rowstd, ldeng = ldeng))
            if (isFALSE(fit$vcov_ok)) note <- paste0(note, "<p><b>Warning:</b> the covariance matrix could not be computed (singular Hessian); standard errors are missing.</p>")
            if (isTRUE(fit$hit_bound)) note <- paste0(note, "<p><b>Warning:</b> a spatial parameter is at the edge of its admissible range; check the weights.</p>")
            if (!is.null(fit$ldmethod)) note <- paste0(note, sprintf("<p>Maximum likelihood; log-determinant by %s. Standard errors from the Hessian of the full log-likelihood.</p>", fit$ldmethod))

            # ---- coefficients
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
            tb$setNote("n", sprintf("Dependent variable: %s; %d observations. The statistic is %s.", o$dep, n,
                                    if (identical(fit$stat_name, "t")) "t (ordinary least squares standard errors)" else "z (asymptotic)"))

            # ---- fit
            ll <- fit$ll; k <- fit$k
            r2 <- stats::cor(y, fit$fitted)^2
            ft <- res$fit; ft$deleteRows()
            addf <- function(key, label, v, p = NULL) ft$addRow(rowKey = key, values = list(item = label, value = v, p = if (is.null(p)) NA else p))
            addf("n", "Observations", n); addf("ll", "Log-likelihood", ll); addf("aic", "AIC", -2 * ll + 2 * k); addf("bic", "BIC", -2 * ll + k * log(n))
            addf("sigma", "Residual standard deviation", fit$sigma); addf("r2", "Pseudo R-squared (squared correlation of y and fitted)", r2)
            if (!is.null(fit$lr)) addf("lr", sprintf("Likelihood ratio test against OLS (chi-squared, %d df)", fit$lr_df), fit$lr, fit$lr_p)

            # ---- LM diagnostics
            if (isTRUE(o$diag)) {
                lt <- res$lmTests; lt$deleteRows()
                for (i in seq_len(nrow(fit$diag))) lt$addRow(rowKey = i, values = as.list(fit$diag[i, ]))
                lt$setNote("n", "Tests on the residuals of the OLS regression of the dependent variable on the covariates (no spatial lags).")
            }

            # ---- effects
            if (isTRUE(o$effects)) {
                et <- res$effectsTable; et$deleteRows()
                ev <- private$.effects(fit, dz, W, o)
                if (is.null(ev)) {
                    et$setNote("n", "Direct, indirect and total effects are defined for models with a spatial lag (SAR, SDM, SAC) or spatially lagged covariates (SLX, SDEM).")
                } else {
                    for (v in seq_len(nrow(ev$est))) for (j in 1:3) {
                        kind <- c("Direct", "Indirect", "Total")[j]
                        se <- if (is.null(ev$se)) NA_real_ else ev$se[v, j]
                        z <- ev$est[v, j] / se
                        et$addRow(rowKey = paste(v, kind), values = list(term = ev$names[v], kind = kind, est = ev$est[v, j], se = se, z = z,
                                                                       p = if (is.na(z)) NA_real_ else 2 * stats::pnorm(-abs(z))))
                    }
                    et$setNote("n", ev$note)
                }
            }

            # ---- Moran's I of the residuals
            if (isTRUE(o$moranRes)) {
                mr <- sp_moran(fit$resid, W, "randomisation", "two.sided")
                res$moranTable$setRow(rowNo = 1, values = list(stat = mr$stat, exp = mr$expect, z = mr$z, p = mr$p))
                res$moranTable$setNote("n", "Residuals y - fitted (for models with a spatial lag the fitted values include the spatial lag of y). The test treats the residuals as an observed variable.")
            }

            # ---- comparison
            if (isTRUE(o$compare)) {
                cm <- res$compareTable; cm$deleteRows()
                labs <- c(ols = "OLS", slx = "SLX", sar = "SAR", sem = "SEM", sdm = "SDM", sdem = "SDEM", sac = "SAC")
                for (m in names(labs)) {
                    f <- tryCatch(sp_sreg(y, X, W, m, dz$lagcols, wt$rowstd, ldeng = ldeng), error = function(e) NULL)
                    if (is.null(f)) { cm$addRow(rowKey = m, values = list(model = labs[[m]])); next }
                    cm$addRow(rowKey = m, values = list(model = labs[[m]], rho = f$rho, lambda = f$lambda, ll = f$ll, k = f$k, aic = -2 * f$ll + 2 * f$k, bic = -2 * f$ll + f$k * log(n)))
                }
                cm$setNote("n", "All models use the same weights and covariates; SLX, SDM and SDEM lag the covariates in the first box. A lower AIC or BIC is better.")
            }

            res$notes$setContent(note)

            rn <- rownames(df)
            if (isTRUE(o$saveFitted) && res$saveFitted$isNotFilled()) { res$saveFitted$setRowNums(rn); res$saveFitted$setValues(as.numeric(fit$fitted)) }
            if (isTRUE(o$saveResid) && res$saveResid$isNotFilled()) { res$saveResid$setRowNums(rn); res$saveResid$setValues(as.numeric(fit$resid)) }
        },

        # Average direct / indirect / total effects; NULL when the model has none.
        .effects = function(fit, dz, W, o) {
            m <- fit$model
            if (!m %in% c("slx", "sar", "sdm", "sdem", "sac")) return(NULL)
            p0 <- ncol(dz$X); nm0 <- colnames(dz$X)
            cols <- 2:p0                                         # all but the intercept
            tpos <- match(cols, dz$lagcols)                      # position among the lagged columns, NA when not lagged
            vars <- data.frame(b = cols, t = ifelse(is.na(tpos), NA, p0 + tpos))
            if (m %in% c("sar", "sdm", "sac")) {
                if (m == "sar" || m == "sac") vars$t <- NA
                if (m == "sac") vars$t <- NA
                e <- sp_effects(fit, vars, W, nsim = o$nsim, seed = o$seed)
                return(list(est = e$est, se = e$se, names = nm0[cols],
                            note = paste0(if (is.null(e$se)) "Point estimates only (no simulations)." else
                                   sprintf("Average effects (LeSage and Pace, 2009); standard errors from %d simulations of rho and the coefficients.", e$nsim),
                                   if (is.null(fit$ldeng$ev)) " With more than 1500 units the direct effects use a stochastic trace approximation (error below 1%)." else "")))
            }
            # SLX / SDEM: linear in the coefficients
            V <- fit$vcov
            est <- matrix(NA_real_, nrow(vars), 3); se <- est
            B <- fit$beta
            for (i in seq_len(nrow(vars))) {
                b <- vars$b[i]; t <- vars$t[i]
                th <- if (is.na(t)) 0 else B[t]
                est[i, ] <- c(B[b], th, B[b] + th)
                if (isTRUE(fit$vcov_ok) || m == "slx") {
                    vb <- if (m == "slx") V[b, b] else V[nm0[b], nm0[b]]
                    if (is.na(t)) { se[i, ] <- c(sqrt(vb), NA, sqrt(vb)) }
                    else {
                        nmt <- fit$pnames[t]
                        vt <- if (m == "slx") V[t, t] else V[nmt, nmt]
                        cv <- if (m == "slx") V[b, t] else V[nm0[b], nmt]
                        se[i, ] <- c(sqrt(vb), sqrt(vt), sqrt(vb + vt + 2 * cv))
                    }
                }
            }
            list(est = est, se = se, names = nm0[cols],
                 note = "Direct = own coefficient, indirect = coefficient of the spatially lagged covariate, total = their sum.")
        }
    )
)
