# jamovi analysis class: geographically weighted regression.

gwrClass <- R6::R6Class(
    "gwrClass",
    inherit = gwrBase,
    private = list(

        .run = function() {
            o <- self$options; res <- self$results
            sp_show_text(res, o, "gwr")
            if (is.null(o$dep) || length(o$covs) == 0 || is.null(o$xVar) || is.null(o$yVar)) {
                res$notes$setContent("<p>Select a dependent variable, at least one covariate and the two coordinate variables to begin.</p>")
                return()
            }
            vars <- c(o$dep, o$covs, o$xVar, o$yVar, o$idVar)
            df <- self$data[, unique(vars), drop = FALSE]
            keep <- stats::complete.cases(df); n_drop <- sum(!keep); df <- df[keep, , drop = FALSE]
            n <- nrow(df)
            if (n < 20) jmvcore::reject("Fewer than 20 complete observations.", code = "")
            y <- sp_num(df[[o$dep]])
            Xc <- vapply(o$covs, function(v) sp_num(df[[v]]), numeric(n))
            Xc <- matrix(Xc, n, dimnames = list(NULL, o$covs))
            if (any(apply(Xc, 2, stats::sd) == 0)) jmvcore::reject("A covariate has no variation.", code = "")
            X <- cbind(`(Intercept)` = 1, Xc)
            xy <- cbind(sp_num(df[[o$xVar]]), sp_num(df[[o$yVar]]))
            if (o$wLonLat && (any(abs(xy[, 1]) > 180) || any(abs(xy[, 2]) > 90)))
                jmvcore::reject("The coordinates are not valid longitude / latitude degrees; untick the longitude/latitude option for projected coordinates.", code = "")
            if (anyDuplicated(xy)) {
                dup <- sum(duplicated(xy))
                if (dup > 0.2 * n) jmvcore::reject("Many units share the same location; geographically weighted regression needs distinct locations.", code = "")
            }
            adaptive <- identical(o$bwType, "adaptive")
            bw <- if (identical(o$bwMethod, "manual")) o$bw else NULL
            if (!is.null(bw) && bw <= 0) jmvcore::reject("Enter a positive bandwidth.", code = "")
            fit <- sp_guard(sp_gwr(y, X, xy, o$wLonLat, o$kernel, adaptive, bw, if (identical(o$bwMethod, "cv")) "cv" else "aicc"))

            note <- sprintf("<p>%d observations%s. %s kernel, %s bandwidth %s%s.</p>", n,
                if (n_drop > 0) sprintf(" (%d with missing values dropped)", n_drop) else "",
                switch(o$kernel, bisquare = "Bi-square", gaussian = "Gaussian", exponential = "Exponential", tricube = "Tricube", boxcar = "Box-car"),
                if (adaptive) "adaptive" else "fixed", if (adaptive) sprintf("%d nearest neighbours", round(fit$bw)) else paste0(formatC(fit$bw, digits = 5, format = "g"), if (o$wLonLat) " km" else " (coordinate units)"),
                if (is.null(bw)) sprintf(", selected by %s", if (identical(o$bwMethod, "cv")) "cross-validation" else "AICc") else " (set manually)")
            if (adaptive && fit$bw < 2 * ncol(X) + 5) note <- paste0(note, "<p><b>Warning:</b> the bandwidth is small relative to the number of coefficients; local estimates may be unstable.</p>")
            res$notes$setContent(note)

            nm <- colnames(X); p <- ncol(X)
            st <- res$summary; st$deleteRows()
            for (j in seq_len(p)) {
                b <- fit$beta[, j]; qs <- stats::quantile(b, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
                st$addRow(rowKey = nm[j], values = list(term = nm[j], ols = fit$ols_coef[j], olsse = fit$ols_se[j], min = qs[1], q1 = qs[2], med = qs[3], q3 = qs[4], max = qs[5],
                                                        nsig = 100 * mean(abs(fit$tval[, j]) > 1.96)))
            }
            st$setNote("n", "Quartiles of the local coefficients over the units; the last column is the share of units where the local coefficient is significant at about 5% (|t| > 1.96; not adjusted for multiple testing).")

            dt <- res$diag; dt$deleteRows()
            add <- function(k, a, b = NA) dt$addRow(rowKey = k, values = list(item = k, gwr = a, ols = b))
            add("Observations", n, n); add("Bandwidth", fit$bw); add("Effective number of parameters", fit$enp, p); add("Residual sum of squares", fit$rss, sum((y - X %*% fit$ols_coef)^2))
            add("R-squared", fit$r2, fit$ols_r2); add("AICc", fit$aicc, fit$ols_aicc)
            dt$setNote("n", "A lower AICc than OLS means the spatially varying model fits better after the extra flexibility is penalised. The effective number of parameters is 2 tr(S) - tr(S'S).")

            if (isTRUE(o$showLocal)) {
                lt <- res$local; lt$deleteRows()
                for (j in seq_len(p)) {
                    lt$addColumn(name = paste0("b", j), title = nm[j], type = "number")
                    lt$addColumn(name = paste0("t", j), title = paste("t:", nm[j]), type = "number")
                }
                ids <- if (!is.null(o$idVar)) as.character(df[[o$idVar]]) else rownames(df)
                for (i in seq_len(n)) {
                    v <- list(id = ids[i], x = xy[i, 1], y = xy[i, 2], fitted = fit$yhat[i], r2 = fit$localR2[i])
                    for (j in seq_len(p)) { v[[paste0("b", j)]] <- fit$beta[i, j]; v[[paste0("t", j)]] <- fit$tval[i, j] }
                    lt$addRow(rowKey = i, values = v)
                }
            }

            if (isTRUE(o$maps)) {
                stat <- if (identical(o$mapStat, "t")) fit$tval else fit$beta
                lab <- if (identical(o$mapStat, "t")) "t value" else "coefficient"
                for (j in 2:p) res$mapPlots$get(key = nm[j])$setState(list(x = xy[, 1], y = xy[, 2], v = stat[, j], label = paste(nm[j], lab), lonlat = o$wLonLat, t = identical(o$mapStat, "t")))
                res$mapInt$setState(list(x = xy[, 1], y = xy[, 2], v = stat[, 1], label = paste("Intercept", lab), lonlat = o$wLonLat, t = identical(o$mapStat, "t")))
            }

            rn <- rownames(df)
            if (isTRUE(o$saveR2) && res$saveR2$isNotFilled()) { res$saveR2$setRowNums(rn); res$saveR2$setValues(as.numeric(fit$localR2)) }
            if (isTRUE(o$saveResid) && res$saveResid$isNotFilled()) { res$saveResid$setRowNums(rn); res$saveResid$setValues(as.numeric(fit$resid)) }
        },

        .plotMap = function(image, ggtheme, theme, ...) {
            s <- image$state
            if (is.null(s)) return(FALSE)
            print(sp_value_map(s, ggtheme))
            TRUE
        },

        .plotMapInt = function(image, ggtheme, theme, ...) {
            s <- image$state
            if (is.null(s)) return(FALSE)
            print(sp_value_map(s, ggtheme))
            TRUE
        }
    )
)
