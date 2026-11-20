# jamovi analysis class: gravity model by Poisson pseudo maximum likelihood with high-dimensional fixed effects.

gravityClass <- R6::R6Class(
    "gravityClass",
    inherit = gravityBase,
    private = list(

        .run = function() {
            o <- self$options; res <- self$results
            sp_show_text(res, o, "gravity")
            if (is.null(o$dep) || length(o$covs) == 0 || is.null(o$origin) || is.null(o$dest)) {
                res$notes$setContent("<p>Select the flow, at least one covariate, the origin and the destination to begin.</p>")
                return()
            }
            if ((isTRUE(o$feOrigTime) || isTRUE(o$feDestTime) || isTRUE(o$feTime)) && is.null(o$timeVar))
                jmvcore::reject("Select the period variable for fixed effects that involve time.", code = "")
            vars <- unique(c(o$dep, o$covs, o$origin, o$dest, o$timeVar))
            df <- self$data[, vars, drop = FALSE]
            keep <- stats::complete.cases(df); n_drop <- sum(!keep); df <- df[keep, , drop = FALSE]
            n <- nrow(df)
            if (n < 20) jmvcore::reject("Fewer than 20 complete observations.", code = "")
            y <- sp_num(df[[o$dep]])
            if (any(y < 0)) jmvcore::reject("The flow must not be negative (zero flows are fine).", code = "")
            if (all(y == 0)) jmvcore::reject("All flows are zero.", code = "")
            dz <- sp_guard(sp_design(df, o$covs, NULL))
            X <- dz$X[, -1, drop = FALSE]
            g_o <- as.integer(factor(as.character(df[[o$origin]]))); g_d <- as.integer(factor(as.character(df[[o$dest]])))
            g_t <- if (!is.null(o$timeVar)) as.integer(factor(as.character(df[[o$timeVar]]))) else NULL
            pair <- as.integer(factor(paste(g_o, g_d, sep = "|")))
            if (anyDuplicated(if (is.null(g_t)) paste(g_o, g_d) else paste(g_o, g_d, g_t))) {
                dup <- TRUE
            } else dup <- FALSE
            fe <- list(); fl <- character(0)
            add <- function(g, label) { fe[[length(fe) + 1]] <<- g; fl <<- c(fl, label) }
            if (isTRUE(o$feOrigin)) add(g_o, "origin"); if (isTRUE(o$feDest)) add(g_d, "destination"); if (isTRUE(o$fePair)) add(pair, "origin-destination pair")
            if (isTRUE(o$feOrigTime)) add(as.integer(factor(paste(g_o, g_t))), "origin x period")
            if (isTRUE(o$feDestTime)) add(as.integer(factor(paste(g_d, g_t))), "destination x period")
            if (isTRUE(o$feTime)) add(g_t, "period")
            cl <- switch(o$cluster, robust = NULL, origin = list(g_o), dest = list(g_d), pair = list(pair), twoway = list(g_o, g_d))
            fit <- sp_guard(sp_ppml(y, X, fe, cl))

            note <- sprintf("<p>%d observations%s; %d used%s. Fixed effects: %s. Standard errors: %s.</p>", n,
                if (n_drop > 0) sprintf(" (%d with missing values dropped)", n_drop) else "", fit$n,
                if (fit$dropped > 0) sprintf(" (%d dropped because their fixed-effect group has only zero flows)", fit$dropped) else "",
                if (length(fl)) paste(fl, collapse = ", ") else "none (a constant is included)",
                switch(o$cluster, robust = "heteroskedasticity-robust", origin = "clustered by origin", dest = "clustered by destination", pair = "clustered by origin-destination pair", twoway = "two-way clustered by origin and destination"))
            if (dup && o$cluster == "pair") note <- paste0(note, "<p>Some origin-destination pairs appear in several rows (several periods), which is what clustering by pair is for.</p>")
            if (!fit$converged) note <- paste0(note, "<p><b>Warning:</b> the estimation did not converge in 100 iterations; check for separation or extreme values.</p>")
            if (!is.null(fit$nclusters) && min(fit$nclusters) < 10) note <- paste0(note, "<p><b>Warning:</b> fewer than 10 clusters; the clustered standard errors are unreliable.</p>")
            res$notes$setContent(note)

            q <- stats::qnorm(1 - (1 - o$level / 100) / 2)
            tb <- res$coef; tb$deleteRows()
            for (i in seq_len(nrow(fit$coef))) {
                b <- fit$coef$est[i]; se <- fit$coef$se[i]; z <- b / se
                tb$addRow(rowKey = i, values = list(term = fit$coef$term[i], est = b, se = se, z = z, p = 2 * stats::pnorm(-abs(z)), lo = b - q * se, hi = b + q * se, pct = 100 * (exp(b) - 1)))
            }
            tb$setNote("n", "Poisson pseudo maximum likelihood on the flow in levels (zeros are kept). For a covariate in logs the coefficient is an elasticity; for a dummy or a level covariate exp(b) - 1 is the proportional change in the expected flow (last column).")

            ft <- res$fit; ft$deleteRows()
            add2 <- function(k, v) ft$addRow(rowKey = k, values = list(item = k, value = v))
            add2("Observations used", fit$n); add2("Pseudo log-likelihood", fit$ll); add2("Deviance", fit$deviance); add2("Squared correlation of flow and fitted flow", fit$pr2)
            for (i in seq_along(fit$nfe)) add2(sprintf("Fixed effects: %s (groups)", fl[i]), fit$nfe[i])
            if (!is.null(fit$nclusters)) add2("Clusters (smallest)", min(fit$nclusters))
            add2("Iterations", fit$iterations)
        }
    )
)
