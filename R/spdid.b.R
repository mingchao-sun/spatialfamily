# jamovi analysis class: difference-in-differences with spatial spillovers (distance rings, lagged treatment, Conley standard errors).

cn <- function(v, nm) { m <- matrix(v, ncol = 1); colnames(m) <- nm; m }

spdidClass <- R6::R6Class(
    "spdidClass",
    inherit = spdidBase,
    private = list(

        .run = function() {
            o <- self$options; res <- self$results
            sp_show_text(res, o, "spdid")
            if (is.null(o$dep) || is.null(o$treat) || is.null(o$timeVar) || is.null(o$idVar)) {
                res$notes$setContent("<p>Select the outcome, the treatment indicator, the time period and the unit id to begin.</p>")
                return()
            }
            cuts <- numeric(0)
            if (nzchar(trimws(o$rings))) {
                cuts <- suppressWarnings(as.numeric(trimws(strsplit(gsub(";", ",", o$rings), ",")[[1]])))
                if (anyNA(cuts) || any(cuts <= 0) || is.unsorted(cuts, strictly = TRUE))
                    jmvcore::reject("Enter the upper limits of the distance rings as increasing positive numbers separated by commas, for example 5, 10, 20.", code = "")
            }
            needXY <- length(cuts) > 0 || identical(o$se, "conley")
            if (needXY && identical(o$wSource, "adjacency"))
                jmvcore::reject("Distance rings and Conley standard errors need coordinates or a boundary file (polygon centroids), not a neighbour list.", code = "")
            if (!length(cuts) && !isTRUE(o$lagTreat))
                jmvcore::reject("Enter distance rings and/or tick the spatially lagged treatment; otherwise this is an ordinary two-way fixed effects DID.", code = "")
            if (!sp_ready(o)) { res$notes$setContent(sp_ready_msg(o, "")); return() }
            if (identical(o$se, "conley") && o$conleyBw <= 0) jmvcore::reject("Enter a positive distance bandwidth for the Conley standard errors.", code = "")

            sp <- sp_guard(sp_sample(self$data, c(o$dep, o$treat, o$covs, o$timeVar), o))
            df <- sp$df
            tv <- df[[o$treat]]
            d <- if (is.factor(tv)) {
                if (nlevels(droplevels(tv)) != 2) jmvcore::reject("The treatment must be a 0/1 indicator (or a variable with two levels; the second is treated).", code = "")
                as.numeric(droplevels(tv) == levels(droplevels(tv))[2])
            } else sp_num(tv)
            if (!all(d %in% c(0, 1))) jmvcore::reject("The treatment must be 0 or 1: 1 when the unit is treated in that period.", code = "")
            if (all(d == 0) || all(d == 1)) jmvcore::reject("The treatment has no variation.", code = "")
            y <- sp_num(df[[o$dep]])

            u <- as.character(df[[o$idVar]])
            units <- sort(unique(if (is.numeric(df[[o$idVar]])) as.numeric(u) else u)); units <- as.character(units)
            N <- length(units)
            tt <- df[[o$timeVar]]
            tl <- if (is.factor(tt)) levels(droplevels(tt)) else as.character(sort(unique(tt)))
            T <- length(tl)
            ui <- match(u, units); ti <- match(as.character(tt), tl)
            if (anyDuplicated(paste(ui, ti))) jmvcore::reject("A unit appears more than once in the same period.", code = "")
            if (T < 2) jmvcore::reject("At least two periods are needed.", code = "")

            first <- df[match(seq_len(N), ui), , drop = FALSE]
            wt <- sp_guard(sp_weights(first, o))
            xy <- wt$w$xy
            if (needXY && is.null(xy)) jmvcore::reject("Coordinates or a boundary file are needed for distance rings and Conley standard errors.", code = "")

            TR <- matrix(NA_real_, N, T); TR[cbind(ui, ti)] <- d
            X <- cbind(Treatment = d)
            Dm <- NULL
            if (!is.null(xy) && needXY) Dm <- sp_guard(sp_dist_matrix(xy, wt$w$lonlat))
            ringCount <- NULL
            if (length(cuts)) {
                dmin <- sp_did_nearest(Dm, TR)
                ring <- sp_did_rings(dmin, TR, cuts)
                rv <- ring[cbind(ui, ti)]
                lab <- vapply(seq_along(cuts), function(k) sprintf("Ring %s-%s%s", if (k == 1) "0" else formatC(cuts[k - 1], format = "g"), formatC(cuts[k], format = "g"), if (wt$w$lonlat) " km" else ""), "")
                for (k in seq_along(cuts)) X <- cbind(X, cn(as.numeric(rv == k), lab[k]))
                ringCount <- vapply(seq_along(cuts), function(k) sum(rv == k), 1L)
                if (any(ringCount == 0)) jmvcore::reject(sprintf("No untreated unit-period falls in %s; change the ring limits.", lab[which(ringCount == 0)[1]]), code = "")
            }
            if (isTRUE(o$lagTreat)) {
                T0 <- TR; T0[is.na(T0)] <- 0
                WD <- as.matrix(wt$W %*% T0)
                X <- cbind(X, `W x Treatment` = WD[cbind(ui, ti)])
            }
            if (length(o$covs)) for (v in o$covs) X <- cbind(X, cn(sp_num(df[[v]]), v))

            fit <- sp_guard(sp_did_fit(y, X, ui, ti, o$se, Dm, if (identical(o$se, "conley")) o$conleyBw else NULL, o$conleyLag, o$effects))

            note <- sprintf("<p>%d observations, %d units, %d periods%s. %s fixed effects; standard errors %s.</p>", length(y), N, T,
                if (sp$n_drop > 0) sprintf(" (%d row(s) with missing values dropped)", sp$n_drop) else "",
                switch(o$effects, twoway = "Unit and time", unit = "Unit", time = "Time"),
                switch(o$se, cluster = sprintf("clustered by unit (%d clusters)", N), robust = "heteroskedasticity-robust",
                       conley = sprintf("Conley spatial HAC (distance bandwidth %s%s, serial lag %d)", formatC(o$conleyBw, format = "g"), if (wt$w$lonlat) " km" else "", o$conleyLag)))
            if (length(cuts)) note <- paste0(note, sprintf("<p>Ring effects compare untreated units whose nearest treated unit is within the ring, in that period, with untreated units farther than %s%s. A treated unit has no ring.</p>",
                formatC(max(cuts), format = "g"), if (wt$w$lonlat) " km" else ""))
            if (isTRUE(o$lagTreat)) note <- paste0(note, "<p>W x Treatment is the weighted share of a unit's neighbours that are treated (weights style as chosen).</p>")
            note <- paste0(note, "<p>With staggered adoption or effects that differ over time, two-way fixed effects can be biased (Goodman-Bacon, 2021); see the DiD module for robust estimators. Spillovers also contaminate the control group when the rings are too small.</p>")
            res$notes$setContent(note)

            q <- stats::qt(1 - (1 - o$level / 100) / 2, fit$dfres)
            cf <- fit$coef; tb <- res$coef; tb$deleteRows()
            for (i in seq_len(nrow(cf)))
                tb$addRow(rowKey = i, values = list(term = cf$term[i], est = cf$est[i], se = cf$se[i], z = cf$est[i] / cf$se[i], p = 2 * stats::pt(-abs(cf$est[i] / cf$se[i]), fit$dfres),
                                                    lo = cf$est[i] - q * cf$se[i], hi = cf$est[i] + q * cf$se[i]))
            tb$setNote("n", sprintf("Dependent variable: %s. t statistics with %d degrees of freedom (number of clusters or units minus one).", o$dep, fit$dfres))

            ft <- res$fit; ft$deleteRows()
            add <- function(k, v) ft$addRow(rowKey = k, values = list(item = k, value = v))
            add("Observations", fit$n); add("Units", N); add("Periods", T); add("Treated unit-periods", sum(d == 1)); add("Units ever treated", length(unique(ui[d == 1])))
            if (length(cuts)) for (k in seq_along(cuts)) add(sprintf("Unit-periods in %s", lab[k]), ringCount[k])
            add("R-squared of the transformed model", fit$r2w)
        }
    )
)
