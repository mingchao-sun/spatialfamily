# jamovi analysis class: industrial agglomeration and regional specialisation indices.

aggClass <- R6::R6Class(
    "aggClass",
    inherit = aggBase,
    private = list(

        .run = function() {
            o <- self$options; res <- self$results
            sp_show_text(res, o, "agg")
            if (is.null(o$region) || is.null(o$industry) || is.null(o$value)) {
                res$notes$setContent("<p>Select the region, the industry and the employment (or output) variable to begin.</p>")
                return()
            }
            vars <- unique(c(o$region, o$industry, o$value, o$firms))
            df <- self$data[, vars, drop = FALSE]
            keep <- stats::complete.cases(df); n_drop <- sum(!keep); df <- df[keep, , drop = FALSE]
            v <- sp_num(df[[o$value]])
            if (any(v < 0)) jmvcore::reject("The employment (or output) variable must not be negative.", code = "")
            f <- if (!is.null(o$firms)) sp_num(df[[o$firms]]) else NULL
            if (!is.null(f) && any(f < 0)) jmvcore::reject("The number of establishments must not be negative.", code = "")
            M <- sp_guard(sp_agg_matrix(df[[o$region]], df[[o$industry]], v, f))
            E <- M$E
            if (nrow(E) < 3) jmvcore::reject("At least three regions are needed.", code = "")
            if (ncol(E) < 2) jmvcore::reject("At least two industries are needed.", code = "")
            if (sum(E) <= 0) jmvcore::reject("The total employment is zero.", code = "")
            dup <- sum(duplicated(paste(df[[o$region]], df[[o$industry]])))
            note <- sprintf("<p>%d regions, %d industries, total %s.%s%s</p>", nrow(E), ncol(E), format(sum(E), big.mark = ","),
                if (dup > 0) sprintf(" %d repeated region-industry rows were summed.", dup) else "",
                if (n_drop > 0) sprintf(" %d row(s) with missing values were dropped.", n_drop) else "")
            if (is.null(M$N)) note <- paste0(note, "<p>No establishment counts were given, so the plant Herfindahl H is taken as 0 and gamma = G / (1 - &sum;x<sup>2</sup>); this overstates concentration in industries with few large plants.</p>")
            else note <- paste0(note, "<p>H is approximated from the establishment counts, assuming equal-sized plants within a region-industry cell.</p>")
            res$notes$setContent(note)

            ind <- sp_agg_industry(E, M$N)
            ind <- ind[order(switch(o$sortBy, gamma = -ind$gamma, gini = -ind$gini, total = -ind$total, name = ind$industry)), , drop = FALSE]
            it <- res$industryTable; it$deleteRows()
            for (k in seq_len(nrow(ind))) {
                r <- ind[k, ]
                it$addRow(rowKey = k, values = list(industry = r$industry, total = r$total, share = r$share, nreg = as.integer(r$nreg), hhi = r$hhi, G = r$G,
                                                    H = if (is.na(r$H)) NA_real_ else r$H, gamma = r$gamma, gini = r$gini))
            }
            it$setNote("n", "G = sum of squared differences between the industry's regional shares and the regions' shares of total employment; gamma = Ellison-Glaeser index (about 0 = as concentrated as a random allocation, above 0.05 clearly concentrated, above 0.02 moderately); Gini = locational Gini of Krugman (0 = the industry is spread like total employment, 1 = all in one region).")

            if (isTRUE(o$showRegion)) {
                rg <- sp_agg_region(E)
                rt <- res$regionTable; rt$deleteRows()
                for (k in seq_len(nrow(rg))) rt$addRow(rowKey = k, values = as.list(rg[k, ]))
                rt$setNote("n", "Krugman index = sum over industries of |industry share in the region - industry share in the rest of the economy| (0 = same structure as the rest, 2 = completely different).")
            }

            if (isTRUE(o$saveLQ) && res$saveLQ$isNotFilled()) {
                LQ <- sp_agg_lq(E)
                # cell values for every kept row (repeated rows get the same value)
                lqv <- LQ[cbind(match(as.character(df[[o$region]]), rownames(E)), match(as.character(df[[o$industry]]), colnames(E)))]
                res$saveLQ$setRowNums(rownames(df)); res$saveLQ$setValues(as.numeric(lqv))
            }
        }
    )
)
