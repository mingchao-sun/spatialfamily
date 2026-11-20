# jamovi analysis class: global spatial autocorrelation (Moran's I, Geary's C, Getis-Ord G).

moranClass <- R6::R6Class(
    "moranClass",
    inherit = moranBase,
    private = list(

        .run = function() {
            o <- self$options; res <- self$results
            sp_show_text(res, o, "moran")
            if (length(o$vars) == 0 || !sp_ready(o)) {
                res$notes$setContent(sp_ready_msg(o, if (length(o$vars) == 0) "Select one or more variables and set up the spatial weights. " else ""))
                return()
            }
            sp <- sp_guard(sp_sample(self$data, o$vars, o))
            wt <- sp_guard(sp_weights(sp$df, o))
            res$notes$setContent(sp_weights_note(wt, sp$n_drop))
            S <- sp_S(wt$W)
            set.seed(o$seed)
            tb <- res$global; tb$deleteRows()
            inf <- o$inference
            for (v in o$vars) {
                x <- sp_num(sp$df[[v]])
                if (stats::sd(x) == 0) jmvcore::reject(sprintf("The variable '%s' has no variation.", v), code = "")
                m <- sp_guard(sp_moran(x, wt$W, inf, o$alt, o$nperm, S))
                tb$addRow(rowKey = paste(v, "I"), values = list(var = v, stat = "Moran's I", est = m$stat, exp = m$expect, sd = sqrt(m$var), z = m$z, p = m$p))
                if (isTRUE(o$geary)) {
                    g <- sp_guard(sp_geary(x, wt$W, inf, o$alt, o$nperm, S))
                    tb$addRow(rowKey = paste(v, "C"), values = list(var = v, stat = "Geary's C", est = g$stat, exp = g$expect, sd = sqrt(g$var), z = g$z, p = g$p))
                }
                if (isTRUE(o$getisG)) {
                    gg <- tryCatch(sp_globalG(x, wt$W, if (inf == "permutation") "permutation" else "normality", o$alt, o$nperm, S),
                                   error = function(e) { tb$setNote(paste("G", v), sprintf("Getis-Ord G for '%s': %s", v, conditionMessage(e))); NULL })
                    if (!is.null(gg))
                        tb$addRow(rowKey = paste(v, "G"), values = list(var = v, stat = "Getis-Ord G", est = gg$stat, exp = gg$expect, sd = sqrt(gg$var), z = gg$z, p = gg$p))
                }
                if (isTRUE(o$scatter)) {
                    z <- (x - mean(x)) / stats::sd(x)
                    res$scatterPlots$get(key = v)$setState(list(z = z, Wz = as.numeric(wt$W %*% z), name = v))
                }
            }
            tb$setNote("inf", switch(inf,
                randomisation = "Inference: analytical, under randomisation.",
                normality = "Inference: analytical, under normality.",
                permutation = sprintf("Inference: %d random permutations (pseudo p-values).", o$nperm)))
        },

        .plotScatter = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            print(sp_scatter_plot(st$z, st$Wz, ggtheme, st$name))
            TRUE
        }
    )
)
