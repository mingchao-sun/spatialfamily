# jamovi analysis class: build and describe the spatial weights.

swmClass <- R6::R6Class(
    "swmClass",
    inherit = swmBase,
    private = list(

        .run = function() {
            o <- self$options; res <- self$results
            sp_show_text(res, o, "swm")
            if (!sp_ready(o)) {
                res$notes$setContent(sp_ready_msg(o))
                return()
            }
            sp <- sp_guard(sp_sample(self$data, character(0), o))
            wt <- sp_guard(sp_weights(sp$df, o))
            sm <- sp_summary(wt$w$B)
            note <- sp_weights_note(wt, sp$n_drop)

            tb <- res$summary; tb$deleteRows()
            add <- function(k, v) tb$addRow(rowKey = k, values = list(item = k, value = v))
            add("Units", format(sm$n))
            add("Links (non-zero weights)", format(sm$links))
            add("Share of non-zero weights (%)", sp_fmt(sm$pct, 2))
            add("Neighbours per unit: mean", sp_fmt(sm$mean, 2))
            add("Neighbours per unit: minimum", format(sm$min))
            add("Neighbours per unit: maximum", format(sm$max))
            add("Units without neighbours", format(sm$islands))
            add("Connected components", sprintf("%d (largest has %d units)", sm$ncomp, max(tabulate(sm$comp))))
            add("Links are symmetric", if (sm$symmetric) "yes" else "no")
            add("Sum of the weights (S0)", sp_fmt(sum(wt$W), 3))

            nt <- res$nbTable; nt$deleteRows()
            tab <- table(sm$nnb)
            for (k in names(tab)) nt$addRow(rowKey = k, values = list(k = as.integer(k), n = as.integer(tab[[k]]), pct = 100 * as.numeric(tab[[k]]) / sm$n))

            if (isTRUE(o$linkPlot)) {
                if (!is.null(wt$w$xy)) res$linkMap$setState(sp_map_state(wt$w, links = TRUE))
                else note <- paste0(note, "<p>A map needs coordinates or a boundary file.</p>")
            }
            if (isTRUE(o$histPlot)) res$nbHist$setState(list(nnb = sm$nnb))
            res$notes$setContent(note)
        },

        .plotLinks = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(sp_map_plot(image$state, ggtheme))
            TRUE
        },

        .plotHist = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(sp_hist_plot(image$state$nnb, ggtheme))
            TRUE
        }
    )
)
