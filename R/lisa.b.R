# jamovi analysis class: local spatial autocorrelation (local Moran's I, Getis-Ord Gi / Gi*).

lisaClass <- R6::R6Class(
    "lisaClass",
    inherit = lisaBase,
    private = list(

        .run = function() {
            o <- self$options; res <- self$results
            sp_show_text(res, o, "lisa")
            if (is.null(o$var) || !sp_ready(o)) {
                res$notes$setContent(sp_ready_msg(o, if (is.null(o$var)) "Select a variable and set up the spatial weights. " else ""))
                return()
            }
            sp <- sp_guard(sp_sample(self$data, o$var, o))
            wt <- sp_guard(sp_weights(sp$df, o))
            df <- sp$df; n <- nrow(df)
            x <- sp_num(df[[o$var]])
            if (stats::sd(x) == 0) jmvcore::reject("The variable has no variation.", code = "")
            set.seed(o$seed)
            isG <- o$stat %in% c("gistar", "gi")
            if (!isG) {
                li <- sp_guard(sp_local_moran(x, wt$W, o$nperm, "two.sided"))
                stat <- li$Ii; zst <- li$z; p <- li$p
                Wz <- li$Wz; zc <- li$zs
            } else {
                lg <- sp_guard(sp_local_g(x, wt$w$B, o$wStyle, star = identical(o$stat, "gistar"), nperm = o$nperm, alt = "two.sided"))
                stat <- lg$G; zst <- lg$G; p <- lg$p
                zc <- x - mean(x); Wz <- as.numeric(wt$W %*% zc)
            }
            padj <- sp_adjust_p(p, o$adjust)
            if (!isG) cls <- sp_lisa_class(zc, Wz, padj, o$alpha)
            else cls <- ifelse(is.na(padj) | padj > o$alpha, "Not significant", ifelse(stat > 0, "Hot spot", "Cold spot"))
            levs <- if (isG) c("Hot spot", "Cold spot", "Not significant") else c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant")
            cls <- factor(cls, levels = levs)

            note <- sp_weights_note(wt, sp$n_drop)
            note <- paste0(note, sprintf("<p>%s, %d permutations (conditional on each unit's own value), significance level %s%s.</p>",
                switch(o$stat, localmoran = "Local Moran's I", gistar = "Getis-Ord Gi*", gi = "Getis-Ord Gi"), o$nperm, o$alpha,
                switch(o$adjust, none = ", no multiple-testing adjustment", fdr = ", Benjamini-Hochberg false discovery rate", bonferroni = ", Bonferroni adjustment")))

            ct <- res$counts; ct$deleteRows()
            tab <- table(cls)
            for (l in levs) ct$addRow(rowKey = l, values = list(cls = l, n = as.integer(tab[[l]]), pct = 100 * as.numeric(tab[[l]]) / n))

            if (isTRUE(o$showUnits)) {
                ut <- res$units; ut$deleteRows()
                ids <- if (!is.null(o$idVar)) as.character(df[[o$idVar]]) else rownames(df)
                lagx <- as.numeric(wt$W %*% x)
                for (i in seq_len(n))
                    ut$addRow(rowKey = i, values = list(id = ids[i], x = x[i], lag = lagx[i], stat = stat[i], z = zst[i], p = p[i], padj = padj[i], cls = as.character(cls[i])))
            }

            pal <- list(levs = levs, cols = unname(sp_cols_lisa[levs]))     # unnamed vectors: jamovi's state storage may drop names
            if (isTRUE(o$map)) {
                if (!is.null(wt$w$xy)) res$clusterMap$setState(list(st = sp_map_state(wt$w, value = as.character(cls)), pal = pal, label = if (isG) "Gi*" else "LISA"))
                else note <- paste0(note, "<p>A map needs coordinates or a boundary file.</p>")
            }
            if (isTRUE(o$scatter)) {
                zs <- zc / sqrt(sum(zc^2) / n)
                res$scatterPlot$setState(list(z = zs, Wz = as.numeric(wt$W %*% zs), cls = as.character(cls), name = o$var, pal = pal))
            }
            res$notes$setContent(note)

            code <- c("Not significant" = 0, "High-High" = 1, "Hot spot" = 1, "Low-Low" = 2, "Cold spot" = 2, "Low-High" = 3, "High-Low" = 4)
            rn <- rownames(df)
            if (isTRUE(o$saveCluster) && res$saveCluster$isNotFilled()) {
                res$saveCluster$setRowNums(rn); res$saveCluster$setValues(as.integer(code[as.character(cls)]))
            }
            if (isTRUE(o$saveStat) && res$saveStat$isNotFilled()) {
                res$saveStat$setRowNums(rn); res$saveStat$setValues(as.numeric(zst))
            }
            if (isTRUE(o$saveP) && res$saveP$isNotFilled()) {
                res$saveP$setRowNums(rn); res$saveP$setValues(as.numeric(padj))
            }
        },

        .plotMap = function(image, ggtheme, theme, ...) {
            s <- image$state
            if (is.null(s)) return(FALSE)
            print(sp_map_plot(s$st, ggtheme, discrete = TRUE, palette = stats::setNames(s$pal$cols, s$pal$levs), fill_label = s$label))
            TRUE
        },

        .plotScatter = function(image, ggtheme, theme, ...) {
            s <- image$state
            if (is.null(s)) return(FALSE)
            p <- sp_scatter_plot(s$z, s$Wz, ggtheme, s$name, cls = s$cls)
            print(p + ggplot2::scale_colour_manual(values = stats::setNames(s$pal$cols, s$pal$levs), name = NULL, drop = TRUE))
            TRUE
        }
    )
)
