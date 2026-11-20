# Helpers shared by the spatial analyses. Nothing here depends on jamovi except sp_guard (jmvcore::reject).

sp_fmt <- function(x, d = 3) formatC(x, format = "f", digits = d, drop0trailing = TRUE)

# Turn engine errors into readable jamovi messages.
sp_guard <- function(expr) tryCatch(expr, error = function(e) jmvcore::reject(conditionMessage(e), code = ""))

sp_show_text <- function(results, opts, key) {
    if (isTRUE(opts$showGuide)) results$guide$setContent(sp_guide_html(key))
    if (isTRUE(opts$showHelp)) results$helpText$setContent(sp_notes_html(key))
}

sp_bt <- function(v) paste0("`", v, "`")

# The weights options as a plain list for sp_build_weights().
sp_wopts <- function(o) {
    list(source = o$wSource, file = o$wFile, idfield = if (is.null(o$wIdField)) "" else o$wIdField, method = o$wMethod,
         k = o$wK, dband = o$wDist, power = o$wPower, cutoff = o$wCutoff, order = o$wOrder, cumul = o$wCumul,
         symm = o$wSymm, lonlat = o$wLonLat, xvar = o$xVar, yvar = o$yVar, idvar = o$idVar, style = o$wStyle)
}

# Complete cases of the analysis variables, the coordinates and the id. Returns the data frame (jamovi row numbers kept as row names) and the count dropped.
sp_sample <- function(data, vars, o) {
    need <- vars
    if (identical(o$wSource, "coords")) {
        if (is.null(o$xVar) || is.null(o$yVar)) stop("Select the two coordinate variables (x / longitude and y / latitude).", call. = FALSE)
        need <- c(need, o$xVar, o$yVar)
    }
    if (!is.null(o$idVar)) need <- c(need, o$idVar)
    need <- unique(need)
    df <- data[, need, drop = FALSE]
    keep <- stats::complete.cases(df)
    n_drop <- sum(!keep)
    df <- df[keep, , drop = FALSE]
    if (nrow(df) < 10) stop("Fewer than 10 complete observations.", call. = FALSE)
    list(df = df, n_drop = n_drop)
}

# Builds the weights for the sample. Returns list(w = weights object, W = styled matrix, desc).
sp_weights <- function(df, o) {
    s <- sp_wopts(o)
    if (identical(s$method, "auto")) s$method <- if (identical(s$source, "boundary")) "queen" else "knn"
    w <- sp_build_weights(df, s)
    W <- sp_style(w$B, o$wStyle)
    list(w = w, W = W, rowstd = identical(o$wStyle, "W"), s = s)
}

sp_weights_note <- function(wt, n_drop) {
    sm <- sp_summary(wt$w$B)
    txt <- sprintf("Weights: %s; style %s; %d units, mean %.1f neighbours.", wt$w$desc,
                   switch(wt$s$style, W = "row-standardised", B = "as built (binary or distance-based)", C = "globally standardised", wt$s$style), sm$n, sm$mean)
    if (sm$islands > 0)
        txt <- paste0(txt, sprintf(" %d unit(s) have no neighbours; their local statistics are not defined and they add nothing to the global ones.", sm$islands))
    if (n_drop > 0) txt <- paste0(txt, sprintf(" %d observation(s) with missing values were dropped.", n_drop))
    paste0("<p>", txt, "</p>")
}

# ---------------------------------------------------------------- plots ----

sp_cols_lisa <- c("High-High" = "#d7191c", "Low-Low" = "#2c7bb6", "Low-High" = "#8fc3e6", "High-Low" = "#f4a261", "Hot spot" = "#d7191c", "Cold spot" = "#2c7bb6", "Not significant" = "#d9d9d9")

# Polygon table for ggplot2: one row per vertex with unit and ring ids.
sp_poly_df <- function(geom, value = NULL) {
    do.call(rbind, lapply(seq_along(geom), function(u) {
        do.call(rbind, lapply(seq_along(geom[[u]]), function(r) {
            m <- geom[[u]][[r]]
            data.frame(x = m[, 1], y = m[, 2], unit = u, ring = r, value = if (is.null(value)) NA else value[u])
        }))
    }))
}

sp_theme_map <- function() {
    ggplot2::theme(axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
}

sp_aspect <- function(xy, lonlat) {
    if (lonlat) 1 / cos(mean(xy[, 2]) * pi / 180) else 1
}

# state: list(xy, geom, lonlat, B (optional), value, discrete, title)
sp_map_state <- function(w, value = NULL, links = FALSE, maxlinks = 15000) {
    st <- list(xy = w$xy, lonlat = w$lonlat, value = value, geom = NULL, seg = NULL)
    if (!is.null(w$geom) && length(w$geom) <= 20000) st$poly <- sp_poly_df(w$geom, value)
    if (links && !is.null(w$xy)) {
        A <- (w$B != 0); A <- A | Matrix::t(A)
        T <- Matrix::summary(as(Matrix::triu(A, 1), "TsparseMatrix"))
        if (nrow(T) > maxlinks) T <- T[sort(sample.int(nrow(T), maxlinks)), , drop = FALSE]
        st$seg <- data.frame(x = w$xy[T$i, 1], y = w$xy[T$i, 2], xend = w$xy[T$j, 1], yend = w$xy[T$j, 2])
    }
    st
}

sp_map_plot <- function(st, ggtheme, discrete = TRUE, palette = sp_cols_lisa, fill_label = "", point_size = NULL) {
    p <- ggplot2::ggplot()
    asp <- sp_aspect(st$xy, st$lonlat)
    if (!is.null(st$poly)) {
        if (is.null(st$value)) {
            p <- p + ggplot2::geom_polygon(data = st$poly, ggplot2::aes(x = x, y = y, group = interaction(unit, ring)), fill = "#f2f2f2", colour = "#7f7f7f", linewidth = 0.2)
        } else {
            st$poly$value <- if (discrete) factor(st$poly$value, levels = names(palette)) else st$poly$value
            p <- p + ggplot2::geom_polygon(data = st$poly, ggplot2::aes(x = x, y = y, group = interaction(unit, ring), fill = value), colour = "#7f7f7f", linewidth = 0.15)
        }
    } else if (!is.null(st$value)) {
        d <- data.frame(x = st$xy[, 1], y = st$xy[, 2], value = if (discrete) factor(st$value, levels = names(palette)) else st$value)
        p <- p + ggplot2::geom_point(data = d, ggplot2::aes(x = x, y = y, colour = value), size = if (is.null(point_size)) max(0.8, min(3, 120 / sqrt(nrow(d)))) else point_size)
    } else {
        d <- data.frame(x = st$xy[, 1], y = st$xy[, 2])
        p <- p + ggplot2::geom_point(data = d, ggplot2::aes(x = x, y = y), size = max(0.6, min(2.5, 100 / sqrt(nrow(d)))), colour = "#333333")
    }
    if (!is.null(st$seg))
        p <- p + ggplot2::geom_segment(data = st$seg, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), colour = "#e07a1f", linewidth = 0.3, alpha = 0.8)
    p <- p + ggtheme + sp_theme_map()      # jamovi's theme object carries colour scales: add it first so that ours replace them
    if (!is.null(st$value)) {
        if (discrete) {
            p <- p + (if (!is.null(st$poly)) ggplot2::scale_fill_manual(values = palette, drop = FALSE, name = fill_label, na.value = "#ffffff")
                      else ggplot2::scale_colour_manual(values = palette, drop = FALSE, name = fill_label))
        } else {
            p <- p + (if (!is.null(st$poly)) ggplot2::scale_fill_gradient2(low = "#2c7bb6", mid = "#f7f7f7", high = "#d7191c", midpoint = 0, name = fill_label)
                      else ggplot2::scale_colour_gradient2(low = "#2c7bb6", mid = "#f7f7f7", high = "#d7191c", midpoint = 0, name = fill_label))
        }
    }
    p + ggplot2::coord_fixed(ratio = asp) + ggplot2::labs(x = NULL, y = NULL)
}

# Moran scatterplot: standardised value against its spatial lag; points coloured by a class when given.
sp_scatter_plot <- function(z, Wz, ggtheme, varname, slope = NULL, cls = NULL) {
    d <- data.frame(z = z, Wz = Wz)
    p <- ggplot2::ggplot(d, ggplot2::aes(x = z, y = Wz)) +
        ggplot2::geom_hline(yintercept = 0, colour = "#999999", linewidth = 0.3) + ggplot2::geom_vline(xintercept = 0, colour = "#999999", linewidth = 0.3)
    if (!is.null(cls)) {
        d$cls <- factor(cls, levels = names(sp_cols_lisa))
        p <- p + ggplot2::geom_point(data = d, ggplot2::aes(colour = cls), size = 1.8, alpha = 0.85) +
            ggplot2::scale_colour_manual(values = sp_cols_lisa, drop = TRUE, name = NULL)
    } else p <- p + ggplot2::geom_point(size = 1.8, alpha = 0.7, colour = "#2c7bb6")
    b <- if (is.null(slope)) unname(stats::coef(stats::lm(Wz ~ z, d))[2]) else slope
    p + ggplot2::geom_abline(intercept = 0, slope = b, colour = "#d7191c", linewidth = 0.7) +
        ggplot2::labs(x = paste0(varname, " (standardised)"), y = "Spatial lag of the standardised value") + ggtheme
}

# Is the weights source set up well enough to run?
sp_ready <- function(o) {
    if (identical(o$wSource, "coords")) return(!is.null(o$xVar) && !is.null(o$yVar))
    !is.null(o$idVar) && !is.null(o$wFile) && nzchar(trimws(o$wFile))
}

sp_ready_msg <- function(o, extra = "") {
    what <- if (identical(o$wSource, "coords")) "the two coordinate variables (X / longitude and Y / latitude)"
            else "the unit id variable and the path of the file"
    paste0("<p>", extra, "Select ", what, " under <i>Spatial weights</i> to begin.</p>")
}

sp_hist_plot <- function(nnb, ggtheme) {
    d <- data.frame(k = nnb)
    ggplot2::ggplot(d, ggplot2::aes(x = k)) +
        ggplot2::geom_bar(fill = "#2c7bb6", width = 0.8) +
        ggplot2::labs(x = "Number of neighbours", y = "Units") + ggtheme
}

# Map of a numeric value at point locations (GWR coefficients or t values), diverging scale centred at zero.
sp_value_map <- function(s, ggtheme) {
    d <- data.frame(x = s$x, y = s$y, v = s$v)
    asp <- sp_aspect(cbind(s$x, s$y), isTRUE(s$lonlat))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, colour = v)) +
        ggplot2::geom_point(size = max(1.2, min(4.5, 160 / sqrt(nrow(d)))), alpha = 0.95) + ggtheme + sp_theme_map()
    lim <- max(abs(d$v), na.rm = TRUE)
    if (isTRUE(s$t)) {
        p <- p + ggplot2::scale_colour_gradientn(colours = c("#2c7bb6", "#abd9e9", "#f2f2f2", "#fdae61", "#d7191c"), values = c(0, 0.3, 0.5, 0.7, 1),
                                                 limits = c(-max(lim, 2.5), max(lim, 2.5)), name = s$label)
    } else {
        p <- p + ggplot2::scale_colour_gradient2(low = "#2c7bb6", mid = "#f2f2f2", high = "#d7191c", midpoint = 0, name = s$label)
    }
    p + ggplot2::coord_fixed(ratio = asp) + ggplot2::labs(x = NULL, y = NULL)
}
