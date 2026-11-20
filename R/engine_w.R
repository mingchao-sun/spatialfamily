# Spatial weights engine. Base R and Matrix only; nothing here depends on jamovi.
#
# A weights object is a list:
#   B       sparse n x n matrix of raw link weights (binary, or distance decay), zero diagonal
#   ids     character vector of unit ids (row order of the estimation sample)
#   xy      n x 2 matrix of coordinates (or centroids), or NULL
#   geom    list of polygon rings per unit (boundary files only), or NULL
#   lonlat  TRUE when xy are degrees
#   desc    one line that describes how the weights were built
# sp_style() turns B into the matrix used by the statistics.

sp_num <- function(x) if (is.factor(x)) as.numeric(as.character(x)) else as.numeric(x)

# ---------------------------------------------------------------- distances ----

sp_dist_from <- function(xy, i, lonlat) {
    if (lonlat) {
        r <- pi / 180
        lon <- xy[, 1] * r; lat <- xy[, 2] * r
        a <- sin((lat - lat[i]) / 2)^2 + cos(lat[i]) * cos(lat) * sin((lon - lon[i]) / 2)^2
        6371.0088 * 2 * asin(pmin(1, sqrt(a)))
    } else sqrt((xy[, 1] - xy[i, 1])^2 + (xy[, 2] - xy[i, 2])^2)
}

# Links from coordinates: "knn", "dist" (distance band) or "idw" (inverse distance within a cutoff).
sp_links_coords <- function(xy, method, k = 5, dband = 0, power = 1, cutoff = 0, lonlat = FALSE) {
    n <- nrow(xy)
    if (n < 3) stop("At least 3 units are needed to build spatial weights.", call. = FALSE)
    ii <- vector("list", n); jj <- vector("list", n); xx <- vector("list", n)
    desc <- ""
    if (method == "knn") {
        if (k < 1 || k >= n) stop(sprintf("The number of nearest neighbours must be between 1 and %d.", n - 1), call. = FALSE)
        for (i in seq_len(n)) {
            d <- sp_dist_from(xy, i, lonlat); d[i] <- Inf
            nb <- order(d)[seq_len(k)]
            ii[[i]] <- rep.int(i, k); jj[[i]] <- nb; xx[[i]] <- rep.int(1, k)
        }
        desc <- sprintf("%d nearest neighbours", k)
    } else if (method == "dist") {
        thr <- dband
        if (thr <= 0) {
            thr <- 0
            for (i in seq_len(n)) { d <- sp_dist_from(xy, i, lonlat); d[i] <- Inf; thr <- max(thr, min(d)) }
            thr <- thr * (1 + 1e-9)
            desc <- sprintf("distance band, automatic threshold %s (smallest value that leaves no unit without a neighbour)",
                            formatC(thr, format = "g", digits = 6))
        } else desc <- sprintf("distance band, threshold %s", formatC(thr, format = "g", digits = 6))
        for (i in seq_len(n)) {
            d <- sp_dist_from(xy, i, lonlat); d[i] <- Inf
            nb <- which(d <= thr)
            ii[[i]] <- rep.int(i, length(nb)); jj[[i]] <- nb; xx[[i]] <- rep.int(1, length(nb))
        }
    } else if (method == "idw") {
        if (cutoff <= 0 && n > 3000)
            stop("Inverse distance weights without a distance cutoff make a dense matrix; set a cutoff for more than 3000 units.", call. = FALSE)
        if (power <= 0) stop("The distance decay exponent must be positive.", call. = FALSE)
        for (i in seq_len(n)) {
            d <- sp_dist_from(xy, i, lonlat); d[i] <- Inf
            nb <- if (cutoff > 0) which(d <= cutoff) else seq_len(n)[-i]
            if (any(d[nb] <= 0)) stop("Two units have the same location, so inverse distance is undefined; remove duplicates or use nearest neighbours.", call. = FALSE)
            ii[[i]] <- rep.int(i, length(nb)); jj[[i]] <- nb; xx[[i]] <- d[nb]^(-power)
        }
        desc <- sprintf("inverse distance, exponent %s%s", formatC(power, format = "g"), if (cutoff > 0) sprintf(", cutoff %s", formatC(cutoff, format = "g", digits = 6)) else ", no cutoff")
    } else stop("Unknown weights method for coordinates.", call. = FALSE)
    B <- Matrix::sparseMatrix(i = unlist(ii), j = unlist(jj), x = unlist(xx), dims = c(n, n))
    list(B = B, desc = desc)
}

# ------------------------------------------------------------ graph helpers ----

sp_symmetrize <- function(B) {
    # keep the larger weight of the two directions
    Bt <- Matrix::t(B)
    as((B + Bt + abs(B - Bt)) / 2, "CsparseMatrix")
}

# Neighbours of order `order` (the units that need exactly `order` steps); cumulative = also the lower orders.
sp_order <- function(B, order = 1, cumulative = FALSE) {
    if (order <= 1) return(B)
    A <- (B != 0) * 1; A <- as(A, "CsparseMatrix")
    A <- A + Matrix::t(A); A@x[] <- 1; A <- as(A, "CsparseMatrix")
    reach <- A; seen <- A; levels <- list(A)
    for (o in 2:order) {
        nxt <- (reach %*% A); nxt@x[] <- 1
        nxt <- nxt - nxt * seen; nxt <- nxt - Matrix::Diagonal(nrow(A), Matrix::diag(nxt))
        nxt@x[nxt@x < 0] <- 0; nxt <- Matrix::drop0(nxt)
        seen <- seen + nxt; seen@x[] <- 1
        levels[[o]] <- nxt; reach <- nxt
        if (length(nxt@x) == 0) break
    }
    out <- if (cumulative) Reduce(`+`, levels) else levels[[length(levels)]]
    out@x[] <- 1
    as(Matrix::drop0(out), "CsparseMatrix")
}

# Row-standardised ("W"), as supplied ("B") or globally standardised ("C") weights.
sp_style <- function(B, style = "W") {
    B <- as(B, "CsparseMatrix")
    n <- nrow(B)
    if (style == "W") {
        rs <- Matrix::rowSums(B); inv <- ifelse(rs > 0, 1 / rs, 0)
        Matrix::Diagonal(n, inv) %*% B
    } else if (style == "C") {
        B * (n / sum(B))
    } else B
}

sp_neighbours <- function(B) {
    T <- Matrix::summary(as(B, "TsparseMatrix"))
    split(T$j, factor(T$i, levels = seq_len(nrow(B))))
}

# Connected components of the (symmetrised) link graph.
sp_components <- function(B) {
    n <- nrow(B)
    A <- (B != 0) * 1; A <- A + Matrix::t(A)
    nb <- sp_neighbours(A)
    comp <- integer(n); k <- 0L
    for (s in seq_len(n)) {
        if (comp[s] != 0L) next
        k <- k + 1L; comp[s] <- k; stack <- s
        while (length(stack)) {
            v <- stack[length(stack)]; stack <- stack[-length(stack)]
            w <- nb[[v]]; w <- w[comp[w] == 0L]
            if (length(w)) { comp[w] <- k; stack <- c(stack, w) }
        }
    }
    comp
}

sp_summary <- function(B) {
    n <- nrow(B)
    nnb <- as.numeric(Matrix::rowSums(B != 0))
    comp <- sp_components(B)
    list(n = n, links = sum(B != 0), pct = 100 * sum(B != 0) / (n * n),
         mean = mean(nnb), min = min(nnb), max = max(nnb), nnb = nnb,
         islands = sum(nnb == 0), ncomp = max(comp), comp = comp,
         symmetric = isTRUE(all(abs(B - Matrix::t(B)) < 1e-12 * max(1, max(abs(B@x), 0)))))
}

# -------------------------------------------------------------- file readers --

sp_read_dbf <- function(path) {
    con <- file(path, "rb"); on.exit(close(con))
    raw <- readBin(con, "raw", file.info(path)$size)
    nrec <- readBin(raw[5:8], "integer", size = 4, endian = "little")
    hlen <- readBin(raw[9:10], "integer", size = 2, signed = FALSE, endian = "little")
    rlen <- readBin(raw[11:12], "integer", size = 2, signed = FALSE, endian = "little")
    nf <- (hlen - 33) %/% 32
    fn <- character(nf); ft <- character(nf); fl <- integer(nf)
    for (f in seq_len(nf)) {
        o <- 32 + (f - 1) * 32
        fn[f] <- trimws(rawToChar(raw[(o + 1):(o + 11)][raw[(o + 1):(o + 11)] != as.raw(0)]))
        ft[f] <- rawToChar(raw[o + 12]); fl[f] <- as.integer(raw[o + 17])
    }
    out <- vector("list", nf); names(out) <- fn
    start <- cumsum(c(1, fl))[seq_len(nf)]          # offset inside a record, after the deletion flag
    recs <- matrix(raw[hlen + seq_len(nrec * rlen)], nrow = rlen)
    for (f in seq_len(nf)) {
        m <- recs[1 + start[f] - 1 + seq_len(fl[f]), , drop = FALSE]
        s <- vapply(seq_len(nrec), function(r) { b <- m[, r]; b <- b[b != as.raw(0)]; trimws(rawToChar(b)) }, "")
        s <- enc2utf8(iconv(s, "UTF-8", "UTF-8", sub = ""))
        out[[f]] <- if (ft[f] %in% c("N", "F")) suppressWarnings(as.numeric(s)) else s
    }
    as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

# Polygon shapefile: a list of units, each a list of rings (numeric matrices with columns x, y).
sp_read_shp <- function(path) {
    con <- file(path, "rb"); on.exit(close(con))
    raw <- readBin(con, "raw", file.info(path)$size)
    if (length(raw) < 100 || readBin(raw[1:4], "integer", size = 4, endian = "big") != 9994L)
        stop("This does not look like a shapefile (.shp).", call. = FALSE)
    stype <- readBin(raw[33:36], "integer", size = 4, endian = "little")
    if (!stype %in% c(5, 15, 25, 31))
        stop(sprintf("The shapefile has geometry type %d; only polygon files are supported.", stype), call. = FALSE)
    pos <- 101; out <- list()
    while (pos + 8 <= length(raw)) {
        clen <- readBin(raw[(pos + 4):(pos + 7)], "integer", size = 4, endian = "big") * 2
        body <- raw[(pos + 8):(pos + 7 + clen)]
        pos <- pos + 8 + clen
        t <- readBin(body[1:4], "integer", size = 4, endian = "little")
        if (t == 0) { out[[length(out) + 1]] <- list(); next }
        np <- readBin(body[37:40], "integer", size = 4, endian = "little")
        npt <- readBin(body[41:44], "integer", size = 4, endian = "little")
        parts <- readBin(body[45:(44 + 4 * np)], "integer", n = np, size = 4, endian = "little")
        o <- 44 + 4 * np
        pts <- readBin(body[(o + 1):(o + 16 * npt)], "double", n = 2 * npt, size = 8, endian = "little")
        pm <- matrix(pts, ncol = 2, byrow = TRUE)
        ends <- c(parts[-1], npt)
        out[[length(out) + 1]] <- lapply(seq_len(np), function(p) pm[(parts[p] + 1):ends[p], , drop = FALSE])
    }
    dbf <- sub("\\.shp$", ".dbf", path, ignore.case = TRUE)
    if (!file.exists(dbf)) dbf <- sub("\\.shp$", ".DBF", path)
    attrs <- if (file.exists(dbf)) sp_read_dbf(dbf) else NULL
    if (!is.null(attrs) && nrow(attrs) != length(out)) attrs <- NULL
    list(units = out, attrs = attrs)
}

sp_read_geojson <- function(path) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Reading GeoJSON needs the jsonlite package.", call. = FALSE)
    g <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    feats <- if (identical(g$type, "FeatureCollection")) g$features else if (identical(g$type, "Feature")) list(g) else stop("The file is not a GeoJSON FeatureCollection.", call. = FALSE)
    ring <- function(r) { m <- matrix(unlist(lapply(r, function(p) c(p[[1]], p[[2]]))), ncol = 2, byrow = TRUE); m }
    units <- lapply(feats, function(f) {
        gm <- f$geometry
        if (is.null(gm)) return(list())
        if (identical(gm$type, "Polygon")) lapply(gm$coordinates, ring)
        else if (identical(gm$type, "MultiPolygon")) unlist(lapply(gm$coordinates, function(p) lapply(p, ring)), recursive = FALSE)
        else stop(sprintf("GeoJSON geometry '%s' is not supported; use Polygon or MultiPolygon features.", gm$type), call. = FALSE)
    })
    props <- lapply(feats, function(f) { p <- f$properties; p[vapply(p, function(v) length(v) == 1 && !is.list(v), TRUE)] })
    keys <- unique(unlist(lapply(props, names)))
    attrs <- if (length(keys)) {
        as.data.frame(stats::setNames(lapply(keys, function(k) {
            v <- lapply(props, function(p) if (is.null(p[[k]])) NA else p[[k]])
            unlist(v)
        }), keys), stringsAsFactors = FALSE, check.names = FALSE)
    } else NULL
    if (!is.null(attrs) && nrow(attrs) != length(units)) attrs <- NULL
    list(units = units, attrs = attrs)
}

sp_read_boundary <- function(path) {
    if (!file.exists(path)) stop(sprintf("The file was not found: %s", path), call. = FALSE)
    ext <- tolower(tools::file_ext(path))
    b <- if (ext == "shp") sp_read_shp(path) else if (ext %in% c("geojson", "json")) sp_read_geojson(path)
         else stop("The boundary file must be a shapefile (.shp, with its .dbf next to it) or a GeoJSON file (.geojson / .json).", call. = FALSE)
    if (any(lengths(b$units) == 0)) stop("The boundary file has units without geometry.", call. = FALSE)
    b
}

# Edge list (csv / txt: from, to [, weight]), GAL or GWT file. Returns a data frame from, to, w (ids as character).
sp_read_edges <- function(path) {
    if (!file.exists(path)) stop(sprintf("The file was not found: %s", path), call. = FALSE)
    ext <- tolower(tools::file_ext(path))
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    lines <- gsub("[\"']", "", lines[nzchar(trimws(lines))])
    if (ext == "gal") {
        hd <- strsplit(trimws(lines[1]), "[[:space:],]+")[[1]]
        rest <- lines[-1]
        if (length(rest) %% 2 != 0) stop("The GAL file is not in the expected format.", call. = FALSE)
        from <- character(0); to <- character(0)
        for (r in seq(1, length(rest), by = 2)) {
            h <- strsplit(trimws(rest[r]), "[[:space:],]+")[[1]]
            nb <- strsplit(trimws(rest[r + 1]), "[[:space:],]+")[[1]]
            if (as.integer(h[2]) > 0) { from <- c(from, rep(h[1], length(nb))); to <- c(to, nb) }
        }
        return(data.frame(from = from, to = to, w = 1, stringsAsFactors = FALSE))
    }
    if (ext == "gwt") lines <- lines[-1]
    parts <- strsplit(trimws(lines), "[[:space:],;\t]+")
    nc <- vapply(parts, length, 1L)
    if (any(nc < 2)) stop("Each line of the edge list needs at least two ids (from, to).", call. = FALSE)
    first <- parts[[1]]
    f1 <- tolower(first[1]); f2 <- tolower(first[2])
    hasHeader <- ext != "gwt" && f1 %in% c("from", "source", "origin", "orig", "i", "id", "unit", "node", "src", "focal") &&
        f2 %in% c("to", "target", "dest", "destination", "j", "neighbor", "neighbour", "nb", "neighbors", "neighbours")
    if (hasHeader) parts <- parts[-1]
    from <- vapply(parts, `[`, "", 1); to <- vapply(parts, `[`, "", 2)
    w <- vapply(parts, function(p) if (length(p) >= 3) suppressWarnings(as.numeric(p[3])) else 1, 1)
    if (anyNA(w)) stop("A weight in the edge list is not a number.", call. = FALSE)
    data.frame(from = from, to = to, w = w, stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------ polygons --

sp_ring_area <- function(m) {
    x <- m[, 1]; y <- m[, 2]; n <- nrow(m)
    s <- x * c(y[-1], y[1]) - c(x[-1], x[1]) * y
    sum(s) / 2
}

# Area-weighted centroid of all rings (signed areas: outer rings and holes have opposite orientation).
sp_centroids <- function(units) {
    t(vapply(units, function(rings) {
        A <- 0; cx <- 0; cy <- 0
        for (m in rings) {
            if (nrow(m) < 3) next
            x <- m[, 1]; y <- m[, 2]; x1 <- c(x[-1], x[1]); y1 <- c(y[-1], y[1])
            cr <- x * y1 - x1 * y; a <- sum(cr) / 2
            A <- A + a; cx <- cx + sum((x + x1) * cr) / 6; cy <- cy + sum((y + y1) * cr) / 6
        }
        if (abs(A) < 1e-14) { m <- do.call(rbind, rings); return(colMeans(m)) }
        c(cx / A, cy / A)
    }, c(0, 0)))
}

# Queen (shared vertex) or rook (shared edge) contiguity between polygons.
sp_contig <- function(units, type = "queen") {
    n <- length(units)
    allv <- do.call(rbind, lapply(units, function(rings) do.call(rbind, rings)))
    tol <- max(diff(range(allv[, 1])), diff(range(allv[, 2])), 1) * 1e-9
    f <- integer(0); ka <- character(0); kb <- character(0)
    for (u in seq_len(n)) for (m in units[[u]]) {
        if (nrow(m) < 2) next
        key <- paste(round(m[, 1] / tol), round(m[, 2] / tol), sep = "_")
        key <- key[-length(key)]
        if (length(key) < 1) next
        a <- key; b <- c(key[-1], key[1])
        if (identical(type, "queen")) { f <- c(f, rep.int(u, length(a))); ka <- c(ka, a) }
        else { lo <- ifelse(a < b, a, b); hi <- ifelse(a < b, b, a); f <- c(f, rep.int(u, length(a))); ka <- c(ka, paste(lo, hi, sep = "|")) }
    }
    d <- unique(data.frame(f = f, k = ka, stringsAsFactors = FALSE))
    d <- d[d$k %in% d$k[duplicated(d$k)], , drop = FALSE]
    if (nrow(d) == 0) return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n)))
    pr <- merge(d, d, by = "k")
    pr <- pr[pr$f.x != pr$f.y, , drop = FALSE]
    A <- Matrix::sparseMatrix(i = pr$f.x, j = pr$f.y, x = 1, dims = c(n, n))
    A@x[] <- 1
    A
}

# ---------------------------------------------------------- building weights --

# File paths: surrounding quotes are removed, %VARIABLE% is expanded, and "sample:name" points to the example files installed with the module.
sp_resolve_path <- function(p) {
    p <- trimws(gsub('^"|"$', "", trimws(p)))
    if (grepl("^sample:", p, ignore.case = TRUE)) {
        f <- sub("^sample:[[:space:]]*", "", p, ignore.case = TRUE)
        q <- system.file("data", f, package = "spatialfamily")
        return(if (nzchar(q)) q else file.path("data", f))
    }
    m <- regmatches(p, gregexpr("%[A-Za-z_]+%", p))[[1]]
    for (v in unique(m)) p <- gsub(v, Sys.getenv(substr(v, 2, nchar(v) - 1)), p, fixed = TRUE)
    p
}

# Reads the weights options, builds the weights for the rows of `df` (the estimation sample).
# `s` is a plain list: source, file, idfield, method, k, dband, power, cutoff, order, cumul, symm, lonlat, xvar, yvar, idvar.
sp_build_weights <- function(df, s) {
    n <- nrow(df)
    src <- s$source
    desc <- ""; xy <- NULL; geom <- NULL; lonlat <- isTRUE(s$lonlat); ids <- as.character(seq_len(n))
    if (src == "coords") {
        if (is.null(s$xvar) || is.null(s$yvar)) stop("Select the two coordinate variables (x / longitude and y / latitude).", call. = FALSE)
        xy <- cbind(sp_num(df[[s$xvar]]), sp_num(df[[s$yvar]]))
        if (!is.null(s$idvar)) ids <- as.character(df[[s$idvar]])
        if (!s$method %in% c("knn", "dist", "idw")) stop("With coordinates choose nearest neighbours, distance band or inverse distance.", call. = FALSE)
        if (lonlat && (any(abs(xy[, 1]) > 180) || any(abs(xy[, 2]) > 90)))
            stop("The coordinates are not valid longitude / latitude degrees; untick the longitude/latitude option for projected coordinates.", call. = FALSE)
        L <- sp_links_coords(xy, s$method, s$k, s$dband, s$power, s$cutoff, lonlat)
        B <- L$B; desc <- paste0("coordinates: ", L$desc)
    } else {
        if (is.null(s$file) || !nzchar(trimws(s$file))) stop("Enter the path of the file.", call. = FALSE)
        path <- sp_resolve_path(s$file)
        if (is.null(s$idvar)) stop("Select the unit id variable that links the data to the file.", call. = FALSE)
        ids <- as.character(df[[s$idvar]])
        if (anyDuplicated(ids)) stop(sprintf("The id variable has repeated values (e.g. '%s'); the weights need one row per unit.", ids[anyDuplicated(ids)]), call. = FALSE)
        if (src == "boundary") {
            b <- sp_read_boundary(path)
            fid <- if (nzchar(trimws(s$idfield))) {
                if (is.null(b$attrs)) stop("The file has no attribute table to read the id field from.", call. = FALSE)
                if (!trimws(s$idfield) %in% names(b$attrs)) stop(sprintf("The id field '%s' is not in the file. Fields: %s.", trimws(s$idfield), paste(names(b$attrs), collapse = ", ")), call. = FALSE)
                as.character(b$attrs[[trimws(s$idfield)]])
            } else {
                if (length(b$units) != n) stop(sprintf("Without an id field the file's %d units must match the %d data rows in the same order; enter the id field name instead.", length(b$units), n), call. = FALSE)
                ids
            }
            if (anyDuplicated(fid)) stop("The id field in the file has repeated values.", call. = FALSE)
            m <- match(ids, fid)
            if (anyNA(m)) stop(sprintf("%d unit(s) in the data have no match in the file (e.g. '%s'). Check the id variable and the id field (text and numbers must look the same).", sum(is.na(m)), ids[which(is.na(m))[1]]), call. = FALSE)
            units <- b$units[m]
            cent <- sp_centroids(units)
            if (s$method %in% c("queen", "rook")) {
                B <- sp_contig(units, s$method); desc <- sprintf("%s contiguity from the boundary file", s$method)
            } else {
                L <- sp_links_coords(cent, s$method, s$k, s$dband, s$power, s$cutoff, lonlat)
                B <- L$B; desc <- paste0("polygon centroids: ", L$desc)
            }
            xy <- cent; geom <- units
        } else if (src == "adjacency") {
            e <- sp_read_edges(path)
            fi <- match(e$from, ids); ti <- match(e$to, ids)
            if (all(is.na(fi) | is.na(ti))) stop("None of the ids in the file match the id variable in the data. Check that they use the same codes.", call. = FALSE)
            ok <- !is.na(fi) & !is.na(ti)
            B <- Matrix::sparseMatrix(i = fi[ok], j = ti[ok], x = e$w[ok], dims = c(n, n))
            desc <- sprintf("neighbour list from file (%d links read, %d with both units in the data)", nrow(e), sum(ok))
        } else stop("Unknown weights source.", call. = FALSE)
    }
    B <- as(B, "CsparseMatrix"); Matrix::diag(B) <- 0; B <- Matrix::drop0(B)
    if (isTRUE(s$symm)) B <- sp_symmetrize(B)
    if (s$order > 1 && (src != "coords" || s$method %in% c("knn", "dist"))) {
        B <- sp_order(B, s$order, isTRUE(s$cumul))
        desc <- sprintf("%s; order %d%s", desc, s$order, if (isTRUE(s$cumul)) " (cumulative)" else "")
    }
    if (length(B@x) == 0) stop("The weights have no links at all; check the settings (for example a larger distance band).", call. = FALSE)
    list(B = B, ids = ids, xy = xy, geom = geom, lonlat = lonlat, desc = desc)
}

# Weights as used by the statistics for the chosen style.
sp_W <- function(w, style = "W") sp_style(w$B, style)
