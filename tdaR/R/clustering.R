
# The graph-summary row TDA prints while building a valued graph from a
# dissimilarity matrix (mds and cluster share it); stored so the audit
# matches it by content.
.graph_summary <- function(out, res = NULL) {
    nm <- c("graph", "edges", "loops", "max_edge_value", "isolated_nodes")
    ex <- if (!is.null(res) && .use_exports())
        res$exports[["graph.summary"]]
    if (is.matrix(ex) && ncol(ex) == 5L && nrow(ex) >= 1L)
        return(stats::setNames(as.list(ex[1L, ]), nm))
    gi <- grep("^Graph\\s+edges\\s+loops", out)
    if (!length(gi))
        return(NULL)
    gv <- suppressWarnings(as.numeric(
        strsplit(trimws(out[gi[1L] + 1L]), "\\s+")[[1L]]))
    if (length(gv) == 5L && !anyNA(gv))
        stats::setNames(as.list(gv), nm)
}
# clustering, scaling, dendrogram cuts, and conjoint analysis
#


# ---- clustering and scaling, from a dist -----------------------------------

# gdd option 4 reads the upper triangle of a dissimilarity matrix from a single
# variable.  R's dist stores the lower triangle column-wise, which for a
# symmetric matrix is the same sequence, so as.vector() is already in the order
# TDA wants.
.dist_graph <- function(d) {
    if (inherits(d, "dist")) {
        v <- as.vector(d)
        n <- attr(d, "Size")
        lab <- attr(d, "Labels")
    } else {
        m <- as.matrix(d)
        if (nrow(m) != ncol(m))
            stop("`d` must be a dist or a square dissimilarity matrix")
        n <- nrow(m)
        v <- m[lower.tri(m)]
        lab <- rownames(m)
    }
    if (n < 3L)
        stop("at least three objects are needed")
    list(data = data.frame(D = as.numeric(v)), n = n,
         labels = lab %||% as.character(seq_len(n)))
}

.dist_cmd <- function(cmd, d, opts, dir, cls, ...) {
    # one validation point covers every distance-based wrapper
    d <- .tda_check_matrix(d, "d", square = TRUE, symmetric = TRUE)
    g <- .dist_graph(d)
    res <- tda_run(c(tda_nvar(g$data), "gdd(opt=4) = D;",
                     do.call(tda_block, c(list(name = cmd), opts,
                             list(rhs = "out.txt")))),
                   data = g$data, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run ", cmd, ": ", err[1L], call. = FALSE)
    # hcld (the default, "hierarchical") writes a ragged table -- node,
    # member count, then that many member indices -- which read.table()
    # cannot parse (every row a different width) and silently drops via the
    # tryCatch below. Parsed separately into one list entry per cluster.
    clusters <- if (identical(cmd, "hcld"))
        .parse_clusters(res, "out.txt") else NULL
    structure(list(call = sys.call(-1L), run = res, n = g$n,
                   graph = .graph_summary(res$output, res),
                   labels = g$labels,
                   table = if (!identical(cmd, "hcld"))
                       tryCatch(tda_file(res, "out.txt"),
                                error = function(e) NULL),
                   clusters = clusters),
              class = c(cls, "tda_table"))
}

# Reads hcld's "node  count  member1 member2 ..." rows -- ragged, so not a
# data.frame -- into a named list of member-index vectors, one per node.
.parse_clusters <- function(res, file) {
    rows <- .file_rows(res, file)
    if (is.null(rows)) {
        p <- file.path(res$dir, file)
        if (!file.exists(p))
            return(NULL)
        lines <- readLines(p, warn = FALSE)
        lines <- lines[nzchar(trimws(lines))]
        if (!length(lines))
            return(NULL)
        rows <- lapply(lines, function(l)
            suppressWarnings(as.integer(strsplit(trimws(l), "\\s+")[[1L]])))
    }
    if (!length(rows))
        return(NULL)
    out <- lapply(rows, function(tok) {
        tok <- as.integer(tok)
        if (length(tok) < 2L || is.na(tok[1L]) || is.na(tok[2L]))
            return(NULL)
        tok[seq_len(min(tok[2L], length(tok) - 2L)) + 2L]
    })
    nm <- vapply(rows, function(tok) as.character(as.integer(tok[1L])),
                 character(1))
    names(out) <- nm
    out[!vapply(out, is.null, logical(1))]
}

#' @rdname tda_cluster
#' @export
TDA_CLUSTER <- c(hierarchical = "hcld", hierarchical_single = "hcls",
                 unidimensional = "ucl", unidimensional_scaling = "uds",
                 blockmodel = "becl", nearest_neighbour = "nncl",
                 additive = "acl", partition = "clu", partition2 = "clp")

#' Clustering from a dissimilarity matrix
#'
#' TDA's clustering commands, all of which take a dissimilarity and so accept
#' a \code{dist} -- including the one \code{\link{tda_seqm}} returns.
#'
#' The default, \code{"hierarchical"} (\code{hcld}), writes a ragged table --
#' a node, its member count, then that many member indices -- which is not a
#' rectangular \code{data.frame}, so it is parsed separately into
#' \code{clusters}, a named list of member-index vectors, one per node of the
#' divisive tree (the root first, then each split). \code{"unidimensional"}
#' and \code{"additive"} report their result only as text (a final partition
#' or a fitted tree), never through \code{df=}; read \code{$run$output} for
#' those.
#'
#' @param d a \code{dist}, or a square dissimilarity matrix.
#' @param method one of the names in \code{TDA_CLUSTER}, matched partially.
#' @param algorithm for \code{method = "hierarchical"} only: split by
#'   \code{"centers"} (default, maximally different cluster centers) or
#'   \code{"diameter"} (minimal-diameter partitions) -- \code{hcld}'s
#'   \code{alg=}.
#' @param min_size for \code{algorithm = "centers"}, the minimum cluster
#'   size (\code{hcld}'s \code{min=}).
#' @param max_splits for \code{algorithm = "diameter"}, the maximum
#'   number of splits (\code{hcld}'s \code{max=}).
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object carrying a \code{table}, or for \code{"hierarchical"}
#'   \code{clusters} instead (see Details). For \code{"hierarchical"} the
#'   object also carries \code{merges}: one row per split, the two
#'   clusters, their sizes and their diameters.
#' @family clustering
#' @examples
#' d <- dist(matrix(c(1, 2, 3, 8, 9, 10), ncol = 1))
#' tda_cluster(d)$clusters  # the divisive tree: root, then each split
#' @export
tda_cluster <- function(d, method = "hierarchical",
                        algorithm = c("centers", "diameter"),
                        min_size = NULL, max_splits = NULL,
                        options = list(), dir = tempfile("tda"), ...) {
    i <- pmatch(method, names(TDA_CLUSTER))
    if (is.na(i))
        stop("unknown method '", method, "'; one of: ",
             paste(names(TDA_CLUSTER), collapse = ", "))
    opts <- .tda_extra(options)
    # alg=/min=/max= are hcld's (confirmed against its help text);
    # only meaningful for method = "hierarchical", which is what hcld
    # implements -- the other TDA_CLUSTER methods use a different command
    # and do not take these.
    if (!missing(algorithm) || !is.null(min_size) || !is.null(max_splits)) {
        if (!identical(unname(TDA_CLUSTER[i]), "hcld"))
            stop("`algorithm`/`min_size`/`max_splits` are hcld's ",
                 "options, only meaningful for method = \"hierarchical\"")
        if (!missing(algorithm))
            opts$alg <- match(match.arg(algorithm), c("centers", "diameter"))
        if (!is.null(min_size)) opts$min <- min_size
        if (!is.null(max_splits)) opts$max <- max_splits
        # hcld writes a second file with df=: one row per split, the two
        # clusters, their sizes and their diameters (manual 7.5.2.1 and
        # 7.5.2.2, each section's Box 2).  Without df= it is not written.
        if (is.null(opts$df))
            opts$df <- "out.df"
    }
    r <- .dist_cmd(unname(TDA_CLUSTER[i]), d, opts, dir, "tda_cluster", ...)
    if (!is.null(opts$df) && !is.null(r$run$dir)) {
        p <- file.path(r$run$dir, opts$df)
        if (file.exists(p)) {
            m <- tryCatch(utils::read.table(p, header = FALSE),
                          error = function(e) NULL)
            if (is.data.frame(m) && ncol(m) == 6L)
                m <- .name_cols(m, c("ci", "cj", "ni", "nj",
                                     "diameter_i", "diameter_j"))
            r$merges <- m
        }
    }
    r
}

#' @rdname tda_mds
#' @export
TDA_MDS <- c(classical = "mdsc", metric = "mdsm", nonmetric = "mdsn",
             minmax = "mdsx")

#' Multidimensional scaling
#'
#' Represent a dissimilarity matrix as points in a low-dimensional space,
#' so that the distances between the points approximate the
#' dissimilarities.
#'
#' \code{method = "classical"} (TDA's \code{mdsc}) is the eigenvalue
#' decomposition also known as principal coordinates analysis -- the same
#' method as \code{\link[stats]{cmdscale}}. TDA computes the full
#' solution, one dimension per positive eigenvalue, and never reads a
#' dimension count of its own; because the classical solution is nested,
#' the first \code{ndim} columns of the full solution \emph{are} the
#' \code{ndim}-dimensional solution, so this wrapper simply keeps those.
#' The eigenvalues, with the share of the total each accounts for, come
#' back in \code{eigenvalues} -- the usual guide to how many dimensions
#' the data support.
#'
#' The other three methods fit the configuration iteratively.
#' \code{"metric"} (\code{mdsm}) minimizes the raw stress, the sum of
#' squared differences between dissimilarities and fitted distances,
#' with TDA's general minimizer; it fits two dimensions only, and errors
#' for any other \code{ndim}. \code{"nonmetric"} (\code{mdsn}) minimizes
#' Kruskal's stress-1, preserving only the rank order of the
#' dissimilarities, by Kruskal's (1964) gradient method with his step
#' length adaptation; it reaches the same or lower stress than
#' \code{\link[MASS]{isoMDS}} on the usual test data, but its default
#' \code{mxit = 100} iterations is often not enough -- pass
#' \code{options = list(mxit = 300)}, and \code{ns =} random restarts to
#' guard against local minima. \code{"minmax"} (\code{mdsx}) fits a
#' minmax criterion; TDA labels that command experimental and its output
#' is left as it comes. Check \code{stress} before trusting any
#' iterative configuration.
#'
#' @param d a \code{dist}, or a square dissimilarity matrix.
#' @param method one of the names in \code{TDA_MDS}, matched partially.
#' @param ndim number of dimensions to return.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object carrying \code{points}, the coordinates (one row per
#'   object, one column per dimension), and, for \code{"classical"},
#'   \code{eigenvalues} (a data frame of each eigenvalue and its
#'   percentage of the total); for the iterative methods, \code{stress},
#'   the criterion value at the returned configuration.
#' @family clustering
#' @examples
#' set.seed(34)
#' pts <- rbind(matrix(rnorm(10, 0), 5, 2), matrix(rnorm(10, 6), 5, 2))
#' m <- tda_mds(dist(pts), method = "classical", ndim = 2)
#' m$points       # recovers the two clusters as two groups of coordinates
#' m$eigenvalues  # two dominant dimensions, as built
#' @export
tda_mds <- function(d, method = "classical", ndim = 2, options = list(),
                    dir = tempfile("tda"), ...) {
    i <- pmatch(method, names(TDA_MDS))
    if (is.na(i))
        stop("unknown method '", method, "'; one of: ",
             paste(names(TDA_MDS), collapse = ", "))
    cmd <- unname(TDA_MDS[i])
    opts <- .tda_extra(options)
    # mdsc never reads ndim= at all -- confirmed in t_mds.c: it writes
    # one coordinate column per positive eigenvalue, full stop -- so the
    # dimension count is applied here instead (see Details for why that
    # is exact for the classical method); the iterative commands do read
    # it and fit that many dimensions directly.
    # mdsm parametrizes two dimensions only (t_mds.c: NParm = 2 (n - 1))
    # and reads no ndim=; refuse rather than return two columns for three
    if (cmd == "mdsm" && !is.null(ndim) && ndim != 2)
        stop("method = \"metric\" (mdsm) fits two dimensions only",
             call. = FALSE)
    if (!is.null(ndim) && cmd != "mdsc")
        opts$ndim <- ndim
    # the coordinates are written with fmt=, default 10.4; full precision
    # costs nothing and spares round-off against cmdscale()
    if (is.null(opts$fmt))
        opts$fmt <- "24.16"
    r <- .dist_cmd(cmd, d, opts, dir, "tda_mds", ...)
    out <- r$run$output
    if (cmd == "mdsc") {
        pts <- r$table
        if (is.data.frame(pts) && !is.null(ndim) && ncol(pts) > ndim)
            pts <- pts[, seq_len(ndim), drop = FALSE]
        r$eigenvalues <- .mds_eigenvalues(out, r$run)
        # the graph-summary row (edges/loops/max edge value/isolated
        # nodes) TDA prints before the eigenvalues, stored so the audit
        # matches it by content
        r$graph <- .graph_summary(out, r$run)
    } else {
        # mdsm and mdsn print the fitted coordinates on the console and
        # use the output file for other things (per-repeat parameter
        # rows; the fitted distance matrix), so the console block is the
        # coordinates; mdsx does write coordinates to the file, first
        # column the node number.
        pts <- .mds_console_points(out)
        if (is.null(pts) && is.data.frame(r$table)) {
            pts <- r$table[, -1L, drop = FALSE]
        }
        st <- regmatches(out, regexpr(
            "(?<=Final best stress value:)\\s*[0-9.eE+-]+", out, perl = TRUE))
        if (!length(st))
            st <- regmatches(out, regexpr("(?<=\\(stress:)\\s*[0-9.eE+-]+",
                                          out, perl = TRUE))
        if (length(st))
            r$stress <- suppressWarnings(as.numeric(trimws(st[length(st)])))
    }
    if (is.data.frame(pts)) {
        names(pts) <- paste0("dim", seq_len(ncol(pts)))
        if (!is.null(r$labels) && length(r$labels) == nrow(pts))
            rownames(pts) <- r$labels
    }
    r$points <- pts
    r$method <- names(TDA_MDS)[i]
    class(r) <- c("tda_mds", class(r))
    r
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::print
print.tda_mds <- function(x, ...) {
    cat("Multidimensional scaling (", x$method, "), ", x$n, " objects",
        sep = "")
    if (is.data.frame(x$points))
        cat(", ", ncol(x$points), " dimension",
            if (ncol(x$points) != 1L) "s", sep = "")
    cat("\n")
    if (!is.null(x$stress))
        cat("Stress:", format(x$stress, digits = 6), "\n")
    if (is.data.frame(x$eigenvalues)) {
        pc <- x$eigenvalues$percent
        k <- if (is.data.frame(x$points)) ncol(x$points) else 2L
        if (any(is.finite(pc)))
            cat("Variance accounted for by the ", k, " dimension",
                if (k != 1L) "s", " kept: ",
                format(sum(pc[seq_len(min(k, length(pc)))], na.rm = TRUE),
                       digits = 4),
                "% (eigenvalues in $eigenvalues)\n", sep = "")
    }
    if (is.data.frame(x$points)) {
        cat("\n")
        print(utils::head(x$points, 6L))
        if (nrow(x$points) > 6L)
            cat("  ...", nrow(x$points) - 6L, "more rows\n")
    }
    invisible(x)
}

# The "Coordinates" block mdsm/mdsn print: one row per point, the node
# number first, then the fitted coordinates.
.mds_console_points <- function(out) {
    i <- grep("^Coordinates", out)
    if (!length(i))
        return(NULL)
    rows <- list()
    for (l in out[seq.int(i[length(i)] + 1L, length(out))]) {
        v <- suppressWarnings(as.numeric(strsplit(trimws(l), "\\s+")[[1L]]))
        if (!length(v) || anyNA(v))
            break
        rows[[length(rows) + 1L]] <- v[-1L]
    }
    if (!length(rows))
        return(NULL)
    as.data.frame(do.call(rbind, rows))
}

# The eigenvalue table mdsc prints: each eigenvalue and, when all are
# non-negative, its percentage of the total.
.mds_eigenvalues <- function(out, res = NULL) {
    i <- grep("^Eigenvalue\\s+per cent", out)
    if (!length(i))
        return(NULL)
    rows <- list()
    for (l in out[seq.int(i[1L] + 1L, length(out))]) {
        v <- suppressWarnings(as.numeric(strsplit(trimws(l), "\\s+")[[1L]]))
        if (!length(v) || anyNA(v))
            break
        rows[[length(rows) + 1L]] <- c(v, NA_real_)[1:2]
    }
    if (!length(rows))
        return(NULL)
    m <- as.data.frame(do.call(rbind, rows))
    names(m) <- c("value", "percent")
    if (!is.null(res))
        m <- .overlay_cols(m, res$exports[["mds.eigenvalues"]],
                           c("value", "percent"))
    m
}


#' Cut a dendrogram into groups
#'
#' Builds a dendrogram with TDA's \code{hcls} and cuts it with \code{hclsp}.
#'
#' Unlike R's \code{cutree}, which returns one group label per observation,
#' \code{hclsp} reports the cut as a set of nodes: \code{table} has one row
#' per resulting cluster, its TDA node number, how many leaves (original
#' observations) it holds, and its parent node -- not a per-observation
#' vector. Getting a label for each observation means walking \code{merge}
#' (the dendrogram \code{hcls} built, one row per step) down from a node in
#' \code{table} to the leaves under it; this wrapper does not do that walk
#' for you.
#'
#' @param d a \code{dist}, or a square dissimilarity matrix.
#' @param nlev number of levels at which to cut.
#' @param nodes an explicit sequence of node numbers to cut at, TDA's
#'   \code{cn=}.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object carrying \code{merge}, the dendrogram, and \code{table},
#'   one row per cluster in the cut (\code{node}, \code{leaves},
#'   \code{parent}).
#' @family clustering
#' @examples
#' set.seed(35)
#' pts <- rbind(matrix(rnorm(10, 0), 5, 2), matrix(rnorm(10, 8), 5, 2))
#' ct <- tda_cutree(dist(pts), nlev = 2)
#' ct$table    # two clusters of five leaves each
#' ct$merge    # the nine merge steps that built the dendrogram
#' @export
tda_cutree <- function(d, nlev = 2, nodes = NULL, options = list(),
                       dir = tempfile("tda"), ...) {
    g <- .dist_graph(d)
    n <- g$n

    # Step one: build the dendrogram.  hcls with opt=2 writes one row per
    # merge -- index, height, and the two clusters joined.
    r1 <- tda_run(c(tda_nvar(g$data), "gdd(opt=4) = D;",
                    "hcls(opt=2) = dend.txt;"), data = g$data, dir = dir, ...)
    err <- grep("^Error", r1$output, value = TRUE)
    if (length(err))
        stop("TDA could not build the dendrogram: ", err[1L], call. = FALSE)
    ed0 <- if (.use_exports()) r1$exports[["hcls.dend"]]
    dend <- if (is.matrix(ed0)) .export_frame(ed0)
            else tda_file(r1, "dend.txt")
    if (is.null(dend) || ncol(dend) < 4L)
        stop("hcls wrote no dendrogram")
    names(dend)[1:4] <- c("step", "height", "a", "b")

    # Step two: hclsp wants that tree as an edge list, and TDA identifies a
    # merged cluster by one of its members rather than by a fresh number, so
    # the current top node of each cluster has to be tracked.  Edges point
    # from child to parent, which is what leaves the root with outdegree 0.
    top <- seq_len(n)
    edges <- matrix(0L, 0L, 2L)
    for (k in seq_len(nrow(dend))) {
        a <- dend$a[k]
        b <- dend$b[k]
        parent <- n + k
        # Capture the two cluster tops before either is overwritten: setting
        # top[a] first makes the propagation below a no-op, so the rest of
        # each cluster never moves under the new parent and the result is not
        # a tree.
        ta <- top[a]
        tb <- top[b]
        edges <- rbind(edges, c(ta, parent), c(tb, parent))
        top[top == ta | top == tb] <- parent
    }
    ed <- data.frame(I = as.integer(edges[, 1L]), J = as.integer(edges[, 2L]),
                     V = 1L)

    opts <- .tda_extra(options)
    if (!is.null(nodes))
        opts$cn <- paste(nodes, collapse = ",")
    else
        opts$nlev <- nlev
    r2 <- tda_run(c(tda_nvar(ed), "gdd(opt=1) = I,J,V;",
                    do.call(tda_block, c(list(name = "hclsp"), opts,
                            list(rhs = "out.txt")))),
                  data = ed, dir = dir, ...)
    err <- grep("^Error", r2$output, value = TRUE)
    if (length(err))
        stop("TDA could not cut this dendrogram: ", err[1L], call. = FALSE)

    # hclsp describes the partition in words rather than writing a table:
    # "(39,20) <- (37,10)" is the root with twenty leaves splitting into a
    # node with ten.  The (node, size) pairs are pulled out of that.
    rows <- .file_rows(r2, "out.txt")
    flat <- if (!is.null(rows)) {
        # every line is "(node,size) <- (node,size) ...": the numbers,
        # in order, are those pairs
        unlist(lapply(rows, function(v) sprintf("(%d,%d)",
                      as.integer(v[c(TRUE, FALSE)]), as.integer(v[c(FALSE, TRUE)]))))
    } else {
        txt <- tryCatch(readLines(file.path(r2$dir, "out.txt"), warn = FALSE),
                        error = function(e) character())
        unlist(regmatches(txt, gregexpr("\\(\\d+,\\d+\\)", txt)))
    }
    tab <- NULL
    if (length(flat)) {
        v <- do.call(rbind, lapply(flat, function(z)
            as.integer(strsplit(gsub("[()]", "", z), ",")[[1L]])))
        tab <- data.frame(node = v[, 1L], leaves = v[, 2L])
        tab$parent <- c(NA_integer_, rep(tab$node[1L], nrow(tab) - 1L))
    }

    structure(list(call = match.call(), run = r2, dendrogram = r1,
                   graph = .graph_summary(c(r1$output, r2$output), r1),
                   n = n, labels = g$labels, merge = dend, table = tab),
              class = c("tda_cutree", "tda_table"))
}


#' Non-metric conjoint analysis
#'
#' \code{nmca} is conjoint measurement: it takes a rank order over profiles
#' and the categorical factors describing them, and estimates the part-worth
#' of each factor level. It draws a Hasse diagram of the order as part of the
#' output.
#'
#' @param formula \code{rank ~ f1 + f2}, where the response is the rank order
#'   and the right-hand side names the factors.
#' @param data a data frame. Factor levels are numbered from 1 upwards, and a
#'   value outside that range is rejected by TDA with the record it occurred
#'   in.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object carrying the run; the diagram and the estimates are in
#'   \code{$run$output}.
#' @family clustering
#' @examples
#' d <- expand.grid(price = 1:2, brand = 1:2)
#' d$rank <- c(4, 2, 3, 1)
#' cj <- tda_conjoint(rank ~ price + brand, d)
#' cat(tda_payload(cj$run), sep = "\n")   # the nmca section itself
#' @export
tda_conjoint <- function(formula, data, options = list(),
                         dir = tempfile("tda"), ...) {
    p <- .reg_parts(formula, data)
    # Every factor has to be coded 1..k with no gaps, which is what TDA checks
    # record by record.
    X <- as.data.frame(p$X)
    for (j in seq_along(X))
        X[[j]] <- as.integer(factor(X[[j]]))
    d <- cbind(stats::setNames(data.frame(as.numeric(p$y)), .tda_names(p$ylab)),
               stats::setNames(X, .tda_names(p$xlab)))
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = "nmca"), .tda_extra(options),
                             list(rhs = paste(names(d), collapse = ","))))),
                   data = d, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run the analysis: ", err[1L], call. = FALSE)
    structure(list(call = match.call(), run = res, n = nrow(d),
                   xlab = c(p$ylab, p$xlab), value = NA_real_),
              class = c("tda_conjoint", "tda_expr"))
}


#' Principal components, factor analysis, correspondence analysis
#'
#' \code{dma}: seven related methods for reducing a data matrix to a
#' handful of derived scores, selected with \code{alg} -- principal
#' components (1, from the raw data; 2, from a covariance or correlation
#' matrix already computed), factor analysis (3, 4, the same two starting
#' points but with eigenvectors rescaled as loadings), dual scaling and
#' correspondence analysis (5, 6, for a frequency table rather than a
#' data matrix), and a direct SVD-based projection (7). Verified against
#' \code{\link[stats]{prcomp}} for algorithm 1: eigenvalues match exactly,
#' eigenvectors match up to sign (an ordinary, harmless ambiguity in any
#' eigendecomposition).
#'
#' TDA's \code{v=} (choosing a subset of variables from a larger data
#' matrix) needs variable names that start uppercase, the same rule as
#' anywhere else in TDA -- lowercase R column names are displayed
#' correctly when declared but cannot be \emph{referenced} again later in
#' the same script, including by \code{v=} itself, failing with a
#' \dQuote{Syntax error or undefined variables} that gives no hint the
#' name's case is the problem. \code{variables} handles this
#' automatically; there is no need to rename columns yourself first.
#'
#' @param data a data frame or matrix. For algorithms 1--4 and 7, a
#'   standard data matrix (rows are cases); for 2 and 4, a covariance or
#'   correlation matrix (square, symmetric); for 5 and 6, a frequency
#'   table.
#' @param variables optional subset of \code{data}'s column names to
#'   use, instead of all of them.
#' @param alg which of the seven algorithms; see Description. Default 1,
#'   ordinary principal components.
#' @param preprocess preprocessing, for algorithms 1, 5 and 7 only -- for
#'   \code{alg = 1}: \code{1} nothing (TDA's default), \code{2} mean
#'   centre, \code{3} standardize; for \code{alg} 5 or 7: \code{2}/\code{3}
#'   convert to row or table relative frequencies.
#' @param ns keep at most this many components; all of them by default.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#'   Further \code{dma} options can be given the same way, notably
#'   \code{opt=} (algorithm variant) and \code{pcf=} (print
#'   classification frequencies).
#' @return A list with \code{values} (the eigenvalues), \code{percent} (each
#'   eigenvalue's share of the total, as TDA prints it), \code{vectors}
#'   (the eigenvectors or factor loadings, one column per component, row
#'   names taken from \code{data}), and \code{scores} (the derived score
#'   for each case, one column per component) -- \code{scores} is empty
#'   for algorithms where TDA does not calculate it (2 and 4).
#' @family clustering
#' @examples
#' set.seed(1)
#' d <- data.frame(x1 = rnorm(30), x2 = rnorm(30))
#' d$x3 <- d$x1 + d$x2 + rnorm(30, sd = 0.3)   # correlated with both
#' pc <- tda_dma(d)
#' pc$values                 # most of the variance is in the first PC
#' pc$vectors                # loadings: x3 dominates it
#' head(pc$scores)
#'
#' # a subset of the variables, instead of all of them
#' tda_dma(d, variables = c("x1", "x2"))$values
#' @export
tda_dma <- function(data, alg = 1, preprocess = NULL, ns = NULL, variables = NULL,
                    options = list(), dir = tempfile("tda"), ...) {
    d <- as.data.frame(data)
    orig_names <- names(d)
    names(d) <- .tda_names(orig_names)
    if (length(alg) != 1L || is.na(alg) || !alg %in% 1:7)
        stop("`alg` must be one of 1..7; TDA itself would silently run ",
             "the default algorithm instead")
    opts <- list(alg = alg, df = "eigvec.out", fmt = "18.12")
    if (!is.null(preprocess))
        opts$opt <- preprocess
    if (!is.null(ns))
        opts$ns <- ns
    if (!is.null(variables)) {
        vi <- match(variables, orig_names)
        if (anyNA(vi))
            stop("no such variable: ", paste(variables[is.na(vi)],
                                             collapse = ", "))
        opts$v <- paste(names(d)[vi], collapse = ",")
    }
    opts <- c(opts, .tda_extra(options))
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = "dma"), opts,
                                          list(rhs = "scores.out")))),
                   data = d, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run dma: ", err[1L], call. = FALSE)
    i <- grep("^Eigenvalue", res$output)
    ev <- pc <- numeric()
    if (length(i)) {
        rest <- trimws(res$output[(i[1L] + 1L):length(res$output)])
        rest <- rest[nzchar(rest)]
        keep <- character()
        for (r in rest) {
            if (!grepl("^-?[0-9.]+\\s", r))
                break
            keep <- c(keep, r)
        }
        ev <- as.numeric(sub("\\s+.*", "", keep))
        pc <- suppressWarnings(as.numeric(sub("^\\S+\\s+", "", keep)))
    }
    used_names <- if (!is.null(variables)) variables else orig_names
    # export first for both matrices: dma.vectors comes from dma_prn3()
    # (the eigvec.out writer) and dma.scores from the score loop, so
    # neither file is read when the producers cover the run.
    ev0 <- if (.use_exports()) res$exports[["dma.vectors"]]
    sc0 <- if (.use_exports()) res$exports[["dma.scores"]]
    vec <- if (is.matrix(ev0)) as.data.frame(ev0)
           else tda_file(res, "eigvec.out")
    sco <- if (is.matrix(sc0)) as.data.frame(sc0)
           else tda_file(res, "scores.out")
    if (!is.null(vec)) {
        vec <- as.matrix(vec)
        colnames(vec) <- paste0("PC", seq_len(ncol(vec)))
        if (ncol(vec) && nrow(vec) == length(used_names))
            rownames(vec) <- used_names
    }
    if (!is.null(sco)) {
        sco <- as.matrix(sco)
        colnames(sco) <- paste0("PC", seq_len(ncol(sco)))
    }
    # the per cent column at %.2f: TDA's is the share of the total, so
    # it is recomputed from the doubles rather than read from the text
    if (length(ev) && all(is.finite(ev)) && sum(ev) != 0)
        pc <- 100 * ev / sum(ev)
    structure(list(call = match.call(), run = res, values = ev,
                   percent = pc, vectors = vec, scores = sco),
              class = "tda_dma")
}

#' Distance matrix from raw variables
#'
#' \code{pdatd}: builds a case-by-case distance matrix directly from a
#' data matrix, the way \code{\link[stats]{dist}} does -- Euclidean or
#' city-block distance, the count of variables on which two cases differ,
#' or a dissimilarity index. The result feeds directly into
#' \code{\link{tda_cluster}} or \code{\link{tda_mds}}. Verified against
#' \code{\link[stats]{dist}} for Euclidean and city-block (Manhattan)
#' distance: exact match.
#'
#' As with \code{\link{tda_dma}}, TDA's \code{v=} needs variable names
#' that start uppercase to be referenced correctly; \code{variables}
#' handles this automatically.
#'
#' Per-variable weights work through \code{weights}, with one genuine
#' TDA quirk absorbed for you: \code{wt=} is parsed by TDA's
#' \emph{time-points} parser, which accepts only a strictly increasing
#' positive list -- \code{wt=2,1,1} is a syntax error while
#' \code{wt=1,2,4} is fine, which is why this option long looked dead.
#' Since a distance does not care about variable order, the wrapper
#' sends the variables sorted by ascending weight. The weighting
#' convention, verified against R directly: euclidean uses
#' \eqn{d^2 = \sum_j w_j \Delta_j^2} (so it matches
#' \code{dist(sweep(x, 2, sqrt(w), "*"))}), city-block
#' \eqn{d = \sum_j w_j |\Delta_j|}. \emph{Tied} weights cannot pass
#' TDA's parser at all; scale those columns yourself
#' (\eqn{\sqrt{w} x} for euclidean, \eqn{w x} for city-block) and
#' call this unweighted -- the result is identical.
#'
#' @param data a data frame or matrix, one row per case.
#' @param weights optional positive per-variable weights, one per
#'   (selected) variable, all distinct -- see Details for TDA's parser
#'   constraint and the exact weighting convention.
#' @param measure the distance measure: \code{1} Euclidean (default), \code{2}
#'   city-block, \code{3} number of variables on which the two cases
#'   differ, \code{4} a dissimilarity index.
#' @param variables optional subset of \code{data}'s column names to
#'   use, instead of all of them.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#'   \code{opt=} (the command's variant switch) passes through as
#'   well.
#' @return A \code{\link[stats]{dist}} object.
#' @family clustering
#' @examples
#' set.seed(1)
#' d <- data.frame(x1 = rnorm(5), x2 = rnorm(5))
#' tda_pdatd(d)
#' tda_pdatd(d, measure = 2)   # city-block, not the euclidean default
#' tda_pdatd(d, variables = "x1")   # a subset of the variables
#' @export
tda_pdatd <- function(data, measure = 1, variables = NULL, weights = NULL,
                      options = list(), dir = tempfile("tda"), ...) {
    d <- as.data.frame(data)
    orig_names <- names(d)
    names(d) <- .tda_names(orig_names)
    opts <- list(opt = measure, fmt = "18.12")
    if (!is.null(variables)) {
        vi <- match(variables, orig_names)
        if (anyNA(vi))
            stop("no such variable: ", paste(variables[is.na(vi)],
                                             collapse = ", "))
        opts$v <- paste(names(d)[vi], collapse = ",")
    } else
        vi <- seq_along(d)
    if (!is.null(weights)) {
        # TDA's wt= is parsed by get_tp(), the time-points parser,
        # which demands a strictly increasing positive list -- reusing
        # that parser is why every wt= with equal or descending weights
        # dies as a syntax error, and why this option looked dead for so
        # long.  The distance itself is order-free, so the wrapper sorts:
        # variables are sent in ascending-weight order (v=), the weights
        # ascending with them.  Ties still cannot pass the parser; TDA's
        # own weighted-euclidean convention is d^2 = sum w_j * delta_j^2
        # (linear in w for city-block), so scale a tied variable yourself
        # (sqrt(w) * x for euclidean, w * x for city-block) and call this
        # unweighted.
        w <- as.numeric(weights)
        if (length(w) != length(vi))
            stop("need one weight per variable (", length(vi), ")",
                 call. = FALSE)
        if (any(w <= 0))
            stop("weights must be positive", call. = FALSE)
        if (anyDuplicated(w))
            stop("TDA's wt= parser cannot express tied weights ",
                 "(it requires a strictly increasing list); ",
                 "scale the tied columns instead -- see ?tda_pdatd",
                 call. = FALSE)
        o <- order(w)
        opts$v <- paste(names(d)[vi][o], collapse = ",")
        opts$wt <- paste(format(w[o], scientific = FALSE, trim = TRUE),
                         collapse = ",")
    }
    opts <- c(opts, .tda_extra(options))
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = "pdatd"), opts,
                                          list(rhs = "dist.out")))),
                   data = d, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run pdatd: ", err[1L], call. = FALSE)
    em0 <- if (.use_exports()) res$exports[["pdatd.table"]]
    m <- if (is.matrix(em0)) .export_frame(em0)
         else tda_file(res, "dist.out")
    if (is.null(m))
        stop("pdatd produced no output")
    stats::as.dist(as.matrix(m))
}
