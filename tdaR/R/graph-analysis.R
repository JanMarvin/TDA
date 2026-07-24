# Wrappers for TDA's graph-analysis commands.
# House pattern: edge data goes through gdd, then the command runs on
# the loaded graph; each wrapper is a thin, documented interface.

.tda_edge_setup <- function(edges, directed = TRUE, valued = FALSE) {
    edges <- as.data.frame(edges)
    stopifnot(ncol(edges) >= 2L)
    if (valued && ncol(edges) < 3L) edges$v <- 1
    list(dat = edges,
         vars = c("FROM[6.0]=c1", "TO[6.0]=c2",
                  if (valued) "V[12.4]=c3"),
         rhs = paste(c("FROM", "TO", if (valued) "V"), collapse = ","),
         # gdd's gt encodes the graph type: 1 undirected unvalued,
         # 2 undirected valued, 3 directed unvalued, 4 directed valued.
         # The original mapping sent directed graphs to gt=1 --
         # undirected, unvalued -- and dropped valuedness entirely.
         gt = 1L + (valued) + 2L * (directed))
}

.tda_graph_run <- function(edges, cmds, directed = TRUE, valued = FALSE,
                           gdd_opt = NULL, ...) {
    # gdd's opt encodes the graph type (1 undirected unvalued ..
    # 4 directed valued); derive it from the flags unless a caller
    # overrides.  Previously this was fixed at 1, so every wrapper
    # analysed an undirected, unvalued version of its input.
    if (is.null(gdd_opt))
        gdd_opt <- 1L                      # edge list input
    es <- .tda_edge_setup(edges, directed, valued)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(es$dat, file.path(dr, "e.dat"),
                       row.names = FALSE, col.names = FALSE)
    res <- tda_run(c(sprintf("nvar(dfile=e.dat, %s);",
                             paste(es$vars, collapse = ", ")),
                     sprintf("gdd(opt=%d, gt=%d) = %s;",
                             gdd_opt, es$gt, es$rhs),
                     cmds), dir = dr, ...)
    res$dir <- dr
    res
}

#' Undirected graph from node lists
#'
#' TDA's \code{gde}: builds an undirected graph from a data matrix
#' where \code{node} numbers each unit and \code{id} identifies the
#' objects they share; units sharing an id become adjacent.
#'
#' @param node,id integer vectors of equal length.
#' @param output "edges" for an edge list, "matrix" for the square
#'   adjacency matrix.
#' @param loops FALSE removes loops.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the values TDA wrote, as a matrix (one row per line of the
#'   output file), from the export channel.
#' @examples
#' # three persons (node) and the organisations they belong to (id):
#' # person 1 in 10 and 20, person 2 in 20 and 30, person 3 in 30.
#' # Sharing an organisation makes two persons adjacent: 1-2 and 2-3.
#' tda_gde(node = c(1, 1, 2, 2, 3), id = c(10, 20, 20, 30, 30))
#' @export
tda_gde <- function(node, id, output = c("edges", "matrix"),
                    loops = TRUE, ...) {
    output <- match.arg(output)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(data.frame(node, id), file.path(dr, "n.dat"),
                       row.names = FALSE, col.names = FALSE)
    res <- tda_run(c("nvar(dfile=n.dat, N[6.0]=c1, ID[6.0]=c2);",
              sprintf("gde(df=g.out, prn=%d, ni=%d) = N,ID;",
                      match(output, c("edges", "matrix")) - 1L,
                      as.integer(!loops))), dir = dr, ...)
    fv <- .file_values(res)
    if (!is.null(fv))
        return(fv)
    readLines(file.path(dr, "g.out"))
}

#' Union of two graphs
#'
#' TDA's \code{gdu}: graph 1 is loaded from \code{edges1} (with the
#' \code{perm} option of \code{gdd}), then the union with the edge
#' list \code{edges2} is written to the output file.
#'
#' @param edges1,edges2 data frames: from, to, value.
#' @param missing_value substitute written for missing values.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the output file's lines.
#' @examples
#' e1 <- data.frame(i = c(1, 2), j = c(2, 3), v = c(1, 1))
#' e2 <- data.frame(i = c(1, 3), j = c(3, 1), v = c(2, 2))
#' tda_gdu(e1, e2)   # value1/value2 side by side, -1 where absent
#' @export
tda_gdu <- function(edges1, edges2, missing_value = -1, ...) {
    es1 <- as.data.frame(edges1); es2 <- as.data.frame(edges2)
    if (ncol(es1) < 3L) es1$v <- 1
    if (ncol(es2) < 3L) es2$v <- 1
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(es1, file.path(dr, "e1.dat"),
                       row.names = FALSE, col.names = FALSE)
    utils::write.table(es2, file.path(dr, "e2.dat"),
                       row.names = FALSE, col.names = FALSE)
    res <- tda_run(c("nvar(dfile=e1.dat, I[6.0]=c1, J[6.0]=c2, V[12.4]=c3);",
                     "gdd(opt=1, perm=1) = I,J,V;",
                     "clear;",
                     "nvar(dfile=e2.dat, I[6.0]=c1, J[6.0]=c2, V[12.4]=c3);",
                     sprintf("gdu(df=u.out, sc=%g) = I,J,V;", missing_value)),
                   dir = dr, ...)
    u <- tda_file(res, "u.out")
    names(u) <- c("from", "to", "value1", "value2")
    u
}

#' Forward control flows in a directed graph
#'
#' TDA's \code{gfcf} on a directed graph loaded from an edge list.
#'
#' @param edges data frame: from, to.
#' @param threshold control threshold, default 0.5.
#' @param variant 1 or 2 (the two forms the command documents).
#' @param ... passed to \code{\link{tda_run}}.
#' @return one row per node: the node, how many nodes it controls,
#'   and the controlled nodes as a comma-joined string; the pair table
#'   TDA writes (node, controlled node, path length L, share S, M, R --
#'   the manual's Box 2 of 7.6.1.1) is the \code{"pairs"} attribute.
#' @examples
#' # 1 -> 2 -> 4 and 1 -> 3 -> 4: node 1 controls everything downstream,
#' # 2 and 3 control node 4 only, node 4 controls nothing
#' fc <- tda_gfcf(data.frame(i = c(1, 1, 2, 3), j = c(2, 3, 4, 4)))
#' fc
#' attr(fc, "pairs")   # per pair: path length and the flow shares
#' @export
tda_gfcf <- function(edges, threshold = 0.5, variant = 1, ...) {
    res <- .tda_graph_run(edges,
        sprintf("gfcf(gn=1, sc=%g, opt=%d, df=f.out) = g.out;", threshold, variant),
        directed = TRUE, valued = TRUE, ...)
    # g.out is the pair file of the manual's Box 2 (7.6.1.1): i, N(i),
    # N(j), L (path length), S (share), M, R -- one row per controlled
    # node; f.out summarises it per controlling node
    pairs <- .file_rows(res, "g.out") %||% .file_rows_text(res, "g.out")
    v <- .file_rows(res, "f.out") %||% .file_rows_text(res, "f.out")
    if (is.null(v) || !length(v))
        return(res$output)
    v <- lapply(v, as.integer)
    out <- data.frame(node = vapply(v, `[`, 0L, 1L),
                      n_controlled = vapply(v, `[`, 0L, 3L),
                      controlled = vapply(v, function(r)
                          paste(r[-(1:3)], collapse = ","), ""))
    if (!is.null(pairs) && length(pairs)) {
        p <- do.call(rbind, lapply(pairs, function(r) r[seq_len(7L)]))
        attr(out, "pairs") <- data.frame(i = as.integer(p[, 1L]),
                                         node = as.integer(p[, 2L]),
                                         controlled = as.integer(p[, 3L]),
                                         L = as.integer(p[, 4L]), S = p[, 5L],
                                         M = as.integer(p[, 6L]),
                                         R = as.integer(p[, 7L]))
    }
    out
}

#' Hasse diagram of multivariate patterns
#'
#' TDA's \code{ghd}: the rows of \code{x} are (assumed unique)
#' patterns; the command derives the dominance order's Hasse diagram.
#'
#' @param x matrix or data frame of patterns.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the edge list lines of the diagram.
#' @examples
#' # five patterns over two criteria; a row dominates another when it
#' # is at least as large in both.  The Hasse diagram keeps only the
#' # covering relations: (0,0) < (1,0), (0,1) < (1,1) < (2,1)
#' tda_ghd(rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1), c(2, 1)))
#' @export
tda_ghd <- function(x, ...) {
    x <- as.matrix(x)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(x, file.path(dr, "p.dat"),
                       row.names = FALSE, col.names = FALSE)
    vn <- sprintf("X%d[6.0]=c%d", seq_len(ncol(x)), seq_len(ncol(x)))
    tda_run(c(sprintf("nvar(dfile=p.dat, %s);", paste(vn, collapse = ", ")),
              sprintf("ghd(df=h.out) = %s;",
                      paste(sprintf("X%d", seq_len(ncol(x))), collapse = ","))),
            dir = dr, ...) -> res
    fv <- .file_values(res)
    if (!is.null(fv))
        return(fv)
    readLines(file.path(dr, "h.out"))
}

#' Quadratic assignment (GRASP)
#'
#' TDA's \code{gqap}: given a flow graph and a distance graph (both
#' undirected, integer-valued), search for the assignment minimizing
#' total cost, with the GRASP heuristic of CACM algorithm 754.
#'
#' @param flow,dist symmetric integer matrices of equal dimension.
#' @param max_iterations iteration cap.
#' @param ... passed to \code{\link{tda_run}}.
#' @return with direct exports enabled, a list with \code{cost},
#'   \code{assignment} and \code{output}; otherwise the printed
#'   output lines.
#' @examples
#' f <- rbind(c(0, 3, 1), c(3, 0, 2), c(1, 2, 0))
#' d <- rbind(c(0, 1, 4), c(1, 0, 2), c(4, 2, 0))
#' tda_gqap(f, d)
#' @export
tda_gqap <- function(flow, dist, max_iterations = 100, ...) {
    flow <- .tda_check_matrix(flow, "flow", square = TRUE,
                              symmetric = TRUE)
    dist <- .tda_check_matrix(dist, "dist", square = TRUE,
                              symmetric = TRUE, same_dim_as = flow,
                              other_arg = "flow")
    # canonical usage from the shipped co2 example: both matrices go
    # in as the two value columns of ONE lower-triangle gdd (opt=3),
    # forming a 2-valued multigraph
    lt <- function(m) {
        n <- nrow(m)
        unlist(lapply(2:n, function(i) m[i, 1:(i - 1)]))
    }
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(cbind(lt(flow), lt(dist)), file.path(dr, "q.dat"),
                       row.names = FALSE, col.names = FALSE)
    res <- tda_run(c("nvar(dfile=q.dat, F[8.0]=c1, D[8.0]=c2);",
                     "gdd(opt=3) = F,D;",
                     sprintf("gqap(mxit=%d, fmt=3.0, df=q.aux) = q.out;",
                             max_iterations)),
                   dir = dr, ...)
    ex <- res$exports
    # df=q.aux is asked for above and was then discarded. It is the
    # permuted multigraph -- the manual's second output file
    # (7.3.1.2, Box 2): each pair, its nodes under the permutation, and
    # the flow and distance between them.
    perm <- NULL
    aux <- file.path(res$dir, "q.aux")
    if (file.exists(aux)) {
        perm <- tryCatch(utils::read.table(aux, header = FALSE),
                         error = function(e) NULL)
        if (is.data.frame(perm) && ncol(perm) == 6L)
            perm <- .name_cols(perm, c("i", "j", "node_i", "node_j",
                                       "flow", "distance"))
    }
    if (.use_exports() && !is.null(ex[["gqap.cost"]])) {
        return(.tda_structured(
            list(cost = as.numeric(ex[["gqap.cost"]])[1L],
                 assignment = as.integer(ex[["gqap.assignment"]]),
                 permuted = perm,
                 output = res$output), "tda_assignment"))
    }
    # text mode: the cost from the console line, the assignment rows
    # from the file
    cost <- suppressWarnings(as.numeric(sub(".*: *([-0-9.eE+]+).*", "\\1",
        grep("Best cost value", res$output, value = TRUE)[1L])))
    rows <- .file_rows(res, "q.out") %||% .file_rows_text(res, "q.out")
    .tda_structured(list(cost = cost, assignment = unname(rows),
                         output = res$output), "tda_assignment")
}

#' Triangulation of points
#'
#' TDA's \code{triang} over an x,y point set.
#'
#' @param x,y coordinate vectors.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the output file's lines (the triangle list).
#' @examples
#' # a Delaunay triangulation of five points: the edges of the triangles
#' tda_triang(x = c(0, 2, 1, 3, 1), y = c(0, 0, 1, 1, 2))
#' @export
tda_triang <- function(x, y, ...) {
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(data.frame(x, y), file.path(dr, "t.dat"),
                       row.names = FALSE, col.names = FALSE)
    res <- tda_run(c("nvar(dfile=t.dat, X[12.4]=c1, Y[12.4]=c2);",
              "triang(prn=1, df=t.out) = X,Y;"), dir = dr, ...)
    fv <- .file_values(res)
    if (!is.null(fv))
        return(fv)
    readLines(file.path(dr, "t.out"))
}

#' Distances made metric
#'
#' TDA's \code{dmet1}: adjusts a distance matrix (given as an
#' undirected valued graph) toward metricity.
#'
#' @param d symmetric distance matrix.
#' @param tolerance,max_iterations deviation tolerance and cap.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the output file's lines.
#' @examples
#' # 5 breaks the triangle inequality (1 + 2 < 5), so the
#' # metricized matrix shortens it to 3
#' d <- rbind(c(0, 1, 5), c(1, 0, 2), c(5, 2, 0))
#' r <- tda_dmet1(d)
#' r$modified
#' @export
tda_dmet1 <- function(d, tolerance = 1e-4, max_iterations = 20, ...) {
    d <- .tda_check_matrix(d, "d", square = TRUE, symmetric = TRUE)
    ut <- do.call(rbind, lapply(seq_len(nrow(d) - 1), function(i)
        cbind(i, (i + 1):nrow(d), d[i, (i + 1):nrow(d)])))
    res <- .tda_graph_run(ut,
        sprintf("dmet1(gn=1, tolfd=%g, mxit=%d) = m.out;", tolerance, max_iterations),
        directed = FALSE, valued = TRUE, ...)
    # m.out holds titled sections -- the title lines carry no numbers,
    # so the tap's rows fall into blocks at them; the modified matrix
    # only appears when the input needed metricization
    bl <- .file_blocks(res, "m.out")
    if (!length(bl)) return(res$output)
    mat <- function(v) do.call(rbind, v)
    r <- list(original = mat(bl[[1L]]),
              modified = if (length(bl) > 1L) mat(bl[[2L]])
                         else mat(bl[[1L]]))
    .tda_structured(c(r, list(output = res$output)), "tda_dmet1")
}

#' Separable clusters of a valued graph
#'
#' TDA's \code{scla} on a distance matrix.
#'
#' @param d symmetric distance matrix.
#' @param ... passed to \code{\link{tda_run}}.
#' @return \code{clusters}, the rows TDA wrote (one numeric vector per
#'   line of the separable-cluster listing), and \code{output}.
#' @examples
#' d <- as.matrix(dist(c(0, 0.1, 5, 5.1)))
#' r <- tda_scla(d)
#' r$clusters
#' @export
tda_scla <- function(d, ...) {
    d <- .tda_check_matrix(d, "d", square = TRUE, symmetric = TRUE)
    ut <- do.call(rbind, lapply(seq_len(nrow(d) - 1), function(i)
        cbind(i, (i + 1):nrow(d), d[i, (i + 1):nrow(d)])))
    res <- .tda_graph_run(ut, "scla(gn=1) = s.out;",
                          directed = FALSE, valued = TRUE, ...)
    rows <- .file_rows(res, "s.out") %||% .file_rows_text(res, "s.out")
    .tda_structured(list(clusters = unname(rows), output = res$output),
                    "tda_scla")
}

#' Column permutations toward graph agreement
#'
#' TDA's \code{gap} (CACM algorithm 548) on two integer-valued
#' graphs, given as matrices the way \code{\link{tda_gqap}} takes
#' them.
#'
#' @param g1,g2 symmetric integer matrices of equal dimension.
#' @param ... passed to \code{\link{tda_run}}.
#' @return \code{permuted}, the permuted matrix, and \code{output}.
#' @examples
#' f <- rbind(c(0, 3, 1), c(3, 0, 2), c(1, 2, 0))
#' d <- rbind(c(0, 1, 4), c(1, 0, 2), c(4, 2, 0))
#' r <- tda_gap_permute(f, d)
#' r$permuted
#' @export
tda_gap_permute <- function(g1, g2, ...) {
    g1 <- .tda_check_matrix(g1, "g1", square = TRUE, symmetric = TRUE)
    g2 <- .tda_check_matrix(g2, "g2", square = TRUE, symmetric = TRUE,
                            same_dim_as = g1, other_arg = "g1")
    lt <- function(m) {
        n <- nrow(m)
        unlist(lapply(2:n, function(i) m[i, 1:(i - 1)]))
    }
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(cbind(lt(g1), lt(g2)), file.path(dr, "g.dat"),
                       row.names = FALSE, col.names = FALSE)
    res <- tda_run(c("nvar(dfile=g.dat, F[8.0]=c1, D[8.0]=c2);",
                     "gdd(opt=3) = F,D;",
                     "gap(df=g.aux) = g.out;"), dir = dr, ...)
    rows <- .file_rows(res, "g.out") %||% .file_rows_text(res, "g.out")
    .tda_structured(list(permuted = if (length(rows)) do.call(rbind, unname(rows)),
                         output = res$output), "tda_gap")
}

#' The assignment problem
#'
#' TDA's \code{gap} on a cost matrix: the permutation that assigns each
#' row to a column at minimal total cost (Carpaneto and Toth's
#' algorithm), the manual's 7.3.1.1. The matrix goes to TDA as a full
#' matrix graph (\code{gdd(opt=7)}), so it need not be symmetric.
#'
#' @param cost a square numeric cost matrix.
#' @param ... passed to \code{\link{tda_run}}.
#' @return \code{assignment}, the permutation as a two-column matrix
#'   (i, p(i)); \code{permuted}, the cost matrix with its rows permuted
#'   accordingly; and \code{output}.
#' @examples
#' cost <- rbind(c(60, 0, 0, 76, 0, 0), c(0, 40, 18, 0, 60, 24),
#'               c(60, 16, 2, 4, 0, 40), c(0, 27, 18, 3, 55, 75),
#'               c(0, 40, 62, 16, 11, 53), c(28, 4, 10, 84, 0, 16))
#' tda_gap(cost)$assignment
#' @export
tda_gap <- function(cost, ...) {
    cost <- .tda_check_matrix(cost, "cost", square = TRUE)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(as.vector(t(cost)), file.path(dr, "g.dat"),
                       row.names = FALSE, col.names = FALSE)
    res <- tda_run(c("nvar(dfile=g.dat, D[8.0]=c1);",
                     "gdd(opt=7) = D;",
                     "gap(df=g.aux, fmt=3.0) = g.out;"), dir = dr, ...)
    rows <- .file_rows(res, "g.out") %||% .file_rows_text(res, "g.out")
    perm <- .file_rows(res, "g.aux") %||% .file_rows_text(res, "g.aux")
    .tda_structured(list(assignment = if (length(rows)) do.call(rbind, unname(rows)),
                         permuted = if (length(perm)) do.call(rbind, unname(perm)),
                         output = res$output), "tda_gap")
}

#' Projection of proximities
#'
#' TDA's \code{gpro} on a proximity matrix given as an undirected
#' valued graph.
#'
#' @param d symmetric proximity matrix.
#' @param max_iterations iteration cap.
#' @param ... passed to \code{\link{tda_run}}.
#' @return \code{projection}, the projected coordinates, and
#'   \code{output}.
#' @examples
#' r <- tda_gpro(as.matrix(dist(cbind(c(0, 0, 3, 3), c(0, 2, 0, 2)))))
#' r$projection
#' @export
tda_gpro <- function(d, max_iterations = 50, ...) {
    d <- as.matrix(d)
    ut <- do.call(rbind, lapply(seq_len(nrow(d) - 1), function(i)
        cbind(i, (i + 1):nrow(d), d[i, (i + 1):nrow(d)])))
    res <- .tda_graph_run(ut,
        sprintf("gpro(gn=1, mxit=%d) = p.out;", max_iterations),
        directed = FALSE, valued = TRUE, ...)
    rows <- .file_rows(res, "p.out") %||% .file_rows_text(res, "p.out")
    .tda_structured(list(projection = if (length(rows)) do.call(rbind, unname(rows)),
                         output = res$output), "tda_gpro")
}

#' Temporal-network summaries
#'
#' TDA's \code{tnet} on a directed edge list.
#'
#' @param edges data frame: from, to, then at least three value
#'   columns in the order VALIDITY (>= 0 keeps the edge), then the
#'   edge's valid-from and valid-to dates with from <= to --
#'   established from the source's check; the command requires a
#'   multigraph with three or more subgraphs.
#' @param report "basic", "subgraphs", "indegrees", "outdegrees", or
#'   "layers".
#' @param ... passed to \code{\link{tda_run}}.
#' @return the printed output.
#' @examples
#' r <- tda_tnet(data.frame(i = c(1, 2, 3), j = c(2, 3, 4),
#'                          valid = c(1, 1, 1),
#'                          from = c(10, 12, 14), to = c(20, 22, 24)))
#' cat(head(tda_payload(r), 12), sep = "\n")
#' @export
tda_tnet <- function(edges, report = c("basic", "subgraphs",
                     "indegrees", "outdegrees", "layers"), ...) {
    report <- match.arg(report)
    edges <- as.data.frame(edges)
    stopifnot("tnet needs from, to and >= 3 layer columns" =
              ncol(edges) >= 5L)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(edges, file.path(dr, "e.dat"),
                       row.names = FALSE, col.names = FALSE)
    nl <- ncol(edges) - 2L
    vl <- sprintf("L%d[8.0]=c%d", seq_len(nl), 2L + seq_len(nl))
    res <- tda_run(c(sprintf("nvar(dfile=e.dat, FROM[6.0]=c1, TO[6.0]=c2, %s);",
                             paste(vl, collapse = ", ")),
                     sprintf("gdd(opt=1, gt=1) = FROM,TO,%s;",
                             paste(sprintf("L%d", seq_len(nl)),
                                   collapse = ",")),
                     sprintf("tnet(opt=%d);",
                             match(report, c("basic", "subgraphs",
                                   "indegrees", "outdegrees", "layers")))),
                   dir = dr, ...)
    res$output
}

#' Adjacency matrix of a graph
#'
#' TDA's \code{mdefg}: builds the adjacency matrix of the current
#' graph, with \code{missing} substituted for absent edges (TDA's
#' \code{sc=}, default -1 there; 0 here, which gives the usual
#' weighted adjacency matrix).
#'
#' @param edges edge list (data frame with from, to and optionally a
#'   value column).
#' @param missing value written where no edge exists.
#' @param directed,valued graph type flags, as in \code{\link{tda_g}}.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the numeric adjacency matrix.
#' @examples
#' e <- data.frame(from = c(1, 1, 2), to = c(2, 3, 3),
#'                 value = c(5, 2, 7))
#' tda_mdefg(e)
#' @export
tda_mdefg <- function(edges, missing = 0, directed = FALSE,
                      valued = TRUE, ...) {
    res <- .tda_graph_run(edges,
        c(sprintf("mdefg(gn=1, sc=%g) = M;", missing),
          "mfmt=15.8;", "mpr(M) = m.out;"),
        directed = directed, valued = valued, ...)
    m <- .mpr_from_exports(res)
    if (is.matrix(m))
        return(unname(m))
    mf <- file.path(res$dir, "m.out")
    if (!file.exists(mf))
        stop("TDA could not build the adjacency matrix: ",
             grep("rror", res$output, value = TRUE)[1L], call. = FALSE)
    unname(as.matrix(tda_file(res, mf)))
}

