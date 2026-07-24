# Wrapping TDA's graph commands.
#
# There are forty of them and they share a shape: a graph is defined once with
# gdd, then commands are run against it and each writes a table.  So the
# wrapper is a graph object, built from an edge list or an adjacency matrix,
# with one function per command on top of a generic that reaches all of them.
#
# gdd takes the graph in one of several forms, chosen with opt=, and gt=2
# marks it undirected.  Which form a command needs varies: some want an edge
# list, some a dissimilarity matrix, some insist on an undirected graph, and
# they say so when given the wrong one.

#' A graph for TDA's graph commands
#'
#' @param edges an edge list -- from, to, and one value, or two values when
#'   the graph carries a different one in each direction, as the shipped
#'   \code{gd1.dat} does -- or a square adjacency matrix, or a \code{dist}.
#' @param directed whether the graph is directed. \code{FALSE} emits
#'   \code{gt=2}. Several commands -- \code{\link{tda_g_components}} among
#'   them, and \code{\link{plot.tda_graph}}'s \code{layout = "tree"} --
#'   only work on an undirected graph, and TDA itself refuses a directed
#'   one for them; \code{directed = TRUE} (the default here) is right for
#'   ordinary, arrow-drawn graphs, but a tree needs \code{FALSE}, since a
#'   tree's layout comes from which nodes are connected, not which way
#'   an edge happens to point.
#' @param form how TDA is to read it: \code{"edges"} for an edge list, which
#'   is \code{gdd} option 1, or \code{"dissimilarity"} for the upper triangle
#'   of a dissimilarity matrix, option 4.
#' @param options a named list of further \code{gdd} options.
#' @return A \code{tda_graph}.
#' @family graph analysis
#' @examples
#' e <- data.frame(from = c(1, 1, 2, 3), to = c(2, 3, 4, 4), value = 1)
#' # gcon works on an undirected graph only, so `directed = FALSE`
#' tda_g_components(tda_graph(e, directed = FALSE))
#' @export
tda_graph <- function(edges, directed = TRUE, form = c("edges",
                      "dissimilarity"), options = list()) {
    form <- match.arg(form)
    if (inherits(edges, "dist")) {
        d <- data.frame(D = as.vector(edges))
        form <- "dissimilarity"
        rhs <- "D"
    } else if (is.matrix(edges) && nrow(edges) == ncol(edges) &&
               form == "dissimilarity") {
        d <- data.frame(D = edges[lower.tri(edges)])
        rhs <- "D"
    } else {
        e <- as.data.frame(edges)
        if (ncol(e) < 2L)
            stop("`edges` needs at least a from and a to column")
        if (ncol(e) == 2L)
            e[[3L]] <- 1
        # An edge list may carry two values, one per direction, which is what
        # the shipped gd1.dat does: "gdd = I,J,V1,V2".  Four columns are kept
        # rather than truncated to three.
        e <- e[, seq_len(min(4L, ncol(e))), drop = FALSE]
        names(e) <- c("I", "J", "V1", "V2")[seq_len(ncol(e))]
        d <- e
        rhs <- paste(names(e), collapse = ",")
    }
    opts <- .tda_extra(options)
    if (form == "dissimilarity")
        opts$opt <- 4
    else
        opts$opt <- 1
    if (!directed)
        opts$gt <- 2
    structure(list(data = d, rhs = rhs, opts = opts, directed = directed,
                   form = form), class = "tda_graph")
}

#' Run a graph command
#'
#' \code{tda_g} reaches any of TDA's forty graph commands; the named functions
#' are the ones with an obvious meaning.
#'
#' @param g a \code{\link{tda_graph}}.
#' @param cmd the command name, e.g. \code{"gsp"}, \code{"gmst"}.
#' @param options a named list of options for the command.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object carrying a \code{table} and the run.
#' @family graph analysis
#' @examples
#' e <- data.frame(from = c(1, 1, 2, 3), to = c(2, 3, 4, 4), value = 1)
#' g <- tda_graph(e, directed = FALSE)
#' tda_g(g, "gmst")$table  # minimum spanning tree
#' @export
tda_g <- function(g, cmd, options = list(), dir = tempfile("tda"), ...) {
    stopifnot(inherits(g, "tda_graph"))
    # tda_block(name, ..., rhs=) took the command name positionally here,
    # via an unnamed list(cmd)/list("gdd") -- fine until an option is
    # itself named something that partially matches "name" (gcd's
    # n=, say: do.call()'s argument matching resolves n= to name=
    # first, leaving the actual command name to fall into ... instead,
    # and TDA gets sent "4(gcd,)" rather than "gcd(n=4,)"). The same
    # class of bug found and fixed elsewhere this session (lt=/lty,
    # r=/rot), caught here by testing tda_g_random's n= directly
    # rather than assumed safe because the pattern looked routine.
    # Naming it explicitly closes this off for every command, not just
    # the one that happened to expose it.
    res <- tda_run(c(tda_nvar(g$data),
                     do.call(tda_block, c(list(name = "gdd"), g$opts,
                                          list(rhs = g$rhs))),
                     do.call(tda_block, c(list(name = cmd),
                                          .tda_extra(options),
                                          list(rhs = "out.txt")))),
                   data = g$data, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run ", cmd, ": ", err[1L], call. = FALSE)
    out <- .read_g_output(res, "out.txt", text = identical(cmd, "gdot"))
    # A command with its own producer overlays the parsed frame: the
    # graph files are written at four decimals, so a ratio column (gnc's
    # density, 2/3) came back as 0.6667.  Only the commands that HAVE a
    # producer are touched; the rest keep the parsed table unchanged.
    if (!is.null(out$table))
        out$table <- .overlay_num(out$table,
                                  res$exports[[paste0(cmd, ".table")]])
    if (is.data.frame(out$table))
        out$table <- .g_name_columns(cmd, options, out$table)
    structure(list(call = sys.call(-1L), run = res, command = cmd,
                   n = nrow(g$data),
                   # the same graph-summary row mds/cluster store; the
                   # helper lives in clustering.R
                   graph = .graph_summary(res$output),
                   table = out$table, text = out$text),
              class = c(paste0("tda_", cmd), "tda_table"))
}

# A handful of graph commands write output tda_file()'s fixed-width reader
# cannot parse: some (gdln, links per node) a ragged table -- a node, its
# link count, then that many linked node ids, the same shape as hcld's
# cluster tree -- and others (gdot) plain text (Graphviz source) that is not
# tabular at all. Both were previously silently dropped ($table stayed NULL
# with no trace of the "N records written" TDA had reported). Falls through
# fixed-width -> ragged numeric -> raw text, so the file's contents always
# end up in $table or $text rather than nowhere.
# Column names for the rectangular graph outputs, from the headers the
# manual prints over each command's output file (chapter 7): i and j
# are TDA's node indices, N(i) and N(j) the node numbers.  Applied only
# when the width is the one the command documents; anything else keeps
# TDA's bare columns.
.g_name_columns <- function(cmd, opts, tab) {
    nc <- ncol(tab)
    opt <- opts$opt %||% 1
    nm <- switch(cmd,
        gsort = if (nc == 3L) c("index", "node", "level"),
        gcon = , gdcon = if (nc == 4L && opt == 1)
            c("component", "size", "index", "node"),
        gcut = if (nc == 2L) c("index", "node"),
        gep = if (opt == 2 && nc == 8L)
            c("i", "j", "node_i", "node_j", "n_paths", "n_simple_paths",
              "min_length", "min_value"),
        # gsp option 1: after i and N(i), (reachable node, distance) pairs
        gsp = if (opt == 1 && nc > 2L && (nc - 2L) %% 2L == 0L)
            c("index", "node", paste0(c("to", "dist"), rep(seq_len((nc - 2L) / 2L), each = 2L))),
        gtcl = if (opt == 2 && nc > 2L)
            c("index", "node", paste0("n", seq_len(nc - 2L))),
        gst = , gmst = if (nc == 4L) c("tree", "node_i", "node_j", "value"),
        gflow = if (opt == 1 && nc == 6L)
            c("index", "i", "j", "node_i", "node_j", "flow"),
        gfc = if (nc >= 6L)
            c("i", "j", "node_i", "node_j", paste0("f", seq_len(nc - 4L))),
        becl = if (opt == 1 && nc == 4L) c("index", "node", "position", "node_at"),
        NULL)
    if (!is.null(nm) && length(nm) == nc)
        names(tab) <- nm
    tab
}

.read_g_output <- function(res, file, text = FALSE) {
    # from the tap: a rectangular table, or ragged rows (a node and its
    # links); only plain text (gdot's Graphviz source) is read from disk
    if (text) {
        lines <- .file_lines(res, file)
        if (is.null(lines)) {
            p <- file.path(res$dir, file)
            if (!file.exists(p) || file.size(p) == 0L)
                return(list(table = NULL, text = NULL))
            lines <- readLines(p, warn = FALSE)
        }
        return(list(table = NULL, text = lines[nzchar(trimws(lines))]))
    }
    tab <- .file_from_exports(res, file)
    if (!is.null(tab))
        return(list(table = tab, text = NULL))
    rows <- .file_rows(res, file)
    if (!is.null(rows)) {
        class(rows) <- "tda_ragged"
        return(list(table = rows, text = NULL))
    }
    p <- file.path(res$dir, file)
    if (!file.exists(p) || file.size(p) == 0L)
        return(list(table = NULL, text = NULL))
    tab <- tryCatch(tda_file(res, file), error = function(e) NULL)
    if (!is.null(tab))
        return(list(table = tab, text = NULL))
    lines <- readLines(p, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]
    if (!length(lines))
        return(list(table = NULL, text = NULL))
    ragged <- lapply(lines, function(l)
        suppressWarnings(as.numeric(strsplit(trimws(l), "\\s+")[[1L]])))
    if (!any(vapply(ragged, anyNA, logical(1)))) {
        class(ragged) <- "tda_ragged"
        return(list(table = ragged, text = NULL))
    }
    list(table = NULL, text = lines)
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::print
print.tda_ragged <- function(x, n = 10L, ...) {
    len <- lengths(x)
    lenspec <- if (length(unique(len)) == 1L) paste0(": ", len[1L])
               else paste0("s ", min(len), "-", max(len))
    cat(sprintf("%d record%s, length%s\n", length(x),
                if (length(x) == 1L) "" else "s", lenspec))
    show <- x[seq_len(min(n, length(x)))]
    for (i in seq_along(show))
        cat(sprintf("[%d] ", i), format(show[[i]]), "\n")
    if (length(x) > n)
        cat("... and", length(x) - n, "more\n")
    invisible(x)
}

#' Graph analyses
#'
#' One function per graph command, each a thin wrapper: the command name, the
#' options it takes and the table TDA writes. Everything TDA offers that has
#' no named wrapper is still reachable with \code{\link{tda_g}}.
#'
#' The descriptions below are Rohwer's, from TDA's manual; the exact
#' wording and the full option list for any command are available with
#' \code{tda_help("gmst")} and so on.
#'
#' @section Connectivity and structure:
#' \describe{
#'   \item{\code{tda_g_components}}{connected components (\code{gcon}).}
#'   \item{\code{tda_g_reachable}}{for each node of a directed graph, the
#'     set of nodes reachable from it by a directed path (\code{gdcon}).}
#'   \item{\code{tda_g_cutpoints}}{cut nodes of an undirected graph
#'     (\code{gcut}).}
#'   \item{\code{tda_g_transitive}}{transitive closure (\code{gtcl}).}
#'   \item{\code{tda_g_symmetric}}{whether every edge from i to j is
#'     matched by one from j to i (\code{gsym}).}
#'   \item{\code{tda_g_toposort}}{topological sort of a directed graph
#'     (\code{gsort}).}
#' }
#'
#' @section Subgraphs:
#' \describe{
#'   \item{\code{tda_g_cliques}}{all cliques -- maximally connected
#'     subgraphs -- of an undirected graph (\code{gcliq}).}
#'   \item{\code{tda_g_independent}}{a maximal independent set, found by a
#'     greedy randomized search procedure (GRASP), so the answer is a good
#'     one rather than a proven optimum (\code{giset}).}
#'   \item{\code{tda_g_compact}}{compact sets in an undirected valued graph
#'     (\code{gcset}).}
#'   \item{\code{tda_g_neighbourhoods}}{compares the neighbourhoods of the
#'     nodes and groups those that agree. Needs an unvalued graph; the
#'     wrapper arranges that (\code{gcni}).}
#' }
#'
#' @section Paths, trees and cycles:
#' \describe{
#'   \item{\code{tda_g_shortest}}{a shortest path for each pair of nodes.
#'     Directed or undirected, valued or unvalued; in an unvalued graph
#'     every edge counts as 1 (\code{gsp}).}
#'   \item{\code{tda_g_paths}}{all paths between any two nodes, with their
#'     minimal or maximal length (\code{gep}).}
#'   \item{\code{tda_g_mst}}{a minimum, or maximum, spanning tree for each
#'     connected component of an undirected valued graph (\code{gmst}).}
#'   \item{\code{tda_g_spantrees}}{enumerates every spanning tree of an
#'     undirected graph (\code{gnst}).}
#'   \item{\code{tda_g_cycles}}{a fundamental set of cycles, and optionally
#'     all cycles, of an undirected graph (\code{gcyc}).}
#'   \item{\code{tda_g_dcycles}}{cycles in a directed graph
#'     (\code{gdcyc}).}
#'   \item{\code{tda_g_links}}{forward and backward links in a directed
#'     graph (\code{gdln}).}
#' }
#'
#' @section Flow, control and ownership:
#' \describe{
#'   \item{\code{tda_g_flow}}{a maximal flow between each pair of nodes
#'     joined by a directed path (\code{gflow}).}
#'   \item{\code{tda_g_flowcontrol}}{the quantities from which measures of
#'     flow control are built, for a directed valued graph (\code{gfc}).}
#'   \item{\code{tda_g_ownership}}{integrated ownership in a directed
#'     valued graph (\code{gio}).}
#'   \item{\code{tda_g_backward}}{direct and indirect backward control
#'     (\code{gbcf}).}
#'   \item{\code{tda_g_centred}}{node-centred networks (\code{gnc}).}
#' }
#'
#' @section Description and output:
#' \describe{
#'   \item{\code{tda_g_nodes}}{information about the nodes of the current
#'     graph (\code{gni}).}
#'   \item{\code{tda_g_eigen}}{eigenvalues and eigenvectors of the
#'     adjacency matrix of an undirected graph (\code{gev}).}
#'   \item{\code{tda_g_aggregate}}{aggregates the nodes of a valued graph
#'     (\code{gda}).}
#'   \item{\code{tda_g_dblocks}}{blocks of a directed graph
#'     (\code{gdcset}).}
#'   \item{\code{tda_g_random}}{generates an unvalued graph with a given
#'     number of nodes and edges (\code{gcd}).}
#'   \item{\code{tda_g_dot}, \code{tda_g_edges}}{write the graph out
#'     as a Graphviz dot file, or an edge list/adjacency matrix,
#'     respectively (\code{gdot}, \code{gdp}).}
#'   \item{\code{tda_g_degrees}}{node degrees from \code{gni}, TDA's
#'     node-information command: in-degree, out-degree, loops, and a
#'     degree column (for an undirected graph each incident edge is
#'     reported on both sides, so the degree is the in-degree, not
#'     the sum). Takes no further arguments.}
#' }
#'
#' @section What each command requires:
#' TDA refuses a command whose graph is of the wrong kind rather than
#' answering wrongly, so the error message is worth reading. Most of the
#' \code{gd*} commands want a directed graph, \code{gcliq}, \code{gcut},
#' \code{gcyc}, \code{gmst}, \code{gnst}, \code{gcset} and \code{gnc}
#' want an undirected one, and \code{gcni} wants an unvalued one.
#'
#' @family graph analysis
#' @examples
#' # a small directed graph: 1 and 2 both feed 3, which feeds 4, which
#' # feeds 5
#' e <- data.frame(from = c(1, 1, 2, 2, 3, 4),
#'                 to   = c(2, 3, 3, 4, 4, 5),
#'                 value = c(1, 2, 1, 3, 2, 1))
#'
#' # undirected questions want an undirected graph
#' g <- tda_graph(e, directed = FALSE)
#' tda_g_components(g)$table
#' tda_g_cliques(g)$table
#' tda_g_mst(g)$table
#' tda_g_cycles(g)$table
#' tda_g_cutpoints(g)$table
#' tda_g_compact(g)$table
#' tda_g_spantrees(g)$table
#' tda_g_spanning(g)$table
#' tda_g_independent(g)$table
#' tda_g_centred(g)$table
#'
#' # directed ones want a directed graph
#' d <- tda_graph(e, directed = TRUE)
#' tda_g_toposort(d)$table
#' tda_g_reachable(d)$table
#' tda_g_transitive(d)$table
#' tda_g_nodes(d)$table
#' tda_g_links(d)$table
#' tda_g_paths(d)$table
#' tda_g_flow(d)$table
#' tda_g_flowcontrol(d)$table
#' tda_g_backward(d)$table
#' tda_g_degrees(d)[, c("node", "in_degree", "out_degree")]
#' tda_g_edges(d)$table         # the raw edge dump (gdp)
#' tda_g_random(d)$table        # a fresh random graph, ignores e's edges
#' tda_g_dot(d)$text            # Graphviz source, not a table -- see tda_g()
#' tda_g_shortest(g)$table      # all-pairs shortest paths
#'
#' # a denser graph, with an actual triangle and edges valued in [0, 1] --
#' # what gcliques (a real clique, not just a connected pair) and gio (an
#' # ownership fraction) both need to find something on this small example
#' e2 <- data.frame(from = c(1, 2, 3, 1, 4), to = c(2, 3, 1, 4, 5),
#'                  value = c(0.3, 0.5, 0.8, 0.6, 0.9))
#' g2 <- tda_graph(e2, directed = FALSE)
#' d2 <- tda_graph(e2, directed = TRUE)
#' tda_g_gcliques(g2)$table
#' tda_g_ownership(d2)$table
#'
#' # aggregate collapses named nodes into new ones -- rcn=, TDA's
#' # recode syntax, is required; without it there is nothing to do
#' tda_g_aggregate(d, rcn = "1[1,2],2[3,4,5]")$table
#'
#' # neighbourhoods wants an unvalued graph (gcni's restriction)
#' e3 <- data.frame(from = c(1, 1, 2, 2, 3, 4), to = c(2, 3, 3, 4, 4, 5))
#' tda_g_neighbourhoods(tda_graph(e3, directed = FALSE))$table
#'
#' # this graph has no cycles (it is a DAG) and every edge already agrees
#' # with itself in one direction, so these two are empty rather than
#' # broken -- not every command has something to report
#' tda_g_dcycles(d)$table
#' tda_g_symmetric(d)$run$output[grepl("Nodes|symmetr", d$run$output)]
#' tda_g_dblocks(d)$table  # 0 compact blocks in this graph, also real
#' tda_g_eigen(g)          # small graphs like these often do not converge
#'
#' # anything without a named wrapper goes through tda_g()
#' tda_g(d, "gdp")$table
#'
#' @param g a \code{\link{tda_graph}}.
#' @param edges for \code{tda_g_components}, report each component as an
#'   edge list instead of a node list (\code{opt=2}). For
#'   \code{tda_g_spantrees}, an edge list instead of one record per tree
#'   (\code{gnst}'s \code{opt=2}).
#' @param format for \code{tda_g_shortest}: \code{"reachable"} (default,
#'   only reachable pairs), \code{"matrix"} (a square matrix with row/
#'   column labels), or \code{"matrix_no_labels"} -- \code{gsp}'s
#'   \code{opt=}. For \code{tda_g_cycles}: \code{"fundamental_nodes"}
#'   (default), \code{"fundamental_edges"}, \code{"all_v1"}, or
#'   \code{"all_v2"} -- \code{gcyc}'s \code{opt=}. For
#'   \code{tda_g_edges}: \code{"edges"} (default), or the adjacency
#'   matrix as \code{"lower_triangle"}, \code{"lower_triangle_diag"},
#'   \code{"full"}, or \code{"full_square"} -- \code{gdp}'s
#'   \code{opt=}. For \code{tda_g_dcycles}: \code{"per_cycle"} (default),
#'   \code{"per_node"}, or \code{"progress"} (report progress to the
#'   console instead of writing anything) -- \code{gdcyc}'s
#'   \code{opt=}. For \code{tda_g_reachable}: \code{"count"} (default,
#'   number of reachable nodes), \code{"list"} (plus the node numbers),
#'   \code{"strong_nodes"}, or \code{"strong_edges"} (strongly connected
#'   components) -- \code{gdcon}'s \code{opt=}. Five different
#'   commands' own \code{opt=}, each with its own meaning.
#' @param algorithm for \code{tda_g_mst}: \code{"kruskal"} (default) or
#'   \code{"prim"}. For \code{tda_g_cliques}: \code{"bron_kerbosch"}
#'   (default) or \code{"harary_ross"}. Different commands' own
#'   \code{alg=}.
#' @param max for \code{tda_g_mst}, find a maximum spanning tree instead
#'   of a minimum one (\code{max=1}).
#' @param sort for \code{tda_g_mst} and \code{tda_g_cliques}, sort the
#'   output by node number.
#' @param matrix for \code{tda_g_flow}, report the result as a square
#'   flow matrix instead of one record per flow (\code{opt=2}).
#' @param min_size for \code{tda_g_cliques}, the minimum clique size
#'   (\code{gcliq}'s \code{min=}, default 3). For
#'   \code{tda_g_independent}, the minimum size of the independent set
#'   to search for (\code{giset}'s \code{min=}).
#' @param per_node for \code{tda_g_cutpoints}, one record per node
#'   instead of one per component (\code{gcut}'s \code{opt=2}). For
#'   \code{tda_g_compact}, the same idea for \code{gcset}'s
#'   \code{opt=2}.
#' @param relabel_graph for \code{tda_g_toposort}, return the graph with
#'   its new labels instead of the list of old/new labels (\code{gsort}'s
#'   own \code{opt=2}).
#' @param backward for \code{tda_g_links}, backward links instead of
#'   forward ones (\code{gdln}'s \code{opt=2}).
#' @param stop_at_first for \code{tda_g_dcycles}, stop after the first
#'   cycle is found (\code{gdcyc}'s \code{dopt=1}).
#' @param use_complement for \code{tda_g_independent}, search the
#'   complementary graph instead (\code{giset}'s \code{cg=1}).
#' @param n for \code{tda_g_eigen}, the number of eigenvalues/vectors to
#'   compute (\code{gev}'s \code{n=}, default 1).
#' @param tol for \code{tda_g_eigen}, the required accuracy
#'   (\code{gev}'s \code{eps=}).
#' @param control for \code{tda_g_eigen}, convergence settings from
#'   \code{\link{tda_control}} (only its \code{maxit} field applies,
#'   as \code{gev}'s \code{mxit=}).
#' @param ... options for the command.
#' @return An object carrying \code{table} -- the records the command
#'   writes, in that command's layout (edges for the tree and path
#'   commands, node rows for the centrality and component ones; see the
#'   command's \code{tda_help()} entry) -- or, for the few commands
#'   whose output is ragged or non-tabular, \code{table} as a ragged
#'   record list or \code{text} as plain lines.
#'
#'   \code{tda_g_backward} also carries \code{summary}, one row per node,
#'   from the second file \code{gbcf} writes; TDA does not produce it
#'   unless a file is named for it, so the wrapper asks by default.
#' @name tda_g_analyses
NULL

#' @rdname tda_g_analyses
#' @export
tda_g_components <- function(g, edges = FALSE, ...) {
    opts <- list(...)
    if (isTRUE(edges)) opts$opt <- 2
    r <- tda_g(g, "gcon", opts)
    # manual 7.2.3.1 Box 1 names both layouts: the node list, and the
    # edge list with each edge's value
    t <- r$table
    if (is.data.frame(t)) {
        if (ncol(t) == 4L)
            r$table <- .name_cols(t, c("component", "n", "i", "node_i"))
        else if (ncol(t) == 7L)
            r$table <- .name_cols(t, c("component", "n", "i", "j",
                                       "node_i", "node_j", "value"))
    }
    r
}

#' @rdname tda_g_analyses
#' @export
tda_g_shortest <- function(g, format = c("reachable", "matrix",
                                        "matrix_no_labels"), ...) {
    opts <- list(...)
    if (!missing(format))
        opts$opt <- match(match.arg(format), c("reachable", "matrix",
                                                "matrix_no_labels"))
    tda_g(g, "gsp", opts)
}

#' @rdname tda_g_analyses
#' @export
tda_g_mst <- function(g, algorithm = c("kruskal", "prim"), max = FALSE,
                      sort = FALSE, ...) {
    opts <- list(...)
    if (!missing(algorithm))
        opts$alg <- match(match.arg(algorithm), c("kruskal", "prim"))
    if (isTRUE(max)) opts$max <- 1
    if (isTRUE(sort)) opts$sort <- ""
    tda_g(g, "gmst", opts)
}

#' @rdname tda_g_analyses
#' @export
tda_g_cliques <- function(g, algorithm = c("bron_kerbosch", "harary_ross"),
                          min_size = NULL, sort = FALSE, ...) {
    opts <- list(...)
    if (!missing(algorithm))
        opts$alg <- match(match.arg(algorithm), c("bron_kerbosch",
                                                   "harary_ross"))
    if (!is.null(min_size)) opts$min <- min_size
    if (isTRUE(sort)) opts$sort <- ""
    tda_g(g, "gcliq", opts)
}

#' @rdname tda_g_analyses
#' @export
tda_g_cycles <- function(g, format = c("fundamental_nodes",
                                      "fundamental_edges", "all_v1",
                                      "all_v2"), ...) {
    opts <- list(...)
    fmt <- if (missing(format)) 1L else
        match(match.arg(format), c("fundamental_nodes", "fundamental_edges",
                                   "all_v1", "all_v2"))
    if (!missing(format))
        opts$opt <- fmt
    r <- tda_g(g, "gcyc", opts)
    # Every option writes the cycle number and an index first; what
    # follows differs (manual 7.2.6.1, Boxes 1 to 4).  Left unnamed the
    # caller gets V1, V2, ... and has to count columns to find the edge.
    #   1  the nodes of each fundamental cycle -- ragged, so left as is
    #   2  the edge, then one indicator per fundamental cycle
    #   3  the edge, then one indicator per cycle found
    #   4  the edge
    t <- r$table
    # option 4 writes one fixed-width row per edge, but it arrives
    # through the ragged tap, so the caller got a list of vectors where
    # a table was meant.  Rectangular rows become a data frame; a run
    # that really is ragged (option 1, a cycle's nodes) is left alone.
    if (inherits(t, "tda_ragged") && length(t) &&
        length(unique(lengths(t))) == 1L && lengths(t)[[1L]] >= 4L) {
        t <- as.data.frame(do.call(rbind, t), row.names = NULL)
        r$table <- t
    }
    if (is.data.frame(t) && ncol(t) >= 4L && fmt %in% c(2L, 3L, 4L)) {
        nm <- c("cycle", "index", "from", "to")
        if (ncol(t) > 4L)
            nm <- c(nm, paste0("in", seq_len(ncol(t) - 4L)))
        r$table <- .name_cols(t, nm)
    }
    r
}

#' @rdname tda_g_analyses
#' @export
tda_g_cutpoints <- function(g, per_node = FALSE, ...) {
    opts <- list(...)
    if (isTRUE(per_node)) opts$opt <- 2
    tda_g(g, "gcut", opts)
}

#' @rdname tda_g_analyses
#' @export
tda_g_flow <- function(g, matrix = FALSE, ...) {
    opts <- list(...)
    if (isTRUE(matrix)) opts$opt <- 2
    tda_g(g, "gflow", opts)
}

#' @rdname tda_g_analyses
#' @export
tda_g_degrees <- function(g) {
    # gni is TDA's node-information command; its columns are index,
    # node, in-degree, out-degree, loops (verified on directed and
    # undirected probes).  gni's remaining options are print formats,
    # which the parsed return makes irrelevant, so no further
    # arguments are accepted.
    r <- tda_g(g, "gni", list())
    d <- r$table
    # index, node, then in-degree, out-degree and loops once per value
    # variable of the graph (gd1.dat's I,J,V1,V2 gives two sets); a
    # loop is counted in its own column, not in the degrees -- the
    # manual's Box 1 of 7.2.1.1 (2002) counted it in both, TDA 6.4p
    # does not
    k <- (ncol(d) - 2L) %/% 3L
    nm <- c("index", "node")
    for (i in seq_len(max(k, 1L)))
        nm <- c(nm, paste0(c("in_degree", "out_degree", "loops"),
                           if (k > 1L) i else ""))
    names(d) <- nm[seq_len(ncol(d))]
    ind <- d[[3L]]; outd <- d[[4L]]
    # an undirected graph reports each incident edge on both sides,
    # so in == out == the degree; summing would double it
    d$degree <- if (isFALSE(g$directed)) ind else ind + outd
    attr(d, "run") <- r
    d
}

#' @rdname tda_g_analyses
#' @export
tda_g_edges <- function(g, format = c("edges", "lower_triangle",
                                       "lower_triangle_diag", "full",
                                       "full_square"), ...) {
    opts <- list(...)
    if (!missing(format))
        opts$opt <- match(match.arg(format), c("edges", "lower_triangle",
                                               "lower_triangle_diag",
                                               "full", "full_square"))
    tda_g(g, "gdp", opts)
}

#' @rdname tda_g_analyses
#' @export
tda_g_transitive <- function(g, ...) tda_g(g, "gtcl", list(...))

# The rest of TDA's graph family.  Each is the same shape as the twelve
# above: the command name, whatever options the caller passes, and the table
# TDA writes.  What differs is what each needs of the graph, which is noted
# where it matters -- TDA refuses the command otherwise rather than giving a
# wrong answer.

#' @rdname tda_g_analyses
#' @export
tda_g_aggregate <- function(g, ...) tda_g(g, "gda", list(...))

#' @rdname tda_g_analyses
#' @export
tda_g_symmetric <- function(g, ...) tda_g(g, "gsym", list(...))

#' @rdname tda_g_analyses
#' @export
tda_g_toposort <- function(g, relabel_graph = FALSE, ...) {
    opts <- list(...)
    if (isTRUE(relabel_graph)) opts$opt <- 2
    r <- tda_g(g, "gsort", opts)
    # manual 7.2.2.1 Box 2 names both layouts: option 1 the new label of
    # each node, option 2 the graph's edges relabelled, with values
    t <- r$table
    if (is.data.frame(t)) {
        if (ncol(t) == 3L)
            r$table <- .name_cols(t, c("i", "node_i", "label_i"))
        else if (ncol(t) == 5L)
            r$table <- .name_cols(t, c("label_i", "label_j", "node_i",
                                       "node_j", "value"))
    }
    r
}

#' @rdname tda_g_analyses
#' @export
tda_g_nodes <- function(g, ...) tda_g(g, "gni", list(...))

#' @rdname tda_g_analyses
#' @export
tda_g_dot <- function(g, ...) tda_g(g, "gdot", list(...))

#' @rdname tda_g_analyses
#' @export
tda_g_links <- function(g, backward = FALSE, ...) {
    opts <- list(...)
    if (isTRUE(backward)) opts$opt <- 2
    r <- tda_g(g, "gdln", opts)
    # gdln's records are ragged -- index, node, number of links, then
    # the linked nodes -- and come back as one row per node with the
    # links as a list column, the manual's Box 1 of 7.2.1.2 as a frame
    rows <- r$table
    if (inherits(rows, "tda_ragged") || (is.list(rows) && !is.data.frame(rows))) {
        m <- lapply(unclass(rows), as.numeric)
        r$table <- data.frame(
            index = vapply(m, function(v) as.integer(v[1L]), 1L),
            node = vapply(m, function(v) as.integer(v[2L]), 1L),
            n_links = vapply(m, function(v) as.integer(v[3L]), 1L),
            links = I(lapply(m, function(v) as.integer(v[-(1:3)]))),
            row.names = NULL)
    }
    r
}

#' @rdname tda_g_analyses
#' @export
tda_g_dcycles <- function(g, format = c("per_cycle", "per_node",
                                       "progress"), stop_at_first = FALSE,
                          ...) {
    opts <- list(...)
    if (!missing(format))
        opts$opt <- match(match.arg(format), c("per_cycle", "per_node",
                                                "progress"))
    if (isTRUE(stop_at_first)) opts$dopt <- 1
    tda_g(g, "gdcyc", opts)
}

#' @rdname tda_g_analyses
#' @export
tda_g_dblocks <- function(g, ...) tda_g(g, "gdcset", list(...))

#' @rdname tda_g_analyses
#' @export
tda_g_compact <- function(g, per_node = FALSE, ...) {
    opts <- list(...)
    if (isTRUE(per_node)) opts$opt <- 2
    tda_g(g, "gcset", opts)
}

#' @rdname tda_g_analyses
#' @export
tda_g_gcliques <- function(g, ...) tda_g(g, "ggcliq", list(...))

#' @rdname tda_g_analyses
#' @export
tda_g_random <- function(g, ...) tda_g(g, "gcd", list(...))

#' @rdname tda_g_analyses
#' @export
tda_g_ownership <- function(g, ...) {
    r <- tda_g(g, "gio", list(...))
    # manual 7.6.1.2 Box 6 names them: the pair, the direct share
    # A(i,j), the share including indirect holdings Y(i,j), the number
    # of nodes in the component, and the iterations taken
    if (is.data.frame(r$table) && ncol(r$table) == 7L)
        r$table <- .name_cols(r$table, c("i", "node_i", "node_j", "a_ij",
                                         "y_ij", "c", "iterations"))
    r
}

#' @rdname tda_g_analyses
#' @export
tda_g_backward <- function(g, ...) {
    opts <- list(...)
    # gbcf writes two files: the right-hand side takes the pair table
    # (one row per controlling/controlled pair), and df= a summary with
    # one row per node.  Without df= the summary is not written at all,
    # so it is asked for by default and returned as $summary.
    if (is.null(opts$df))
        opts$df <- "out.df"
    r <- tda_g(g, "gbcf", opts)
    # manual 7.6.1.1 Box 4 names the columns of both files
    if (is.data.frame(r$table) && ncol(r$table) == 8L)
        r$table <- .name_cols(r$table,
                              c("i", "node_i", "n_a", "node_j", "a_ji",
                                "n_y", "node_k", "y_ki"))
    p <- file.path(r$run$dir, opts$df)
    if (!is.null(r$run$dir) && nzchar(r$run$dir[1L]) && file.exists(p)) {
        d <- tryCatch(utils::read.table(p, header = FALSE),
                      error = function(e) NULL)
        if (is.data.frame(d) && ncol(d) == 3L)
            r$summary <- .name_cols(d, c("i", "node_i", "m"))
        else
            r$summary <- d
    }
    r
}

#' @rdname tda_g_analyses
#' @export
tda_g_paths <- function(g, ...) tda_g(g, "gep", list(...))

#' @rdname tda_g_analyses
#' @export
tda_g_flowcontrol <- function(g, ...) tda_g(g, "gfc", list(...))

#' @rdname tda_g_analyses
#' @export
tda_g_centred <- function(g, ...) tda_g(g, "gnc", list(...))

#' @rdname tda_g_analyses
#' @export
tda_g_spantrees <- function(g, edges = FALSE, ...) {
    opts <- list(...)
    if (isTRUE(edges)) opts$opt <- 2
    tda_g(g, "gnst", opts)
}

#' @rdname tda_g_analyses
#' @export
tda_g_reachable <- function(g, format = c("count", "list",
                                         "strong_nodes", "strong_edges"),
                            ...) {
    opts <- list(...)
    if (!missing(format))
        opts$opt <- match(match.arg(format), c("count", "list",
                                                "strong_nodes",
                                                "strong_edges"))
    tda_g(g, "gdcon", opts)
}

# gcni needs an unvalued graph.  tda_graph() always builds a valued one --
# gt = 2 and a V1 column, defaulting the values to 1 -- so the block is
# rebuilt here as gt = 1 with only the two node columns, rather than leaving
# the caller to work out why TDA refused.
#' @rdname tda_g_analyses
#' @export
tda_g_neighbourhoods <- function(g, ...) {
    v <- g$data$V1
    if (!is.null(v) && length(unique(v)) > 1L)
        stop("gcni needs an unvalued graph, but the edges carry ",
             length(unique(v)), " distinct values")
    # gdd's edge-list option always wants three columns -- the two nodes and
    # a graph indicator -- so V1 stays; gt is what makes the graph unvalued.
    g$opts$gt <- 1
    tda_g(g, "gcni", list(...))
}

#' @rdname tda_g_analyses
#' @export
tda_g_independent <- function(g, min_size = NULL, use_complement = FALSE,
                              ...) {
    opts <- list(...)
    if (!is.null(min_size)) opts$min <- min_size
    if (isTRUE(use_complement)) opts$cg <- 1
    tda_g(g, "giset", opts)
}

#' @rdname tda_g_analyses
#' @export
tda_g_spanning <- function(g) {
    # TDA's gst builds the MAXIMUM spanning tree of each component --
    # the strongest-ties convention for similarity-valued networks
    # (verified against igraph::mst on negated weights).  gst's
    # remaining options are print formats, which the parsed return
    # makes irrelevant, so no further arguments are accepted.
    tda_g(g, "gst", list())
}

#' @rdname tda_g_analyses
#' @export
tda_g_eigen <- function(g, n = NULL, tol = NULL, control = NULL, ...) {
    opts <- list(...)
    if (!is.null(n)) opts$n <- n
    if (!is.null(tol)) opts$eps <- tol
    opts <- c(opts, .control_opts(control))
    tda_g(g, "gev", opts)
}

#' @rdname tdaR-methods
#' @keywords internal
#' @exportS3Method base::print
print.tda_graph <- function(x, ...) {
    d <- x$data
    if (x$form == "edges" && ncol(d) >= 2L) {
        nodes <- unique(c(d[[1L]], d[[2L]]))
        cat("TDA graph: ", length(nodes), " nodes, ", nrow(d),
            if (x$directed) " directed" else " undirected",
            " edge", if (nrow(d) != 1L) "s",
            if (ncol(d) >= 3L) ", valued", "\n", sep = "")
        cat("\n")
        print(utils::head(d, 6L))
        if (nrow(d) > 6L)
            cat("  ...", nrow(d) - 6L, "more edges\n")
    } else {
        # adjacency / dissimilarity form: the data are the matrix rows
        cat("TDA graph: ", nrow(d), " nodes,",
            if (x$directed) " directed, " else " undirected, ",
            x$form, " form\n", sep = "")
        cat("\n")
        print(utils::head(d, 6L))
        if (nrow(d) > 6L)
            cat("  ...", nrow(d) - 6L, "more rows\n")
    }
    invisible(x)
}


#' Draw a graph
#'
#' TDA's graph-drawing command, \code{pltree}, only draws trees. For
#' anything else -- a cycle, more than one component, or a directed graph --
#' this lays the nodes out on a circle and draws it with
#' \code{\link{tda_pl_graph}} instead, which takes explicit coordinates and
#' so can draw any graph.
#'
#' @param x a \code{\link{tda_graph}}.
#' @param width,height size of the plotting area in millimetres.
#' @param file name for the PostScript file.
#' @param newpage start a new page before drawing.
#' @param layout \code{"auto"} (the default) uses \code{pltree} when the
#'   graph is actually a tree (connected, one fewer edge than nodes) and a
#'   circular layout otherwise; \code{"circle"} always uses the circular
#'   layout, \code{"tree"} always tries \code{pltree} and lets TDA refuse a
#'   non-tree graph with its own error.
#' @param ... further options for \code{pltree} (\code{layout = "tree"}) or
#'   \code{\link{tda_pl_graph}} (\code{layout = "circle"}).
#' @return The parsed drawing operations, invisibly.
#' @family graph analysis
#' @examples
#' e <- data.frame(from = c(1, 1, 2), to = c(2, 3, 4), value = 1)
#' g <- tda_graph(e, directed = FALSE)  # 4 nodes, 3 edges: already a tree
#' pdf(NULL); plot(g); dev.off()
#'
#' # a graph with a cycle is not a tree, and this one is directed besides --
#' # both are fine with the default circular layout
#' e2 <- data.frame(from = c(1, 1, 2, 2, 3, 4), to = c(2, 3, 3, 4, 4, 5))
#' g2 <- tda_graph(e2, directed = TRUE)
#' pdf(NULL); plot(g2); dev.off()
#' @exportS3Method base::plot
plot.tda_graph <- function(x, width = 110, height = 80, file = "graph.ps",
                           newpage = TRUE,
                           layout = c("auto", "circle", "tree"), ...) {
    layout <- match.arg(layout)
    e <- x$data
    n_nodes_v <- unique(c(e[[1L]], e[[2L]]))
    is_tree <- !x$directed && nrow(e) == length(n_nodes_v) - 1L
    if (layout == "tree" || (layout == "auto" && is_tree))
        return(.plot_tree(x, width, height, file, newpage, ...))
    .plot_graph_circle(x, width, height, newpage, ...)
}

# pltree: TDA's tree layout. Used when the graph is (or is asserted to
# be) a tree; TDA itself refuses anything else with "current graph is not
# a tree" if layout = "tree" was forced on a non-tree graph.
.plot_tree <- function(x, width, height, file, newpage, ...) {
    dir <- tempfile("tda")
    res <- tda_run(c(tda_nvar(x$data),
                     do.call(tda_block, c(list(name = "gdd"), x$opts,
                                          list(rhs = x$rhs))),
                     sprintf("psfile = %s;", file),
                     # pltree needs a logical coordinate system, so pxa and
                     # pya have to be given even though the layout chooses
                     # its positions within them.
                     tda_block("psetup", pxlen = width, pylen = height,
                               pxa = "0,100", pya = "0,100"),
                     do.call(tda_block, c(list(name = "pltree"), .tda_extra(list(...)))),
                     "psclose;"),
                   data = x$data, dir = dir)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err)) {
        # TDA's pltree specifically requires an undirected graph
        # (gt=2 on gdd) -- a tree's parent/child structure comes
        # from which nodes are connected at all, not from which way
        # any edge points, and TDA's tree-layout algorithm simply
        # does not accept a directed one, checked. This is
        # a real constraint of TDA's pltree command, not
        # something this package adds on top of it -- but TDA's
        # error for it ("graph must be undirected") does not say why
        # or what to do about it, so a directed graph (tda_graph()'s
        # own default) built for a tree layout is given a clearer
        # pointer back to directed = FALSE here instead of just
        # relaying that raw message.
        if (isTRUE(x$directed) && any(grepl("undirected", err)))
            stop("TDA's pltree only draws undirected trees: build ",
                "this graph with tda_graph(..., directed = FALSE) ",
                "instead (a tree's layout comes from which nodes ",
                "are connected, not which way an edge points, so ",
                "direction is not meaningful here regardless).",
                call. = FALSE)
        stop("TDA could not draw this graph: ", err[1L], call. = FALSE)
    }
    tda_plot_ps(tda_read_ps(res, which = file), newpage = newpage)
}

# A general-purpose fallback for any graph pltree cannot draw: place the
# nodes evenly around a circle and draw it with tda_pl_graph, which takes
# explicit coordinates rather than laying a graph out itself.  A circular
# layout is not the prettiest possible drawing, but unlike a distance-based
# one (multidimensional scaling on shortest-path distances, say) it never
# depends on the graph being connected, undirected, or of a particular
# shape, so it always produces something.
.plot_graph_circle <- function(x, width, height, newpage, ...) {
    e <- x$data
    ids <- sort(unique(c(e[[1L]], e[[2L]])))
    n <- length(ids)
    theta <- if (n > 1L) seq(0, 2 * pi, length.out = n + 1L)[seq_len(n)] else 0
    nodes <- data.frame(id = ids, x = 50 + 40 * cos(theta),
                        y = 50 + 40 * sin(theta))
    edges <- data.frame(from = e[[1L]], to = e[[2L]])
    if (x$directed)
        edges$arrow <- "1.5,1.0"
    p <- tda_ps(xlim = c(0, 100), ylim = c(0, 100), width = width,
               height = height)
    p <- tda_pl_graph(p, nodes, edges, ...)
    p <- plot(p, newpage = newpage)
    invisible(tda_read_ps(p$run, which = p$file))
}

#' Convert a distance matrix into a metric
#'
#' \code{dmet}: an undirected valued graph's edge values, read as a
#' distance matrix, do not always satisfy the triangle inequality (a
#' direct edge can be "longer" than going via a third node). This finds
#' the smallest constant that, added to every distance, makes the whole
#' matrix satisfy it -- \eqn{c = \max(0, \max_{i,j,k} d_{ij} - d_{ik} -
#' d_{kj})} -- and returns both the original and the corrected matrix.
#'
#' @param g a \code{\link{tda_graph}}, undirected and valued.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A list: \code{original} and \code{modified}, both distance
#'   matrices, and \code{constant}, the value added to get from one to
#'   the other.
#' @family graph analysis
#' @examples
#' e <- data.frame(from = c(1, 2, 3), to = c(2, 3, 1), value = c(1, 2, 5))
#' g <- tda_graph(e, directed = FALSE)   # 1-3 direct (5) exceeds
#'                                       # 1-2-3 via node 2 (1+2=3)
#' tda_dmet(g)
#' @export
tda_dmet <- function(g, options = list(), dir = tempfile("tda"), ...) {
    stopifnot(inherits(g, "tda_graph"))
    opts <- c(list(fmt = "18.12"), .tda_extra(options))
    res <- tda_run(c(tda_nvar(g$data),
                     do.call(tda_block, c(list(name = "gdd"), g$opts,
                                          list(rhs = g$rhs))),
                     do.call(tda_block, c(list(name = "dmet"), opts,
                                          list(rhs = "out.txt")))),
                   data = g$data, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run dmet: ", err[1L], call. = FALSE)
    # export first: dmet writes both matrices and reports the additive
    # constant on the console, and all three are exported, so out.txt
    # is only the tdaR.use_exports = FALSE path.
    eo <- if (.use_exports()) res$exports[["dmet.original"]]
    em <- if (.use_exports()) res$exports[["dmet.modified"]]
    ec <- if (.use_exports()) res$exports[["dmet.constant"]]
    # the same SHAPE the parser path returns -- a plain three-element
    # list, no class, no run -- so the two are interchangeable
    if (is.matrix(eo) && is.matrix(em) && is.matrix(ec))
        return(list(original = eo, modified = em,
                    constant = as.vector(ec)[1L]))
    p <- file.path(res$dir, "out.txt")
    if (!file.exists(p))
        stop("dmet produced no output")
    lines <- readLines(p, warn = FALSE)
    i1 <- grep("^Original", lines)
    i2 <- grep("^Modified", lines)
    if (!length(i1) || !length(i2))
        stop("dmet's output was not in the expected format")
    read_block <- function(from, to) {
        body <- lines[(from + 1L):to]
        do.call(rbind, lapply(strsplit(trimws(body), "\\s+"), as.numeric))
    }
    orig <- read_block(i1, i2 - 1L)
    mod <- read_block(i2, length(lines))
    const <- suppressWarnings(as.numeric(regmatches(res$output,
        regexpr("(?<=Additive constant: )[0-9.eE+-]+", res$output,
               perl = TRUE))))
    list(original = orig, modified = mod,
        constant = if (length(const)) const[1L] else NA_real_)
}

#' Print a tree's data structure
#'
#' \code{ptree}: the same tree \code{\link{plot.tda_graph}} draws (via
#' \code{pltree}), but as data rather than a picture -- one row per node,
#' its parent, its first child, and its next sibling, the linked-list
#' representation TDA itself builds internally to draw one.
#'
#' TDA's \code{gn=} (which graph, when more than one is open) is not
#' optional here in practice: \code{ptree}, unlike \code{pltree}, never
#' calls the setup step that gives it a safe default, so leaving it out
#' reliably fails with \dQuote{current graph is not a tree} even for a
#' graph that unambiguously is one (\code{pltree} accepts the identical
#' graph) -- \code{gn = 1} is always sent explicitly to work around
#' it.
#'
#' @param g a \code{\link{tda_graph}}, undirected and shaped like a tree.
#' @param root which node to treat as the root; the first node in
#'   \code{g}'s data by default.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A data frame: \code{node}, \code{edge_from}, \code{edge_to},
#'   \code{value}, \code{parent}, \code{first_child}, \code{next_sibling}
#'   (\code{0} where a node has none).
#' @family graph analysis
#' @examples
#' e <- data.frame(from = c(1, 1, 2), to = c(2, 3, 4), value = 1)
#' g <- tda_graph(e, directed = FALSE)   # 4 nodes, 3 edges: a tree
#' tda_ptree(g, root = 1)
#' @export
tda_ptree <- function(g, root = NULL, options = list(), dir = tempfile("tda"),
                      ...) {
    stopifnot(inherits(g, "tda_graph"))
    opts <- list(gn = 1)
    if (!is.null(root))
        opts$rt <- root
    opts <- c(opts, .tda_extra(options))
    res <- tda_run(c(tda_nvar(g$data),
                     do.call(tda_block, c(list(name = "gdd"), g$opts,
                                          list(rhs = g$rhs))),
                     do.call(tda_block, c(list(name = "ptree"), opts,
                                          list(rhs = "out.txt")))),
                   data = g$data, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run ptree: ", err[1L], call. = FALSE)
    etab0 <- if (.use_exports()) res$exports[["ptree.table"]]
    tab <- if (is.matrix(etab0)) .export_frame(etab0)
           else tda_file(res, "out.txt")
    if (is.null(tab))
        stop("ptree produced no output")
    names(tab) <- c("index", "node", "edge_from", "edge_to", "value",
                    "parent", "first_child", "next_sibling")
    tab$index <- NULL
    tab
}

#' Split an ordering into contiguous clusters
#'
#' \code{tda_ucl} splits an ordering of a graph's nodes into
#' \code{clusters} contiguous groups so that the largest cluster
#' diameter is as small as possible, by the dynamic program of Alpert
#' and Kahng (1997). The diameter of a cluster is the largest edge value
#' between any two of its nodes.
#'
#' The clusters are runs of the \emph{ordering}, not arbitrary subsets:
#' this answers "where do I cut this sequence" rather than "which things
#' group together". Give the ordering with \code{order}; the default is
#' the nodes in their internal order, which is rarely what you want
#' unless the graph was built that way. A seriation or a first principal
#' coordinate is the usual source.
#'
#' The graph must be undirected and valued, and no edge value may be
#' negative or missing.
#'
#' @section A correctness note:
#' TDA 6.4's \code{ucl} did not work. It optimised over the identity
#' ordering with its first two nodes swapped, counted a state with no
#' legal split as costing nothing, searched only the first few split
#' positions, and read past the end of its own cost matrix. On the
#' example shipped with TDA it returned a partition of diameter 7 where
#' one of 4 exists. The number of clusters and the size bounds were
#' fixed at 3, 1 and 3 whatever you asked for.
#'
#' All of that is repaired here, and the result is checked in this
#' package's tests against a brute force over every contiguous
#' partition. The fixes are ours, not Rohwer's, and are listed in the
#' package's \code{README}.
#'
#' @param g a \code{\link{tda_graph}}, undirected and valued.
#' @param clusters number of clusters, at least 2 (default 2).
#' @param min_size,max_size smallest and largest number of nodes a
#'   cluster may have. Defaults: 1, and the number of nodes.
#' @param order the ordering to split, as node indices \code{1..n} in
#'   the sequence they should be considered. \code{NULL}, the default,
#'   uses the graph's node order.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A data frame with one row per node: \code{position} in the
#'   ordering, \code{node} (internal index), \code{cluster}, and
#'   \code{label} (the node number the graph was built with).
#'   \code{attr(x, "diameter")} is the largest cluster diameter, the
#'   quantity minimised, and \code{attr(x, "run")} carries the run.
#' @references Alpert, C. J. and Kahng, A. B. (1997). Splitting an
#'   ordering into a partition to minimize diameter.
#'   \emph{Journal of Classification} \strong{14}, 51--74.
#' @family graph analysis
#' @examples
#' # six points on a line, distance growing with the gap between them
#' d <- as.dist(outer(1:6, 1:6, function(i, j) 3 * abs(i - j) + 1))
#' g <- tda_graph(d)
#' p <- tda_ucl(g, clusters = 3)
#' p
#' attr(p, "diameter")   # the largest cluster diameter, minimised
#' @export
tda_ucl <- function(g, clusters = 2, min_size = NULL, max_size = NULL,
                    order = NULL, options = list(), dir = tempfile("tda"),
                    ...) {
    stopifnot(inherits(g, "tda_graph"))
    if (!is.numeric(clusters) || length(clusters) != 1L || clusters < 2)
        stop("`clusters` must be a single number, 2 or more", call. = FALSE)
    o <- .tda_extra(options)
    o$nc <- as.integer(clusters)
    if (!is.null(min_size)) o$min <- as.integer(min_size)
    if (!is.null(max_size)) o$max <- as.integer(max_size)
    if (!is.null(order)) {
        if (anyNA(order) || !is.numeric(order) ||
            !setequal(order, seq_along(order)))
            stop("`order` must be a permutation of 1..n", call. = FALSE)
        o$cn <- paste(as.integer(order), collapse = ",")
    }
    res <- tda_g(g, "ucl", o, dir = dir, ...)
    tab <- .ucl_frame(res)
    if (is.null(tab))
        stop("ucl produced no partition", call. = FALSE)
    tab
}

# ucl writes one record per node and exports the same four columns as
# ucl.part, so the frame is built from the export and only falls back to
# the printed file if that is missing.
.ucl_frame <- function(res) {
    r <- attr(res, "run")
    if (is.null(r)) r <- res$run
    m <- if (!is.null(r) && .use_exports()) r$exports[["ucl.part"]] else NULL
    tab <- if (is.matrix(m) && ncol(m) == 4L)
        as.data.frame(m) else if (is.data.frame(res) && ncol(res) == 4L)
        res else NULL
    if (is.null(tab))
        return(NULL)
    names(tab) <- c("position", "node", "cluster", "label")
    if (!is.null(r)) {
        # the tap holds the double behind the diameter line; the printed
        # text is the fallback, and is rounded
        dv <- .tap_values(r, "Maximal cluster diameter", 1L)
        if (is.null(dv)) {
            d <- grep("Maximal cluster diameter", r$output, value = TRUE)
            if (length(d))
                dv <- as.numeric(sub(".*: ", "", d[1L]))
        }
        if (!is.null(dv))
            attr(tab, "diameter") <- dv
        attr(tab, "run") <- r
    }
    tab
}

#' One-dimensional multifacility location
#'
#' Place k new (variable) points among m existing (fixed) points on a
#' line so that the weighted interaction costs are minimal -- TDA's
#' \code{gloc}. Each variable point ends up at one of the fixed points;
#' only the order of the fixed points on the line matters, not their
#' spacing (verified against exhaustive search in the test suite).
#'
#' @param fixed_weights a k x m matrix: interaction weight between each
#'   variable point (row) and each fixed point (column, in line order).
#' @param variable_weights optional symmetric k x k matrix of weights
#'   between the variable points themselves (default: none).
#' @param options a named list of further TDA options for \code{gloc}.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#'   The graph-setup options \code{opt=}, \code{gt=}, \code{perm=}
#'   and \code{sc=} pass through to the underlying edge load, as in
#'   \code{\link{tda_g}}.
#' @return An integer vector of length k: for each variable point, the
#'   index of the fixed point it is placed at.
#' @family graph analysis
#' @examples
#' # two facilities, four sites; facility 1 mostly serves site 2,
#' # facility 2 site 3, and the two facilities interact
#' vf <- rbind(c(8, 1, 1, 2), c(1, 1, 7, 3))
#' tda_locate_line(vf, variable_weights = rbind(c(0, 5), c(5, 0)))
#' @export
tda_locate_line <- function(fixed_weights, variable_weights = NULL,
                            options = list(), dir = tempfile("tda"),
                            ...) {
    vf <- as.matrix(fixed_weights)
    k <- nrow(vf); m <- ncol(vf)
    if (k < 1L || m < 2L)
        stop("`fixed_weights` needs at least one row (variable point) ",
             "and two columns (fixed points)", call. = FALSE)
    if (anyNA(vf) || any(vf < 0))
        stop("weights must be non-negative and complete", call. = FALSE)
    ed <- data.frame(from = integer(), to = integer(),
                     wvv = numeric(), wvf = numeric())
    if (!is.null(variable_weights)) {
        vv <- as.matrix(variable_weights)
        if (!all(dim(vv) == k) || any(abs(vv - t(vv)) > 1e-12))
            stop("`variable_weights` must be a symmetric k x k matrix",
                 call. = FALSE)
        for (i in seq_len(k - 1L)) for (j in (i + 1L):k)
            if (vv[i, j] > 0)
                ed <- rbind(ed, data.frame(from = i, to = j,
                                           wvv = vv[i, j], wvf = 0))
    }
    for (i in seq_len(k)) for (j in seq_len(m))
        if (vf[i, j] > 0)
            ed <- rbind(ed, data.frame(from = i, to = k + j,
                                       wvv = 0, wvf = vf[i, j]))
    if (!nrow(ed))
        stop("all weights are zero: nothing constrains the locations",
             call. = FALSE)
    # every node must appear in the edge list so the graph has k+m nodes;
    # tie any isolated fixed point in with a zero-weight edge
    seen <- unique(c(ed$from, ed$to))
    for (j in setdiff(seq_len(k + m), seen))
        ed <- rbind(ed, data.frame(from = min(seen), to = j,
                                   wvv = 0, wvf = 0))
    names(ed) <- c("FROM", "TO", "WVV", "WVF")
    res <- tda_run(c(tda_nvar(ed),
                     tda_block(name = "gdd", opt = 1,
                               rhs = "FROM,TO,WVV,WVF"),
                     do.call(tda_block,
                             c(list(name = "gloc", gn = "1,2", n = k),
                               .tda_extra(options)))),
                   data = ed, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not place the points: ", err[1L], call. = FALSE)
    em <- if (.use_exports()) res$exports[["gloc.assignment"]]
    if (is.matrix(em) && ncol(em) == 3L)
        return(as.integer(em[, 2L]))
    ln <- grep("^Variable point [0-9]+ -> fixed point", res$output,
               value = TRUE)
    if (length(ln) != k)
        stop("gloc produced no readable locations", call. = FALSE)
    as.integer(sub(".*fixed point ([0-9]+).*", "\\1", ln))
}
