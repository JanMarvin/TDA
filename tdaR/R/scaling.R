# Measurement and scaling wrappers: MDS,
# relation fitting, scalogram, conjoint, clustering, combinatorics,
# rank orders, Rasch, unfolding, Procrustes.

.tda_valued_graph <- function(d, cmds, ...) {
    d <- as.matrix(d)
    ut <- do.call(rbind, lapply(seq_len(nrow(d) - 1), function(i)
        cbind(i, (i + 1):nrow(d), d[i, (i + 1):nrow(d)])))
    .tda_graph_run(ut, cmds, directed = FALSE, valued = TRUE, ...)
}

#' Nonmetric multidimensional scaling
#'
#' TDA's \code{mdsn1} (Kruskal) on a distance matrix.
#'
#' @param d symmetric distance matrix.
#' @param ties tie handling. Only \code{"ignore"} is implemented in
#'   \code{mdsn1} (its primary and secondary projections are commented
#'   out in the C, and the option would silently do nothing); the other
#'   two values error and point to \code{\link{tda_mds}}, whose
#'   \code{mdsn} has them.
#' @param max_iterations iteration cap.
#' @param restarts random starting configurations, TDA's \code{ns=}; the
#'   one with the lowest stress is returned. Nonmetric MDS on few points
#'   has degenerate local minima (clusters of near-coincident points with
#'   small stress), and restarts are the guard against them.
#' @param ... passed to \code{\link{tda_run}}.
#' @return with direct exports enabled, a list with \code{stress},
#'   \code{configuration} (n x 2) and \code{output}; otherwise the
#'   printed output lines.
#' @examples
#' d <- as.matrix(dist(cbind(c(0, 0, 3, 3), c(0, 2, 0, 2))))
#' tda_mdsn(d)
#' @export
tda_mdsn <- function(d, ties = c("ignore", "primary", "secondary"),
                     max_iterations = 100, restarts = 5, ...) {
    ties <- match.arg(ties)
    # mdsn1_proj() in t_mds.c implements opt=1 only; the primary and
    # secondary projections are in the source but commented out, and
    # opt=2/3 fall through to no projection at all
    if (ties != "ignore")
        stop("mdsn1 handles ties as \"ignore\" only; for Kruskal's primary ",
             "or secondary approach use tda_mds(method = \"nonmetric\", ",
             "options = list(opt = 2))  (or opt = 3)", call. = FALSE)
    d <- .tda_check_matrix(d, "d", square = TRUE, symmetric = TRUE)
    res <- .tda_valued_graph(d,
        sprintf("mdsn1(gn=1, opt=%d, mxit=%d, ns=%d) = m.out;",
                match(ties, c("ignore", "primary", "secondary")),
                max_iterations, restarts), ...)
    fv <- .file_values(res)
    st <- grep("Final best stress", res$output, value = TRUE)
    stress <- suppressWarnings(as.numeric(
        sub(".*: *([0-9.eE+-]+).*", "\\1", st[1L])))
    if (is.null(fv)) {                       # flag-off: same shape
        out <- file.path(res$dir, "m.out")   # from the file
        if (file.exists(out) && file.size(out) > 0)
            fv <- as.matrix(utils::read.table(out))
    }
    if (!is.null(fv) && is.matrix(fv) && nrow(fv) >= 1L) {
        # one row per restart: its number, its stress, then x2..xn and
        # y2..yn -- node 1 is pinned at the origin (2n-1 free parameters;
        # established from the output layout and mdsm's coordinate
        # scheme). The row with the lowest stress is the result.
        v <- fv[which.min(fv[, 2L]), -(1:2)]
        n <- (length(v) + 2L) %/% 2L
        cfg <- rbind(c(0, 0), cbind(v[1:(n - 1)], v[n:(2L * n - 2L)]))
        return(.tda_structured(
            list(stress = stress, configuration = cfg,
                 output = res$output), "tda_mds_nonmetric"))
    }
    res$output
}

#' Axes in an MDS configuration
#'
#' TDA's \code{mdsr}: regresses external variables onto a
#' two-dimensional configuration to draw interpretable axes.
#'
#' @param config two-column matrix, the configuration.
#' @param x matrix of external variables (same row count).
#' @param ... passed to \code{\link{tda_run}}.
#' @return the printed output.
#' @examples
#' cfg <- cbind(c(0, 1, 2, 3), c(0, 1, 0, 1))
#' r <- tda_mdsr(cfg, x = cbind(v = c(1, 2, 3, 4)))
#' cat(head(tda_payload(r), 12), sep = "\n")
#' @export
tda_mdsr <- function(config, x, ...) {
    if (inherits(config, "tda_mds_nonmetric"))
        config <- config$configuration   # accept tda_mdsn's result
    config <- as.matrix(config); x <- as.matrix(x)
    if (ncol(config) != 2L)
        stop("`config` must be a two-column configuration ",
             "(or a tda_mdsn result)", call. = FALSE)
    if (nrow(x) != nrow(config))
        stop("`x` must have one row per configuration point",
             call. = FALSE)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(cbind(config, x), file.path(dr, "m.dat"),
                       row.names = FALSE, col.names = FALSE)
    vx <- sprintf("V%d[12.4]=c%d", seq_len(ncol(x)), 2L + seq_len(ncol(x)))
    res <- tda_run(c(sprintf("nvar(dfile=m.dat, X[12.4]=c1, Y[12.4]=c2, %s);",
                             paste(vx, collapse = ", ")),
                     sprintf("mdsr(xv=%s) = X,Y;",
                             paste(sprintf("V%d", seq_len(ncol(x))),
                                   collapse = ","))),
                   dir = dr, ...)
    res$output
}

#' Fit binary relations to a valued graph
#'
#' TDA's \code{rfit}: finds the relations with the requested
#' properties closest to the graph, via linear programming.
#'
#' @param d symmetric matrix (proximities).
#' @param properties character subset of "reflexive", "symmetric",
#'   "antisymmetric", "transitive", "complete".
#' @param max_solutions cap on the number of equally good relations TDA
#'   keeps (\code{max=}, 10). When more exist the fit is refused with an
#'   error naming this argument; a linear-order request
#'   (antisymmetric, transitive, complete) on a symmetric input has 17
#'   optimal relations for 4 nodes and 491 for 6. What comes back is
#'   correct but not always complete: against every relation on 4 nodes
#'   enumerated in R (test-gap-wrappers.R, and 30 directed cases in the
#'   session-52 notes) \code{rfit}'s value is always the optimum and every
#'   relation it lists is optimal with the properties, but in about one
#'   case in five with ties it lists a subset of the tied optima (2 of
#'   4, 4 of 5) -- a limit of the CACM 449 enumeration it uses.
#' @param variant "standard" runs rfit; "experimental" runs the
#'   rfit1 sibling (same syntax, alternative implementation).
#' @param ... passed to \code{\link{tda_run}}.
#' @return \code{relations}, one fitted relation matrix per solution,
#'   and \code{output}.
#' @examples
#' d <- rbind(c(0, 2, 1), c(2, 0, 2), c(1, 2, 0))
#' r <- tda_rfit(d)
#' r$relations[[1]]
#' @export
tda_rfit <- function(d, properties = c("symmetric", "transitive"),
                     max_solutions = 10, variant = c("standard",
                     "experimental"), ...) {
    variant <- match.arg(variant)
    d <- .tda_check_matrix(d, "d", square = TRUE, symmetric = TRUE)
    props <- match(match.arg(properties,
        c("reflexive", "symmetric", "antisymmetric", "transitive",
          "complete"), several.ok = TRUE),
        c("reflexive", "symmetric", "antisymmetric", "transitive",
          "complete"))
    res <- .tda_valued_graph(d,
        sprintf("%s(gn=1, rel=%s, max=%d) = r.out;",
                if (variant == "standard") "rfit" else "rfit1",
                paste(props, collapse = ","), max_solutions), ...)
    # r.out holds each fitted relation twice: a flattened line, then
    # its n rows; solutions are blank-separated, which the tap shows as
    # a gap in the line numbers
    # lpi() returns -2 when more optimal relations exist than max=
    # allows; rfit then prints this and writes nothing.  A relation set
    # like antisymmetric + transitive + complete (linear orders) on a
    # symmetric input has many equally good fits (17 for 4 nodes, 491 for
    # 6), so the default 10 is soon too small.
    if (any(grepl("exceeded maximal number of solutions", res$output)))
        stop("more than ", max_solutions, " equally good relations; ",
             "raise `max_solutions`", call. = FALSE)
    bl <- .file_blocks(res, "r.out")
    if (!length(bl)) return(res$output)
    sols <- lapply(bl, function(b) {
        if (length(b) < 2L) return(NULL)
        do.call(rbind, lapply(b[-1L], as.integer))
    })
    sols <- Filter(Negate(is.null), sols)
    names(sols) <- NULL
    .tda_structured(list(relations = sols, output = res$output),
                    "tda_rfit")
}

#' Simple Rasch model
#'
#' TDA's \code{rmod} on binary item response data.
#'
#' @param x matrix of 0/1 responses, persons in rows.
#' @param ... passed to \code{\link{tda_run}}.
#' @return centered item parameters in \code{$items} (joint ML is
#'   identified only up to a constant, and its item estimates carry
#'   the known upward bias of about k/(k-1) for k items), the
#'   pattern parameters in \code{$patterns}, and the maximised log
#'   likelihood in \code{$logLik}.
#' @examples
#' set.seed(1)
#' ab <- rnorm(80)
#' x <- sapply(c(-1, 1), function(d) rbinom(80, 1, plogis(ab - d)))
#' tda_rmod(x)$items   # easier item first, harder second
#' @export
tda_rmod <- function(x, ...) {
    x <- .tda_check_matrix(x, "x", binary = TRUE)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(x, file.path(dr, "r.dat"),
                       row.names = FALSE, col.names = FALSE)
    vn <- sprintf("X%d[4.0]=c%d", seq_len(ncol(x)), seq_len(ncol(x)))
    res <- tda_run(c(sprintf("nvar(dfile=r.dat, %s);", paste(vn, collapse = ", ")),
                     sprintf("rmod(m=1, prot=p.txt) = %s;",
                             paste(sprintf("X%d", seq_len(ncol(x))),
                                   collapse = ","))),
                   dir = dr, ...)
    k <- ncol(x)
    par <- res$exports[["rmod.par"]]
    par <- if (.use_exports() && is.matrix(par) && ncol(par) >= k)
        as.numeric(par[1L, ])
    else {
        # flag-off: the protocol file prints the final parameter
        # vector after the algorithm ends
        pf <- file.path(dr, "p.txt")
        pl <- if (file.exists(pf)) readLines(pf) else character()
        i <- grep("^Parameter vector", pl)
        if (length(i))
            suppressWarnings(as.numeric(strsplit(trimws(pl[i[1L] + 1L]),
                                                 " +")[[1L]]))
        else numeric()
    }
    if (length(par) < k)
        return(res$output)
    # joint (unconditional) ML: identified only up to a constant, so
    # the item parameters are reported centered; they carry the known
    # upward bias of joint ML (about k/(k-1) for k items)
    items <- par[seq_len(k)] - mean(par[seq_len(k)])
    names(items) <- colnames(x) %||% sprintf("item%d", seq_len(k))
    ll <- grep("Maximum of log likelihood:", res$output, value = TRUE)
    .tda_structured(list(
        items = items,
        patterns = if (length(par) > k) par[-seq_len(k)] else NULL,
        logLik = .tap_values(res, "Maximum of log likelihood:", 1L) %||%
            (if (length(ll))
                 suppressWarnings(as.numeric(sub(".*: *", "", ll[1L])))
             else NA_real_),
        output = res$output), "tda_rmod")
}

#' Unfolding of preference data
#'
#' TDA's \code{unf}: each row ranks the alternatives (larger is
#' preferred).
#'
#' @param x matrix of preference values.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the printed output.
#' @examples
#' r <- tda_unf(rbind(c(3, 2, 1), c(1, 3, 2), c(2, 3, 1)))
#' cat(tail(tda_payload(r), 4), sep = "\n")   # the best permutation
#' @export
tda_unf <- function(x, ...) {
    x <- as.matrix(x)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(x, file.path(dr, "u.dat"),
                       row.names = FALSE, col.names = FALSE)
    vn <- sprintf("X%d[8.2]=c%d", seq_len(ncol(x)), seq_len(ncol(x)))
    res <- tda_run(c(sprintf("nvar(dfile=u.dat, %s);", paste(vn, collapse = ", ")),
                     sprintf("unf = %s;",
                             paste(sprintf("X%d", seq_len(ncol(x))),
                                   collapse = ","))),
                   dir = dr, ...)
    res$output
}

#' Scalogram (Guttman) analysis
#'
#' TDA's \code{sga}: values > 0 mean success on a task.
#'
#' @param x matrix, persons in rows, tasks in columns.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the printed output.
#' @examples
#' r <- tda_sga(rbind(c(1, 1, 0), c(1, 0, 0), c(1, 1, 1)))
#' cat(head(tda_payload(r), 12), sep = "\n")
#' @export
tda_sga <- function(x, ...) {
    x <- as.matrix(x)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(x, file.path(dr, "s.dat"),
                       row.names = FALSE, col.names = FALSE)
    vn <- sprintf("X%d[4.0]=c%d", seq_len(ncol(x)), seq_len(ncol(x)))
    res <- tda_run(c(sprintf("nvar(dfile=s.dat, %s);", paste(vn, collapse = ", ")),
                     sprintf("sga = %s;",
                             paste(sprintf("X%d", seq_len(ncol(x))),
                                   collapse = ","))),
                   dir = dr, ...)
    res$output
}

#' Combinatorial patterns
#'
#' TDA's \code{com}: enumerate tuples, sets, permutations or
#' partitions into an output file.
#'
#' @param pattern one of "tuples", "subsets", "m_tuples",
#'   "permutations", "partitions", "partitions_m", "bipartitions".
#' @param n,m dimensions (n the base size, m where the pattern
#'   needs it).
#' @param ... passed to \code{\link{tda_run}}.
#' @return the output file's lines.
#' @examples
#' # the six 2-element subsets of {0, 1, 2, 3}: a row number, then the members
#' tda_com("subsets", n = 4, m = 2)
#' # the five partitions of 4: 4, 3+1, 2+2, 2+1+1, 1+1+1+1
#' tda_com("partitions", n = 4)
#' @export
tda_com <- function(pattern = c("tuples", "subsets", "m_tuples",
                    "permutations", "partitions", "partitions_m",
                    "bipartitions"), n = 1, m = 1, ...) {
    pattern <- match.arg(pattern)
    dr <- tempfile("tda"); dir.create(dr)
    res <- tda_run(sprintf("com(opt=%d, n=%d, m=%d) = c.out;",
                           match(pattern, c("tuples", "subsets", "m_tuples",
                                 "permutations", "partitions",
                                 "partitions_m", "bipartitions")), n, m),
                   dir = dr, ...)
    fv <- .file_values(res)
    if (!is.null(fv))
        return(fv)
    out <- file.path(dr, "c.out")
    { fl <- if (file.exists(out)) readLines(out) else character()
      if (length(fl)) fl else res$output }
}

#' Conjoint analysis
#'
#' TDA's \code{conj} (alternating least squares) on a preference
#' variable and factor columns.
#'
#' @param y preference values.
#' @param x factor matrix.
#' @param algorithm "monotone" regression or "lp" (linear
#'   programming).
#' @param ... passed to \code{\link{tda_run}}.
#' @return the printed output.
#' @examples
#' x <- expand.grid(a = 1:2, b = 1:2)
#' r <- tda_conj(y = c(4, 3, 2, 1), x = x)
#' cat(tail(tda_payload(r), 8), sep = "\n")   # the part-worths
#' @export
tda_conj <- function(y, x, algorithm = c("monotone", "lp"), ...) {
    algorithm <- match.arg(algorithm)
    x <- as.matrix(x)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(cbind(y, x), file.path(dr, "c.dat"),
                       row.names = FALSE, col.names = FALSE)
    vn <- sprintf("X%d[6.0]=c%d", seq_len(ncol(x)), 1L + seq_len(ncol(x)))
    res <- tda_run(c(sprintf("nvar(dfile=c.dat, Y[8.2]=c1, %s);",
                             paste(vn, collapse = ", ")),
                     sprintf("conj(alg=%d) = Y,%s;",
                             match(algorithm, c("monotone", "lp")),
                             paste(sprintf("X%d", seq_len(ncol(x))),
                                   collapse = ","))),
                   dir = dr, ...)
    res$output
}

#' Pyramidal clustering
#'
#' TDA's \code{clpyr} on a distance matrix.
#'
#' @param d symmetric distance matrix.
#' @param linkage "single" or "complete" link.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the printed output.
#' @examples
#' r <- tda_clpyr(as.matrix(dist(c(0, 1, 5, 6))))
#' # the top of the pyramid joins everything: 4 objects spanning 1..4
#' tail(r$pyramid, 1)
#' @export
tda_clpyr <- function(d, linkage = c("single", "complete"), ...) {
    linkage <- match.arg(linkage)
    d <- .tda_check_matrix(d, "d", square = TRUE, symmetric = TRUE)
    res <- .tda_valued_graph(d,
        sprintf("clpyr(opt=%d) = clpyr.out;",
                match(linkage, c("single", "complete"))), ...)
    m <- tda_file(res, "clpyr.out")
    names(m) <- c("cluster", "left", "right", "size", "min", "max",
                  "next_left", "next_right", "index")
    .tda_structured(list(pyramid = m, output = res$output), "tda_pyramid")
}

#' Context dependencies in Boolean functions
#'
#' TDA's \code{bfc}: Y as a Boolean function of binary X's.
#'
#' @param y binary outcome.
#' @param x binary matrix.
#' @param undefined how undefined arguments count: "dontcare",
#'   "as_one", or "as_zero".
#' @param ... passed to \code{\link{tda_run}}.
#' @return the printed output (implicants).
#' @examples
#' x <- expand.grid(a = 0:1, b = 0:1)
#' r <- tda_bfc(y = as.integer(x$a & x$b), x = x)
#' cat(head(tda_payload(r), 12), sep = "\n")
#' @export
tda_bfc <- function(y, x, undefined = c("dontcare", "as_one",
                    "as_zero"), ...) {
    undefined <- match.arg(undefined)
    x <- .tda_check_matrix(x, "x", binary = TRUE)
    if (!all(y %in% c(0, 1)))
        stop("`y` must contain only 0 and 1", call. = FALSE)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(cbind(y, x), file.path(dr, "b.dat"),
                       row.names = FALSE, col.names = FALSE)
    vn <- sprintf("X%d[4.0]=c%d", seq_len(ncol(x)), 1L + seq_len(ncol(x)))
    res <- tda_run(c(sprintf("nvar(dfile=b.dat, Y[4.0]=c1, %s);",
                             paste(vn, collapse = ", ")),
                     sprintf("bfc(opt=%d) = Y,%s;",
                             match(undefined, c("dontcare", "as_one", "as_zero")),
                             paste(sprintf("X%d", seq_len(ncol(x))),
                                   collapse = ","))),
                   dir = dr, ...)
    res$output
}

#' Procrustes rotation
#'
#' TDA's \code{mproc}: rotates configuration Y to best match X.
#' Verified in the suite against the closed-form solution from the
#' singular value decomposition.
#'
#' @param x,y n-by-m matrices.
#' @param ... passed to \code{\link{tda_run}}.
#' @return the rotated configuration as a matrix.
#' @examples
#' x <- cbind(c(0, 1, 0), c(0, 0, 1))
#' th <- 0.7; R <- rbind(c(cos(th), -sin(th)), c(sin(th), cos(th)))
#' round(tda_mproc(x, x %*% R), 6)
#' @export
tda_mproc <- function(x, y, ...) {
    x <- as.matrix(x); y <- as.matrix(y)
    stopifnot(all(dim(x) == dim(y)))
    dr <- tempfile("tda"); dir.create(dr)
    mdef <- function(nm, m)
        sprintf("mdef(%s,%d,%d) = %s;", nm, nrow(m), ncol(m),
                paste(t(m), collapse = ","))
    res <- tda_run(c(mdef("X", x), mdef("Y", y), "mfmt = 20.10;",
              "mproc(X,Y,Z);", "mpr(Z) = z.out;"), dir = dr, ...)
    m <- .mpr_from_exports(res)
    if (is.null(m))
        m <- as.matrix(tda_file(res, "z.out"))
    unname(m)
}

#' Rank-order data utilities
#'
#' TDA's \code{rod}: rank orders per record, with options for
#' standardization, distances, medians, and tabulation.
#'
#' @param x matrix of rank values, one order per row.
#' @param statistic one of "standard", "distances", "central",
#'   "graph", "graph5", "graph6", "between", "table".
#' @param ... passed to \code{\link{tda_run}}.
#' @return the rank orders as a data frame in \code{$orders}
#'   (plus a count column for the table statistic), with the full
#'   protocol in \code{$output}.
#' @examples
#' r <- tda_rod(rbind(c(1, 2, 3), c(2, 3, 1), c(3, 1, 2)))
#' r$orders
#' @export
tda_rod <- function(x, statistic = c("standard", "distances",
                    "central", "graph", "graph5", "graph6", "between",
                    "table"), ...) {
    statistic <- match.arg(statistic)
    x <- as.matrix(x)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(x, file.path(dr, "r.dat"),
                       row.names = FALSE, col.names = FALSE)
    vn <- sprintf("X%d[4.0]=c%d", seq_len(ncol(x)), seq_len(ncol(x)))
    res <- tda_run(c(sprintf("nvar(dfile=r.dat, %s);", paste(vn, collapse = ", ")),
                     sprintf("rod(opt=%d, df=rod.out) = %s;",
                             match(statistic, c("standard", "distances",
                                   "central", "graph", "graph5", "graph6",
                                   "between", "table")),
                             paste(sprintf("X%d", seq_len(ncol(x))),
                                   collapse = ","))),
                   dir = dr, ...)
    f <- file.path(dr, "rod.out")
    if (!file.exists(f))
        return(res$output)
    d <- tda_file(res, f)
    names(d)[seq_len(ncol(x))] <- sprintf("r%d", seq_len(ncol(x)))
    if (ncol(d) > ncol(x))
        names(d)[ncol(x) + 1L] <- "count"
    .tda_structured(list(orders = d, output = res$output), "tda_rod")
}

#' Rates in age-period form
#'
#' TDA's \code{rap}: builds age-period tables from starting time,
#' cohort, ending time and destination variables.
#'
#' @param cohort,entry_year,exit_year,destination integer
#'   vectors, one entry per case: birth year, year of entry
#'   into observation, year of exit, and the destination code
#'   (0 censored).  The engine requires
#'   cohort <= entry_year <= exit_year.
#' @param years,ages length-2 table ranges.
#' @param ... passed to \code{\link{tda_run}}.
#' @return \code{table}, the rate table TDA wrote, and \code{output}.
#' @examples
#' # three subjects born 1960/1961, observed 1980-1986: the risk
#' # table fills along the Lexis diagonals
#' r <- tda_rap(cohort = c(60, 61, 60), entry_year = c(80, 82, 83),
#'              exit_year = c(85, 86, 86), destination = c(1, 1, 0),
#'              years = c(80, 86), ages = c(18, 26))
#' r$risk
#' @export
tda_rap <- function(cohort, entry_year, exit_year, destination,
                    years, ages, ...) {
    if (any(cohort > entry_year) || any(entry_year > exit_year))
        stop("rap needs cohort <= entry_year <= exit_year per case",
             call. = FALSE)
    dr <- tempfile("tda"); dir.create(dr)
    utils::write.table(data.frame(cohort, entry_year, exit_year, destination),
                       file.path(dr, "a.dat"),
                       row.names = FALSE, col.names = FALSE)
    res <- tda_run(c("nvar(dfile=a.dat, TS[6.0]=c1, TC[6.0]=c2, TF[6.0]=c3, D[4.0]=c4);",
                     sprintf("rap(year=%d,%d, age=%d,%d, df=a.out) = TS,TC,TF,D;",
                             years[1], years[2], ages[1], ages[2])),
                   dir = dr, ...)
    # a.out holds three titled tables -- risk set, events, rates -- each
    # a year header row over one row per age; the titles leave gaps in
    # the tap's line numbers, which is where the blocks split
    bl <- .file_blocks(res, "a.out")
    {
        tab <- lapply(bl, function(b) {
            hdr <- b[[1L]]
            m <- do.call(rbind, b[-1L])
            if (ncol(m) == length(hdr) + 1L)
                colnames(m) <- c("age", paste0("y", hdr))
            m
        })
        nm <- c("risk", "events", "rates")
        names(tab) <- nm[seq_along(tab)]
        return(.tda_structured(c(tab, list(output = res$output, run = res)),
                               "tda_rap"))
    }
}
