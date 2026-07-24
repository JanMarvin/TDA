# sequence analysis
#


# ---- sequence analysis -----------------------------------------------------

# seqdef defines the sequence structure and seqm computes proximities between
# sequences by optimal matching.  TDA writes the result as a triple list --
# case i, case j, their lengths, and the distance -- which becomes a dist
# object, so hclust(), cmdscale() and the rest work on it directly.

# TDA's step-by-step dump of the alignment algorithm's internal
# dynamic-programming matrix (tst=2,3 + df=), one block per pair
# actually compared -- format checked against a real run's
# output (examples/exam/seqm3.tst): seq_pps() in t_seqm.c
# writes "Sequence A (case number ..)"/"Sequence B (case number ..)" and
# the two sequences themselves (tst=2), then seq_dmat() writes one
# "<D|E|F> Matrix" block per matrix the method uses (tst=3) -- D always,
# E and F too for method=2's affine costs. Row/column names on the
# matrices this returns are the actual state values at each position of
# the two sequences, prefixed with "A"/"B" for the empty-prefix
# row/column -- matching TDA's labels there exactly, not TDA's
# internal 0-based state rank.
# The same list, built from the exports: seqm.dp.pair carries each
# pair's two case numbers, seqm.dp.seqA / .seqB their state sequences
# (one block per pair), and seqm.dp.D/.E/.F the matrices themselves.
# Nothing here reads dp.tst.
.seq_dp_from_exports <- function(res) {
    if (is.null(res) || !.use_exports())
        return(NULL)
    pr <- res$exports[["seqm.dp.pair"]]
    if (!is.matrix(pr) || ncol(pr) != 2L)
        return(NULL)
    sa <- .exports_blocks(res$exports, "seqm.dp.seqA")
    sb <- .exports_blocks(res$exports, "seqm.dp.seqB")
    if (length(sa) != nrow(pr) || length(sb) != nrow(pr))
        return(NULL)
    out <- vector("list", nrow(pr))
    nms <- character(nrow(pr))
    for (k in seq_len(nrow(pr))) {
        seqA <- as.vector(sa[[k]])
        seqB <- as.vector(sb[[k]])
        mats <- list()
        for (typ in c("D", "E", "F")) {
            m <- .exkey_n(res, paste0("seqm.dp.", typ), k)
            if (!is.matrix(m))
                next
            if (!identical(dim(m),
                           c(length(seqA) + 1L, length(seqB) + 1L)))
                return(NULL)
            dimnames(m) <- list(c("A", as.character(seqA)),
                                c("B", as.character(seqB)))
            mats[[typ]] <- m
        }
        if (!length(mats))
            return(NULL)
        nms[k] <- paste0(pr[k, 1L], "_", pr[k, 2L])
        out[[k]] <- c(list(caseA = as.integer(pr[k, 1L]),
                           caseB = as.integer(pr[k, 2L]),
                           seqA = seqA, seqB = seqB), mats)
    }
    names(out) <- nms
    out
}

.seq_parse_dp <- function(path, res = NULL) {
    e <- .seq_dp_from_exports(res)
    if (!is.null(e))
        return(e)
    if (!file.exists(path))
        return(NULL)
    lines <- readLines(path, warn = FALSE)
    starts <- grep("^Sequence A ", lines)
    if (!length(starts))
        return(list())
    ends <- c(starts[-1L] - 1L, length(lines))
    out <- vector("list", length(starts))
    nms <- character(length(starts))
    for (k in seq_along(starts)) {
        blk <- lines[starts[k]:ends[k]]
        caseA <- as.integer(sub(".*case number (\\d+).*", "\\1", blk[1L]))
        seqA <- as.numeric(strsplit(trimws(blk[2L]), "\\s+")[[1L]])
        bline <- grep("^Sequence B ", blk)
        caseB <- as.integer(sub(".*case number (\\d+).*", "\\1", blk[bline]))
        seqB <- as.numeric(strsplit(trimws(blk[bline + 1L]), "\\s+")[[1L]])
        nrow_i <- length(seqA) + 1L
        mats <- list()
        mstarts <- grep("^[DEF] Matrix", blk)
        for (m in seq_along(mstarts)) {
            typ <- substr(blk[mstarts[m]], 1L, 1L)
            rows <- blk[(mstarts[m] + 4L):(mstarts[m] + 4L + nrow_i - 1L)]
            vals <- lapply(rows, function(r) {
                after_bar <- strsplit(r, "\\|")[[1L]][2L]
                as.numeric(strsplit(trimws(after_bar), "\\s+")[[1L]])
            })
            mat <- do.call(rbind, vals)
            # dp.tst prints these at %7.2f; the export carries the
            # float TDA actually held in SeqMD/SeqME/SeqMF, which is
            # what $distance is now taken from, so the two agree.
            ex <- if (!is.null(res))
                .exkey_n(res, paste0("seqm.dp.", typ), k)
            if (.use_exports() && is.matrix(ex) && all(dim(ex) == dim(mat)))
                mat <- ex
            rownames(mat) <- c("A", as.character(seqA))
            colnames(mat) <- c("B", as.character(seqB))
            mats[[typ]] <- mat
        }
        nms[k] <- paste0(caseA, "_", caseB)
        out[[k]] <- c(list(caseA = caseA, caseB = caseB, seqA = seqA,
                           seqB = seqB), mats)
    }
    names(out) <- nms
    out
}

#' Sequence analysis by optimal matching
#'
#' Computes proximities between sequences and returns a \code{dist}, so
#' \code{hclust}, \code{cmdscale} and the rest work on the result directly.
#'
#' @section Costs:
#' With the default costs a substitution costs twice an indel, so it never
#' beats a delete plus an insert, and the result is the indel distance
#' \code{|a| + |b| - 2 * LCS(a, b)}. Setting \code{indel = 2} makes
#' substitution the cheaper move and the result becomes twice the Levenshtein
#' distance.
#'
#' @param sequences one row per case, one column per time point -- a data
#'   frame or matrix, or a list of vectors of differing lengths, one per
#'   case, padded with \code{NA} to a common width internally (TDA's
#'   \code{seqdef} accepts \code{NA} for a case's unused trailing
#'   positions). States may be numeric, character or
#'   factor; the alphabet is taken over the whole frame.
#' @param method optimal matching method, TDA's \code{m=}.
#' @param indel indel cost, TDA's \code{icost=}.
#' @param subcost substitution costs, TDA's \code{scost=}: either a
#'   single number selecting one of TDA's built-in cost schemes, or
#'   a full \code{k x k} matrix of costs (\code{k} = number of distinct
#'   states in \code{sequences}), which gets written out as an
#'   \code{mdef()} block and referenced by name (TDA itself only
#'   accepts a matrix \emph{name} for \code{scost=}, never a literal
#'   list of numbers). Rows and columns follow the states in ascending
#'   order, TDA's internal ordering -- see \code{states} in the
#'   return value for exactly what that order is for a given call.
#' @param common_length require all sequences to share one length before
#'   aligning them (\code{seqm}'s \code{rr=1}); TDA's manual text
#'   for this is brief ("use common sequence length"), so this is
#'   documented no further than that.
#' @param max_restrict restrict the alignment search (\code{seqm}'s
#'   \code{max=}); TDA's manual text for this is brief ("alignment
#'   restriction"), so this is documented no further than that.
#' @param compare_with compare every sequence against specific sequences
#'   rather than every pair (\code{seqm}'s \code{cn=}, one or more
#'   case numbers); TDA's manual text for this is brief ("compare
#'   with specified sequences"), so this is documented no further than
#'   that.
#' @param preprocess one or both of \code{"skip_gaps"} (ignore internal
#'   gaps -- missing values -- during alignment) and
#'   \code{"skip_identical"} (collapse runs of identical states before
#'   aligning) -- \code{seqm}'s \code{sm=}.
#' @param n_random request an approximate random sample of this many
#'   sequence pairs rather than comparing every pair -- \code{seqm}'s
#'   own \code{r=}.
#' @param print add further columns to \code{attr(x, "pairs")} beyond
#'   the default index/length/distance: \code{"sequential"} replaces
#'   the single \code{distance} column with a running, cumulative-so-far
#'   cost at each time point (\code{dist_t1..dist_tn}, \code{n} = number
#'   of time points in \code{sequences}), also built into
#'   \code{attr(x, "dist_by_time")}, an \code{n x n x n_time} array --
#'   \code{[, , t]} is the full pairwise distance matrix using only the
#'   alignment cost accumulated up to time point \code{t}, the same
#'   symmetric, zero-diagonal shape \code{x} itself has. \code{distance}
#'   itself becomes the final time point, \code{dist_tn}, kept for
#'   compatibility with the default shape. \code{"lcs"} adds, alongside
#'   the usual single \code{distance}, the longest common subsequence's
#'   own length (\code{lcs_length}) and its state sequence, one column
#'   per time point (\code{lcs_t1..lcs_tn}, \code{-1}-padded past the
#'   LCS's length).
#'   \code{"sequential"} errors if \code{compare_with} is also given, its
#'   shape there not being one this wrapper parses; \code{"lcs"} with
#'   \code{compare_with} runs but leaves its extra columns unnamed
#'   for the same reason (one distance/LCS group per reference sequence,
#'   not the single group above) -- read \code{attr(x, "pairs")} (or
#'   \code{table} for the \code{compare_with} shape) directly for these.
#' @param dp_matrix capture TDA's step-by-step dump of the alignment
#'   algorithm's internal dynamic-programming matrix (\code{tst=2,3} +
#'   \code{df=}), one block per pair actually compared -- an
#'   illustration of how the distance was arrived at, not part of
#'   TDA's regular output at all. When \code{TRUE}, \code{dp_matrix} in
#'   the return value is a list keyed \code{"caseA_caseB"}, each holding
#'   \code{seqA}/\code{seqB} (the two sequences, states in order) and a
#'   \code{D} matrix (plus \code{E}/\code{F} for \code{method}s with
#'   affine costs) -- \code{D[i, j]} is the alignment cost using only
#'   the first \code{i} states of \code{seqA} and first \code{j} of
#'   \code{seqB}, so \code{D}'s bottom-right corner is that pair's final
#'   distance. Off by default: turning it on can produce a large file,
#'   one block per pair compared.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return Ordinarily a \code{dist}, with the pairwise table kept in its
#'   \code{pairs} attribute, and, with \code{print = "sequential"}, the
#'   per-time-point breakdown also in \code{dist_by_time} (see
#'   \code{print}). Either shape also carries \code{states}, the
#'   distinct states \code{sequences} actually has, in the ascending
#'   order \code{subcost} (if given as a matrix) must follow, and, with
#'   \code{dp_matrix = TRUE}, the per-pair alignment breakdown in
#'   \code{dp_matrix} (see \code{dp_matrix}). Also carries \code{run},
#'   the underlying \code{\link{tda_run}} result -- \code{run$dir} is
#'   where the generated TDA commands (\code{$run$commands}) and, with
#'   \code{dp_matrix = TRUE}, the raw, unparsed debug file TDA itself
#'   wrote (\code{dp.tst}) both sit, if what \code{dp_matrix} parsed out
#'   of it needs checking directly. With \code{compare_with} given,
#'   there is no full pairwise matrix to return -- \code{cn=} writes
#'   each sequence's distance to the reference set instead -- so
#'   the return value is instead an object carrying that as \code{table}
#'   (\code{case}, \code{length}, \code{distance}), and \code{states}
#'   the same way.
#' @family sequence analysis
#' @examples
#' s <- data.frame(t1 = c(1, 1, 2), t2 = c(1, 2, 2), t3 = c(2, 2, 3))
#' d <- tda_seqm(s)
#' hclust(d)
#' @export
tda_seqm <- function(sequences, method = 1, indel = NULL, subcost = NULL,
                     common_length = FALSE, max_restrict = NULL,
                     compare_with = NULL, preprocess = NULL,
                     n_random = NULL, print = NULL, dp_matrix = FALSE,
                     options = list(), dir = tempfile("tda"), ...) {
    if (is.list(sequences) && !is.data.frame(sequences)) {
        # A list of vectors is the natural R shape for sequences of
        # differing lengths -- the whole reason an edit-distance method
        # like this exists at all -- but as.data.frame() on such a list
        # fails outright ("differing number of rows"), so this used to
        # only ever work for sequences that already happened to share
        # one length. TDA's seqdef (m=1) expects one column per time
        # point with a fixed width, but does accept NA for a case's
        # unused trailing positions -- checked:
        # a real run reports "Sequences with zero length or internal
        # gaps: 0" and completes normally with NA-padded rows.
        len <- max(lengths(sequences))
        sequences <- as.data.frame(do.call(rbind, lapply(sequences,
            function(s) { length(s) <- len; s })))
    }
    d <- as.data.frame(sequences)
    if (!ncol(d))
        stop("no sequence columns")
    # Work out the state alphabet once, over the whole frame, before
    # converting anything: computing it inside the loop would see columns
    # already turned into integers and shift the mapping.
    chr <- vapply(d, function(z) is.character(z) || is.factor(z), logical(1))
    if (any(chr)) {
        lev <- sort(unique(unlist(lapply(d[chr], as.character))))
        for (j in which(chr))
            d[[j]] <- as.integer(factor(as.character(d[[j]]), levels = lev))
    }
    for (j in seq_along(d))
        d[[j]] <- as.integer(d[[j]])
    nm <- paste0("Y", seq_len(ncol(d)))
    names(d) <- nm
    # TDA's substitution-cost matrix (scost=NAME) is indexed by each
    # state's 0-based rank among the sorted distinct states it sees --
    # SCost[i*SeqNS+j] = MatVal[i*smc+j+1] in t_seqm.c, SeqNS being the
    # count of distinct states, not their raw values -- so this is the
    # order a `subcost` matrix's rows/columns need to follow, confirmed
    # directly against the source rather than assumed from scost='s name.
    # Any negative value, not just -1 specifically, is TDA's gap
    # marker here (seq_sget/get_sdata in t_seq.c: `s < 0` is "no state",
    # checked as a plain sign test, not a comparison to -1) -- so it's
    # excluded the same way NA already is, whether the gap arrived as
    # R's NA or as a literal negative value straight from the data
    # (seq.d1's -1 padding, say, if `sequences` wasn't converted to
    # NA first).
    states <- sort(unique(unlist(d, use.names = FALSE)))
    states <- states[!is.na(states) & states >= 0]

    opts <- c(list(m = method), .tda_extra(options))
    if (!is.null(indel))
        opts$icost <- paste(indel, collapse = ",")
    scost_mdef <- NULL
    if (!is.null(subcost)) {
        if (is.matrix(subcost) || is.data.frame(subcost)) {
            sm <- as.matrix(subcost)
            storage.mode(sm) <- "double"
            k <- length(states)
            if (nrow(sm) < k || ncol(sm) < k)
                stop("`subcost` must be at least ", k, " x ", k,
                     " -- the number of distinct states in `sequences`, ",
                     "in ascending order: ", paste(states, collapse = ", "))
            if (any(sm < 0, na.rm = TRUE))
                stop("`subcost` must not contain negative values")
            # mdef takes its values row by row, hence the transpose --
            # matches the existing mdef() emitter in smoothing.R.
            scost_mdef <- sprintf("mdef(SCOST,%d,%d) = %s;",
                                  nrow(sm), ncol(sm),
                                  paste(format(t(sm), digits = 15),
                                       collapse = ","))
            opts$scost <- "SCOST"
        } else
            opts$scost <- paste(subcost, collapse = ",")
    }
    if (isTRUE(common_length)) opts$rr <- 1
    if (!is.null(max_restrict)) opts$max <- max_restrict
    if (!is.null(compare_with)) opts$cn <- paste(compare_with, collapse = ",")
    # sm=/r=/s= were left unexposed earlier this session as too vague to
    # commit to a meaning from TDA's manual text alone
    # ("preprocessing", "random selection", "print ... distances, or
    # LCS"); read directly from t_seqm.c/t_parm.c instead of guessed
    # further, and confirmed against real output rather than the source
    # alone:
    # - sm=i1,i2,... sets ctx->PMSM[i] flags -- 1 skips internal gaps
    #   during alignment, 2 skips runs of identical states (confirmed by
    #   the source's comments directly above where each flag is
    #   read, not inferred from the option name).
    # - r=N requests a random sample of approximately N sequence pairs
    #   rather than every pair (ctx->PMR, printed back as "Random
    #   selection of approximately N sequence pairs" when it takes
    #   effect).
    # - s=1/s=2 (ctx->PMS) both add columns to the pairwise
    #   output file beyond the default (index/length/distance) --
    #   checked from the name: s=1 ("sequential")
    #   *replaces* the single DIST column with one running, cumulative-
    #   so-far cost per time point (DIST1..DISTn), not an addition
    #   alongside it -- confirmed against a real run's .tda
    #   description file. s=2 ("lcs") does add alongside it: the longest
    #   common subsequence's length and its state sequence (-1-padded).
    if (!is.null(preprocess))
        opts$sm <- paste(match(preprocess, c("skip_gaps", "skip_identical")),
                         collapse = ",")
    if (!is.null(n_random)) opts$r <- n_random
    if (!is.null(print)) {
        print <- match.arg(print, c("sequential", "lcs"))
        opts$s <- match(print, c("sequential", "lcs"))
    }
    if (isTRUE(dp_matrix)) {
        # tst=2,3 + df= is TDA's step-by-step dump of the alignment
        # algorithm's internal dynamic-programming matrix, one D matrix
        # (plus E and F, for method=2's affine costs) per pair actually
        # compared -- checked against t_seqm.c (seq_pps() for
        # tst=2's "Sequence A"/"Sequence B" lines, seq_dmat() for tst=3's
        # matrix itself) and against a real run's output, not
        # assumed from TDA's brief manual text for either option. This
        # is TDA's *internal* working matrix, not part of its regular
        # output at all -- turning it on can produce a large file, one
        # block per pair, so it's opt-in.
        opts$tst <- "2,3"
        opts$df <- "dp.tst"
    }

    res <- tda_run(c(tda_nvar(d), scost_mdef,
                     sprintf("seqdef = %s;", paste(nm, collapse = ",")),
                     # dtda= is TDA's description of the file it
                     # writes -- an nvar block naming every column, its
                     # width and its meaning. The manual prints it beside
                     # the data (6.7.2.3, Box 6), and it is the record of
                     # what out.d's columns are, so it is asked for and
                     # kept rather than left unwritten.
                     do.call(tda_block, c(list(name = "seqm"),
                                          list(dtda = "out.tda"), opts,
                                          list(rhs = "out.d")))),
                   data = d, dir = dir, ...)
    dpm <- if (isTRUE(dp_matrix))
        .seq_parse_dp(file.path(res$dir, "dp.tst"), res)

    # seqm.pairs carries every value the out.d row holds, so the file is
    # only read when the producer does not cover this layout
    epairs0 <- if (.use_exports()) res$exports[["seqm.pairs"]]
    tab <- if (is.matrix(epairs0)) as.data.frame(epairs0)
           else tda_file(res, "out.d")
    if (is.null(tab))
        stop("seqm produced no distances")
    seqfl <- identical(print, "sequential")
    if (!is.null(compare_with)) {
        # cn= writes a different shape: one row per sequence
        # against the reference set (case, length, distance), not the
        # all-pairs (i, j, len.i, len.j, distance) shape the default
        # comparison uses -- checked against real output, not
        # assumed to share it: a first attempt at cn= support kept the
        # default parsing and errored outright ("produced no distances")
        # on cn='s own, narrower column count.
        if (seqfl)
            # print="sequential" replaces DIST with one running-cost
            # column per time point *per reference sequence* here (t_seqm.c's
            # DIST%d_%d block), a layout this wrapper doesn't parse -- rather
            # than mislabel column 3 as "distance" the way the plain cn=
            # case below would, error instead of guessing at the shape.
            stop("`print` is not supported together with `compare_with` -- ",
                 "seqm's output widens per reference sequence in a ",
                 "shape this wrapper doesn't parse; use `tda_run` directly ",
                 "if you need it")
        names(tab)[1:3] <- c("case", "length", "distance")
        return(structure(list(call = match.call(), run = res, n = nrow(d),
                              table = tab, states = states, dp_matrix = dpm),
                         class = c("tda_seqm", "tda_table")))
    }
    if (ncol(tab) < 5L)
        stop("seqm produced no distances")
    if (seqfl) {
        # print="sequential" (s=1) doesn't add a column alongside the
        # usual single DIST -- it replaces it outright with one running,
        # cumulative-so-far cost per time point (t_seqm.c's dtda
        # writer: "DIST1 ... # distance, time 1", "DIST2 ... time 2", and
        # so on, one per time point in the sequence) -- confirmed against
        # a real run's .tda description file from the
        # option's name. The final total (what `distance` means
        # everywhere else in this package) is therefore the *last* of
        # these, not the column in the same position the default layout
        # would put "distance" at.
        nt <- ncol(d)
        if (ncol(tab) < 4L + nt)
            stop("seqm's output has fewer per-time-point distance ",
                 "columns than `sequences` has time points -- can't line ",
                 "up DIST1..DIST", nt)
        names(tab)[1:(4L + nt)] <- c("i", "j", "len.i", "len.j",
                                     paste0("dist_t", seq_len(nt)))
        tab$distance <- tab[[4L + nt]]
    } else {
        names(tab)[1:5] <- c("i", "j", "len.i", "len.j", "distance")
        if (identical(print, "lcs")) {
            # s=2 adds alongside the usual single DIST, rather than
            # replacing it (t_seqm.c's dtda writer: "LCSL ... #
            # length of LCS" then "LCS1 ... # LCS t=1", "LCS2 ...", one
            # per time point) -- confirmed against the source directly,
            # same as the sequential case above.
            nt <- ncol(d)
            if (ncol(tab) >= 6L + nt)
                names(tab)[6:(6L + nt)] <- c("lcs_length",
                                             paste0("lcs_t", seq_len(nt)))
        }
    }

    # seqm.pairs carries every value the out.d row holds, in the order
    # written -- the four identifiers, the s= distance columns, then
    # v='s -- so the whole frame is switched at once rather than
    # only the columns the default layout happens to have.
    desc <- file.path(res$dir, "out.tda")
    if (file.exists(desc))
        attr(tab, "description") <- readLines(desc, warn = FALSE)
    if (!is.matrix(epairs0))
        tab <- .overlay_num(tab, res$exports[["seqm.pairs"]])
    if (seqfl)
        tab$distance <- tab[[4L + nt]]
    n <- nrow(d)
    m <- matrix(0, n, n)
    m[cbind(tab$i, tab$j)] <- tab$distance
    m[cbind(tab$j, tab$i)] <- tab$distance
    byt <- NULL
    if (seqfl) {
        # dist_t1..dist_tn sit right next to each other in seqm's
        # output, one running-cost matrix's worth of pairs per time
        # point -- built into an n x n x nt array here rather than left
        # as flat columns in `pairs`, so a single time slice is one
        # `[, , t]` away and already has the usual symmetric,
        # zero-diagonal dist-matrix shape the default output has.
        byt <- array(0, dim = c(n, n, nt),
                     dimnames = list(NULL, NULL, paste0("t", seq_len(nt))))
        for (t in seq_len(nt)) {
            byt[, , t][cbind(tab$i, tab$j)] <- tab[[paste0("dist_t", t)]]
            byt[, , t][cbind(tab$j, tab$i)] <- tab[[paste0("dist_t", t)]]
        }
    }
    structure(stats::as.dist(m), call = match.call(),
              pairs = tab, dist_by_time = byt, states = states,
              dp_matrix = dpm, run = res, class = c("dist"))
}

#' Sequence pattern matching
#'
#' Counts, per sequence, how many times each of up to 20 patterns occurs --
#' TDA's \code{seqpm}. A pattern is a vector of states and wildcards, one
#' element per position: \code{c(3, 3)} looks for state 3 immediately
#' followed by state 3; \code{"?"} matches any one state, \code{"*"} any
#' run of states (including none), \code{"+"} any repeat of the state
#' before it, and \code{"-"} any run of \emph{identical} states. Several
#' patterns are given as a list, \code{list(c(3, 3), c(3, "*", 3))}.
#'
#' @param sequences a data frame or matrix of one column per time point,
#'   in order -- the same shape \code{\link{tda_seqm}} takes -- optionally
#'   with further, non-sequence columns (an id, say) alongside them; use
#'   \code{variables}/\code{id} to say which is which.
#' @param patterns a list of patterns, each a vector of states (numbers)
#'   and wildcards (\code{"?"}, \code{"*"}, \code{"+"}, \code{"-"}, as
#'   strings), as described above. At most 20.
#' @param variables which columns of \code{sequences} are the actual
#'   time-point columns, in order -- TDA's \code{seqdef=Y0,,Y7}, which
#'   likewise names the sequence out of a data matrix holding other
#'   variables too. A character vector of names in \code{sequences};
#'   defaults to every column (so an already-subsetted \code{sequences}
#'   still works unchanged).
#' @param id further columns of \code{sequences} to carry through to the
#'   output as-is, one value per case -- \code{seqpm}'s \code{v=}. A
#'   character vector of names in \code{sequences}.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object carrying a \code{table} with one row per case: its
#'   sequence id, its length, one match-count column per pattern, and any
#'   \code{id} columns requested.
#' @family sequence analysis
#' @examples
#' # six short sequences over three states; does "3,3" (state 3 twice in a
#' # row) occur, and does the sequence ever return to 1 after leaving it?
#' s <- data.frame(Y1 = c(1, 1, 3, 3, 3, 1), Y2 = c(1, 3, 3, 3, 1, 2),
#'                 Y3 = c(1, 3, 3, 1, 2, 2), Y4 = c(2, 3, 1, 2, 2, 2),
#'                 Y5 = c(2, 1, 2, 2, 2, 3))
#' r <- tda_seqpm(s, list(c(3, 3), c(1, "*", 1)))
#' r$table
#'
#' # an id column alongside the sequence, positively selected
#' s2 <- cbind(person = 101:106, s)
#' r2 <- tda_seqpm(s2, list(c(3, 3)), variables = names(s), id = "person")
#' r2$table
#' @export
tda_seqpm <- function(sequences, patterns, variables = NULL, id = NULL,
                      options = list(), dir = tempfile("tda"), ...) {
    picked <- .seq_select(sequences, variables, id, allow_id = TRUE,
                          cmd = "seqpm")

    d <- as.data.frame(picked$sequences)
    if (!ncol(d))
        stop("no sequence columns")
    chr <- vapply(d, function(z) is.character(z) || is.factor(z), logical(1))
    if (any(chr)) {
        lev <- sort(unique(unlist(lapply(d[chr], as.character))))
        for (j in which(chr))
            d[[j]] <- as.integer(factor(as.character(d[[j]]), levels = lev))
    }
    for (j in seq_along(d))
        d[[j]] <- as.integer(d[[j]])
    nm <- paste0("Y", seq_len(ncol(d)))
    names(d) <- nm

    if (!is.list(patterns))
        patterns <- list(patterns)
    if (!length(patterns))
        stop("give at least one pattern")
    if (length(patterns) > 20L)
        stop("seqpm allows at most 20 patterns")
    ps <- paste(vapply(patterns, function(p)
        sprintf("[%s]", paste(p, collapse = ",")), character(1)),
        collapse = ",")

    vnm <- NULL
    if (!is.null(picked$id)) {
        vnm <- .tda_names(names(picked$id))
        d <- cbind(d, stats::setNames(picked$id, vnm))
    }
    opts <- c(list(ps = ps), .tda_extra(options))
    if (!is.null(vnm))
        opts$v <- paste(vnm, collapse = ",")
    res <- tda_run(c(tda_nvar(d),
                     sprintf("seqdef = %s;", paste(nm, collapse = ",")),
                     do.call(tda_block, c(list(name = "seqpm"), opts,
                                          list(rhs = "out.d")))),
                   data = d, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run seqpm: ", err[1L], call. = FALSE)
    # seqpm's output is always case, length, then one match count per
    # pattern in order -- confirmed against its dtda description, not
    # assumed -- so the columns can be named directly without reading it.
    # id= appends after that, in the order given, same as v= does for the
    # .seq_desc-based commands. The first column is named "case", not
    # "id", so it never collides with an id= column of that name.
    etab0 <- if (.use_exports()) res$exports[["seqpm.table"]]
    tab <- if (is.matrix(etab0)) .export_frame(etab0)
           else tda_file(res, "out.d")
    if (is.null(tab))
        stop("seqpm produced no output")
    npat <- length(patterns)
    names(tab)[seq_len(2L + npat)] <- c("case", "length",
                                        paste0("matches", seq_len(npat)))
    if (!is.null(vnm) && ncol(tab) == 2L + npat + length(vnm))
        names(tab)[(3L + npat):ncol(tab)] <- vnm
    tab <- .overlay_num(tab, res$exports[["seqpm.table"]])
    structure(list(call = match.call(), run = res, n = nrow(d),
                   table = tab), class = c("tda_seqpm", "tda_table"))
}


#' Maximum likelihood with a user-defined likelihood
#'
#' \code{fml} maximises a log-likelihood you write yourself, which is what TDA
#' offers that a fixed catalogue of models cannot. The likelihood is written
#' in an expression language, evaluated inside TDA rather than as an R
#' function -- but as of this version, that expression can be written in
#' ordinary R syntax and translated automatically, rather than TDA's.
#'
#' @section Writing the likelihood in R syntax:
#' Give \code{definitions} as an unevaluated \code{\{ \}} block of plain R
#' assignments instead of a character vector, and it is translated to TDA's
#' own expression syntax for you -- never actually run as R code (\code{DES},
#' \code{rate} and the rest only ever exist inside TDA, so it cannot be).
#' \code{ifelse(cond, yes, no)} becomes TDA's \code{if()} (TDA's
#' \code{if} takes a 2- or 3-argument
#' \code{if(cond, then)}/\code{if(cond, then, else)} shape matching
#' \code{ifelse} exactly, not R's control-flow \code{if}), and
#' R's comparison operators (\code{<}, \code{>}, \code{<=}, \code{>=},
#' \code{==}, \code{!=}) become TDA's comparison functions (\code{lt},
#' \code{gt}, \code{le}, \code{ge}, \code{eq}, \code{ne}). Arithmetic
#' (\code{+ - * / ^}), \code{&}/\code{|}, and ordinary function calls
#' (\code{exp}, \code{log}, \code{negbin}, ...) already mean the same thing
#' in both and pass through unchanged.
#'
#' \strong{These six TDA functions are not the strict comparisons their R
#' equivalents are.} Every one of them -- \code{lt}, \code{le}, \code{gt},
#' \code{ge}, \code{eq}, \code{ne} alike -- carries a fixed tolerance,
#' \code{EPSI2 = 1000 * DBL_EPSILON} (about \code{2.22e-13}): \code{le},
#' say, computes \code{x <= y + EPSI2}, not \code{x <= y}. This can
#' matter on real data: a case in \code{examples/exam/rrdat.1}
#' has \code{PRESN/PRES - 1} equal to \code{0.2} exactly in exact arithmetic,
#' and the two land on different sides of \code{0.2} once represented as
#' IEEE 754 doubles -- \code{ge(x, 0.2)} inside TDA still calls it true, R's
#' own \code{x >= 0.2} on the identical computation does not. Applied
#' consistently across all six, at least, but nowhere in TDA's
#' documentation: not in the manual, not in \code{tda_help()}'s
#' operator table, which describes \code{le} as plainly as \dQuote{true if
#' less than, or equal} with no footnote at all. A hand-written R
#' replication of a TDA calculation involving a comparison near a round
#' number (\code{==}, or a boundary like \code{0.2} here) can therefore
#' disagree with TDA's output at that boundary, for a
#' reason with nothing to do with which language computed it -- this is a
#' property of TDA's expression evaluator, present whether the
#' comparison was written by hand in TDA's syntax or reached it through
#' this translation.
#'
#' \preformatted{
#' tda_fml({
#'     rate = exp(a0 + COHO2 * a1 + COHO3 * a2 + W * a3)
#'     l1 = ifelse(DES, log(rate), 0)
#'     fn = l1 - rate * DUR
#' }, data = d, start = c(-4, 0, 0, 0))
#' }
#'
#' The older, direct form -- a character vector of TDA's assignment
#' syntax -- still works exactly as before, for anything the translator
#' does not (yet) cover, or for pasting a \code{.cf} file's text in
#' directly.
#'
#' \preformatted{
#' tda_fml(c("xb = b0 + x * b1",
#'           "fn = -0.5 * (y - xb)^2"), data = d)
#' }
#'
#' estimates \code{b0} and \code{b1} by least squares written as a likelihood.
#' TDA's functions are available, so \code{fn = negbin(alpha, gamma, y)}
#' and the rest work.
#'
#' @param definitions either an unevaluated \code{\{ \}} block of R
#'   assignments (see the section on writing the likelihood) or a
#'   character vector of TDA's
#'   assignments, the last defining \code{fn}, the contribution of one
#'   case to the log-likelihood. Names that appear on the right but are
#'   neither variables in \code{data} nor defined earlier are the
#'   parameters to be estimated.
#' @param data a data frame; its columns are the variables the definitions can
#'   refer to.
#' @param start optional named vector of starting values, e.g.
#'   \code{c(b0 = 0, b1 = 1)}.
#' @param constraints optional linear constraints on the parameters --
#'   \code{fml}'s \code{con=}, the same machinery
#'   \code{\link{tda_qreg}}'s \code{constraints} uses. \strong{Refers to parameters by
#'   their position, \code{b1}/\code{b2}/..., 1-based, in the order they
#'   are first introduced across \code{definitions} -- never by whatever
#'   name was actually given them.} A parameter named \code{b0} in
#'   \code{definitions} is still \code{b1} here if it is the first one
#'   introduced; TDA's error for a literal \code{"b0 = ..."}
#'   constraint is \dQuote{Error in parameter index}, since \code{b0}
#'   is not a valid 1-based position. Run once without
#'   \code{constraints} and check
#'   \code{tda_estimates()}'s \code{Idx} column for which position is
#'   which parameter before writing one. Must be strictly fewer
#'   constraints than parameters (TDA's rule; its error reads
#'   \dQuote{number of constraints should be less than number of
#'   parameters}).
#' @param residuals ask TDA to also compute, per case, its
#'   contribution to the log-likelihood at the converged parameters --
#'   \code{fml}'s \code{pres=}. \strong{Not a classical observed-minus-fitted
#'   residual} -- \code{fml} has no built-in notion of a fitted value to
#'   subtract from, unlike \code{\link{tda_lsreg}}'s
#'   \code{residuals}. For
#'   \code{fn = -0.5 * (y - xb)^2}, say, this is that expression's
#'   value per case, not \code{y - xb} itself.
#' @param residual_vars with \code{residuals = TRUE}, extra columns
#'   from \code{data} to write alongside \code{fn}'s value, one
#'   column per case, in the order given -- \code{fml}'s
#'   \code{v=}. \code{$residuals} then comes back as a data
#'   frame (\code{fn} plus each named column) instead of a plain
#'   vector. Has no effect without \code{residuals = TRUE}.
#' @param protocol ask TDA to also write its iteration-by-iteration
#'   diagnostic log -- \code{fml}'s \code{prot=}: the iteration
#'   table (function value, gradient norm, parameter change per step),
#'   the parameter vector, and the covariance matrix at convergence.
#'   Not one consistent shape across every algorithm and situation, so
#'   returned as \code{$protocol}, the file's text lines, rather
#'   than force-parsed into a table that would not fit every case.
#' @param control convergence settings from \code{\link{tda_control}}.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object of class \code{tda_fit}. With \code{residuals =
#'   TRUE}, also carries \code{$residuals} (see \code{residuals} above).
#' @family sequence analysis
#' @examples
#' set.seed(1)
#' d <- data.frame(x = 1:20, y = 2 + 0.5 * (1:20) + rnorm(20))
#' tda_fml(c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"), d)
#'
#' # the same model, written in R syntax instead
#' tda_fml({
#'     xb = b0 + x * b1
#'     fn = -0.5 * (y - xb)^2
#' }, d)
#'
#' # a compound Poisson (negative binomial) count regression -- TDA has no
#' # separate command for this (unlike its plain Poisson, which is a glm()
#' # family), but negbin(alpha, gamma, k) is a built-in log-likelihood
#' # function in its expression language, reachable through fml like any
#' # other: gamma is the mean, alpha = 1/sigma the dispersion
#' set.seed(7)
#' d2 <- data.frame(x = rnorm(200))
#' d2$ndi <- rnbinom(200, mu = exp(0.5 + 0.8 * d2$x), size = 2)
#' coef(tda_fml({
#'     xb = b0 + x * b1
#'     gamma = exp(xb)
#'     sigma = exp(sig)
#'     alpha = 1 / sigma
#'     fn = negbin(alpha, gamma, ndi)
#' }, d2, start = list(b0 = 0, b1 = 0, sig = 0)))
#' @export
tda_fml <- function(definitions, data, start = NULL, control = NULL,
                    constraints = NULL, residuals = FALSE,
                    residual_vars = NULL, protocol = FALSE,
                    options = list(), dir = tempfile("tda"), ...) {
    definitions <- .fml_check_definitions(definitions, substitute(definitions),
                                          "the log-likelihood contribution of one case")
    d <- as.data.frame(data)
    lab <- names(d)
    names(d) <- .tda_names(lab)
    # The definitions are written against the data frame's column names,
    # so they have to be rewritten to the names TDA was given.
    for (i in seq_along(lab))
        if (lab[i] != names(d)[i])
            definitions <- gsub(sprintf("\\b%s\\b", lab[i]), names(d)[i],
                                definitions)
    rv_tda <- if (!is.null(residual_vars))
        names(d)[match(residual_vars, lab)]

    opts <- .fml_opts(start, constraints, control, options)
    opts <- .fml_residuals_opts(opts, residuals, rv_tda)
    opts <- .fml_protocol_opts(opts, protocol)
    cmd <- do.call(tda_block, c(list(name = "fml"), opts,
                                list(rhs = paste(definitions, collapse = ",\n    "))))
    res <- tda_run(c(tda_nvar(d), cmd), data = d, dir = dir, ...)
    .tda_result(res, list(n = nrow(d), xlab = character(), xname = character()),
                match.call(), "tda_fml",
                list(residuals = .fml_residuals_read(res, residuals, residual_vars),
                     protocol = .fml_protocol_read(res, protocol)))
}

# Translates a { } block of plain R assignments into TDA's expression
# syntax, for tda_fml()'s R-syntax form -- never evaluated as R (rate, DES
# and the rest only ever exist inside TDA), only ever walked as an
# unevaluated call tree and rewritten. ifelse() becomes TDA's if()
# (t_eval.c: 2- or 3-argument, matching ifelse's shape exactly) and R's
# comparison operators become TDA's comparison functions (t_eval.c:
# lt/gt/le/ge/eq/ne) -- both checked against the source, not
# assumed from a naming convention. Arithmetic, &/|, and ordinary function
# calls already mean the same thing in both and are left untouched.
#
# ifelse() is rewritten to a placeholder symbol, .TDA_IF, not literally to
# `if` -- checked: R's deparse() special-cases
# any call whose head is the symbol `if`, printing it back out as
# if (cond) then else else_expr control-flow syntax regardless of how the
# call tree was actually built, which is not TDA's if(cond, then, else)
# function-call syntax at all. The placeholder survives deparse() as an
# ordinary function call and is fixed up textually afterward, once it is
# safely just a string.
.fml_ops <- c("<" = "lt", ">" = "gt", "<=" = "le", ">=" = "ge",
             "==" = "eq", "!=" = "ne", "ifelse" = ".TDA_IF")

.fml_rewrite <- function(e) {
    if (!is.call(e))
        return(e)
    parts <- as.list(e)
    head <- parts[[1L]]
    args <- lapply(parts[-1L], .fml_rewrite)
    if (is.symbol(head)) {
        op <- as.character(head)
        if (op %in% names(.fml_ops))
            head <- as.symbol(.fml_ops[[op]])
    }
    as.call(c(head, args))
}

.fml_translate <- function(block) {
    stmts <- as.list(block)[-1L]
    if (!length(stmts))
        stop("`definitions` block is empty")
    vapply(stmts, function(s) {
        if (!(is.call(s) && (identical(s[[1L]], quote(`<-`)) ||
                             identical(s[[1L]], quote(`=`)))))
            stop("each line in the `definitions` block must be a plain ",
                 "assignment, e.g. `rate = exp(...)` -- got: ",
                 paste(deparse(s), collapse = " "))
        lhs <- deparse(s[[2L]])
        rhs <- paste(deparse(.fml_rewrite(s[[3L]]), width.cutoff = 500L),
                    collapse = " ")
        rhs <- gsub(".TDA_IF(", "if(", rhs, fixed = TRUE)
        paste(lhs, "=", rhs)
    }, character(1L))
}

# Shared plumbing for the whole fml/freg/frml family: all three go
# through TDA's f_min() (typ=0 ML for fml, typ=2 nonlinear
# regression for freg, typ=3 for frml -- t_gmin.c's comment),
# and all three accept the same xp=/con=/pres= options --
# checked for freg and frml too from fml alone,
# by running each against the real binary: con=bN=value correctly
# fixes the Nth parameter (by position, never by whatever name it was
# given -- the same trap documented on tda_fml's constraints), and
# pres= writes each case's (or, for frml, each episode's) own fn value
# at the converged parameters, matching an independent evaluation in R
# exactly, for both.
#
# Definitions in { } block form (plain R, ifelse()/comparisons
# translated automatically -- see .fml_translate() above) or as a
# character vector of TDA's syntax are both accepted by all three,
# through this one shared entry point rather than three separate
# copies of the same substitute()/is.call() dance.
.fml_check_definitions <- function(definitions, sub_def, fn_desc) {
    if (is.call(sub_def) && identical(sub_def[[1L]], quote(`{`)))
        definitions <- .fml_translate(sub_def)
    if (!is.character(definitions) || !length(definitions))
        stop("`definitions` must be a character vector of TDA assignments")
    definitions <- sub(",\\s*$", "", trimws(definitions))
    if (!any(grepl("^fn\\s*=", definitions)))
        stop("the last definition must define fn, ", fn_desc)
    definitions
}

# start=/constraints=/residuals= all become plain f_min() options
# (xp=/con=/pres=) the same way for fml, freg, and frml alike --
# built once here instead of three times. tfmt=24.16 is included by
# default (mfmt=24.16 always is); freg's default control -- Newton
# overflowing exp() on nonlinear-regression-shaped functions -- is
# handled by its caller before this, not here, since it depends on
# `control` already being non-NULL by the time this runs.
#
# ccov= (covariance matrix type) is deliberately not exposed anywhere
# in tda_control() or here: dead for this whole
# family, not just undocumented. t_gmin.c parses it into ctx->CCTyp,
# but the only code that would ever read that value back is commented
# out (`/* if (CCTyp > 2) CCTyp = 2; */`); the `cov` flag actually used
# is set unconditionally by separate logic a few lines later,
# regardless of whatever ccov= was given.
.fml_opts <- function(start, constraints, control, options) {
    opts <- c(list(mfmt = "24.16", tfmt = "24.16"),
              .control_opts(control), .tda_extra(options))
    if (!is.null(start))
        # xp= takes the values in parameter order, not name=value pairs.
        # Names are kept for the caller's benefit and checked against the
        # order the definitions introduce the parameters in.
        opts$xp <- paste(format(unname(start), trim = TRUE), collapse = ",")
    if (!is.null(constraints))
        # con=: TDA's bN numbering refers to a parameter's
        # *position*, 1-based, in the order the definitions first
        # introduce it -- never whatever name it was actually given.
        # Must be strictly fewer constraints than parameters (TDA's
        # own rule). See tda_fml's documentation for the full
        # explanation and a worked example of the naming trap.
        opts <- c(opts, stats::setNames(as.list(constraints),
                                        rep("con", length(constraints))))
    opts
}

# residuals= needs one more option (pres=, plus fmt= for its own
# precision, the same gap already found and fixed for qreg's df=) on
# top of what .fml_opts() already built, and a matching read-back step
# once the run is done. residual_vars adds v= (extra columns alongside
# fn's value in the same file) -- checked to genuinely
# work: a real run with v=Dose,Weight wrote exactly those
# two extra columns. Only meaningful together with residuals = TRUE;
# silently has no effect otherwise, the same as passing an option TDA
# itself would just ignore for a command that never writes the file at
# all.
.fml_residuals_opts <- function(opts, residuals, residual_vars = NULL) {
    if (isTRUE(residuals)) {
        opts$pres <- "res.out"
        opts$fmt <- "24.16"
        if (!is.null(residual_vars))
            opts$v <- paste(residual_vars, collapse = ",")
    }
    opts
}

.fml_residuals_read <- function(res, residuals, residual_vars = NULL) {
    if (!isTRUE(residuals))
        return(NULL)
    # Not a classical observed-minus-fitted residual -- this family has
    # no built-in notion of a fitted value to subtract from -- each
    # case's (or episode's) own fn value at the converged parameters
    # (gm_pres()/get_fival() in t_gmin.c), checked by
    # evaluating fn independently in R and matching it exactly, for
    # fml, freg, and frml alike. Two columns in the file without
    # residual_vars (case index, that value -- the second is what
    # callers want, returned as a plain vector); with residual_vars,
    # the extra columns follow in the order given, and the whole thing
    # comes back as a data frame instead, fn plus each named column.
    rt <- tryCatch(tda_file(res, "res.out"), error = function(e) NULL)
    if (is.null(rt))
        return(NULL)
    if (!is.null(residual_vars) &&
        ncol(rt) >= 2L + length(residual_vars)) {
        names(rt)[seq_len(2L + length(residual_vars))] <-
            c("case", "fn", residual_vars)
        return(rt[c("fn", residual_vars)])
    }
    if (ncol(rt) >= 2L) rt[[2L]] else rt
}

# protocol=: TDA's iteration-by-iteration diagnostic log (prot=) --
# useful (the iteration table, gradient norms, the parameter
# vector and covariance matrix at each stage), checked by
# running it, but not a table with one consistent shape across every
# algorithm and situation, so returned as the plain text lines rather
# than force-parsed into a data frame that would not fit every case.
.fml_protocol_opts <- function(opts, protocol) {
    if (isTRUE(protocol))
        opts$prot <- "prot.out"
    opts
}

.fml_protocol_read <- function(res, protocol) {
    if (!isTRUE(protocol))
        return(NULL)
    pr <- .file_lines(res, "prot.out")
    if (!is.null(pr))
        return(pr)
    p <- file.path(res$dir, "prot.out")
    if (!file.exists(p))
        return(NULL)
    readLines(p, warn = FALSE)
}


# Neither ineq nor spl writes a column header into its output file -- the
# names are only in the console table -- so they are supplied here.
# Rows of a table TDA printed but did not write to a file: every line with
# exactly `n` numeric fields and nothing else.
.numeric_rows <- function(txt, n, nm) {
    rows <- list()
    for (l in txt) {
        tok <- strsplit(trimws(l), "\\s+")[[1L]]
        if (length(tok) != n)
            next
        v <- suppressWarnings(as.numeric(tok))
        if (!anyNA(v))
            rows[[length(rows) + 1L]] <- v
    }
    if (!length(rows))
        return(NULL)
    d <- as.data.frame(do.call(rbind, rows))
    names(d) <- nm
    d
}

.name_cols <- function(tab, nm, labels = NULL) {
    if (is.null(tab))
        return(NULL)
    if (ncol(tab) <= length(nm))
        names(tab) <- nm[seq_len(ncol(tab))]
    if (!is.null(labels) && "index" %in% names(tab) &&
        all(tab$index %in% seq_along(labels)))
        tab$index <- labels[tab$index]
    tab
}


#' Regression models for events
#'
#' Turns sequence data into the discrete-time expanded form a logistic or
#' probit model needs: one row per individual per time point at risk, a
#' binary indicator of whether the target event (a transition from one
#' state to another) occurred there, dummy variables for each time period
#' (so period effects can be estimated freely, the way \code{factor(t)}
#' would in a regular model), and any covariates carried along. Fit the
#' result with \code{\link{tda_qreg}} (\code{model = "logit"} or
#' \code{"probit"}) or \code{\link{tda_glm}(family = binomial)} on the
#' \code{event} column against the period dummies and covariates --
#' TDA's manual (chapter 6.18) calls this approach out explicitly as
#' the discrete-time counterpart to \code{\link{tda_rate}}.
#'
#' \code{tda_seqev} and \code{tda_seqevd} are the exploratory step before
#' committing to a target event: \code{tda_seqev} counts how often each
#' state-to-state transition actually occurs in the data, and
#' \code{tda_seqevd} tabulates event counts over time.
#'
#' \strong{A note on names}: the manual documents this command as
#' \code{seqnd}, but TDA 6.4 itself calls it \code{seqmd};
#' \code{tda_seqmd()} uses the name that works.
#'
#' This wrapper always defines exactly one sequence data structure per
#' call (a fresh \code{seqdef}, every time), so TDA's \code{sn=} --
#' selecting among several already-defined structures -- has nothing to
#' select among here: \code{sn = 1} is always correct, and there is no
#' way to reach any other value through this interface. Not exposed as
#' an argument for that reason, rather than included as a parameter
#' that could never meaningfully differ from its default.
#'
#' \code{tda_seqev} ignores TDA's \code{dtda=} (no description file
#' is written even when asked for one); its output shape
#' (\code{from}/\code{to}/\code{count}) needs no further description to
#' read reliably. \code{tda_seqmd}'s \code{dtda=} is used
#' internally to read TDA's real column structure back, since
#' \code{event_covariates}, in particular, can change how many columns
#' come back and where.
#'
#' @param sequences a data frame or matrix, one row per case, one column
#'   per time point, in order -- the same shape \code{\link{tda_seqm}}
#'   takes -- or, given \code{data} too, a character vector naming
#'   those columns within it directly, rather than requiring a
#'   separately pre-subsetted data frame: real sequence data typically
#'   carries an ID column and covariates alongside the states in the
#'   same data frame (TDA's \code{examples/exam/seq.d4}, e.g.), not
#'   only the states on their own. For \code{tda_seqev}/\code{tda_seqevd},
#'   also a \emph{list} of several such specs, to define several
#'   independent sequence data structures at once -- TDA's
#'   \code{seqdef(sn=1)}, \code{seqdef(sn=2)}, ... -- the pattern TDA's
#'   own manual uses to build a second sequence structure from
#'   \code{seq.d4}'s \code{S0..S5} columns, entirely separate from the
#'   first (\code{Y0..Y5}), not a covariate of it. Use \code{sn=} to
#'   say which one a given call analyzes.
#' @param sn for \code{tda_seqev}/\code{tda_seqevd}, which sequence data
#'   structure to use when \code{sequences} defines more than one (see
#'   above) -- TDA's \code{sn=}, 1 by default, the same as when only
#'   one structure was ever defined.
#' @param data optional data frame \code{sequences} (and, for
#'   \code{tda_seqmd}, \code{covariates}) name columns within, instead
#'   of \code{sequences}/\code{covariates} being the actual data
#'   already.
#' @param event for \code{tda_seqmd}, the target transition, as
#'   \code{c(from, to)} -- the two states such that going from the first
#'   to the second at time \code{t} counts as the event.
#' @param covariates optional further time-independent covariates,
#'   carried into the output unchanged -- a data frame or matrix, one
#'   row per case in the same order as \code{sequences}, or (with
#'   \code{data} given) a character vector naming columns within it.
#' @param event_covariates optional time-varying covariates -- TDA's
#'   \code{xe=}, different from \code{covariates}/\code{v=}:
#'   a separate value for each time point, rather than one fixed value
#'   per case. A named list, needs \code{data}: each element the
#'   per-time-point column names for one such covariate, in the same
#'   order as \code{sequences} itself, named by what the result's
#'   column should be called -- \code{list(price = paste0("Price", 0:5))},
#'   say, for a covariate stored one column per time point the same way
#'   the states themselves are. \strong{Exactly how TDA derives the
#'   value it puts in the output is not documented} -- it is not simply
#'   the value at the event's time point; check against your data
#'   before relying on a particular reading of this column.
#' @param select optional expression (\code{sel=}) restricting which
#'   sequences are used -- written with R's comparison operators
#'   (\code{Y0 == 1}, or \code{Y0 = 1}) or TDA's function form
#'   (\code{eq(Y0,1)}) interchangeably; \code{&}/\code{|} combine several
#'   conditions the same way in both.
#' @param tp optional restriction of the time axis, as a vector of time
#'   points or a single TDA range expression such as \code{"0(1)10"}.
#' @param summary instead of the full per-case, per-period data, just a
#'   \code{time}/\code{riskset}/\code{events} table -- TDA's
#'   \code{seqmd(ev=...)} called with no output file at all, a separate
#'   mode in which TDA prints this table to the console
#'   instead of writing any records. \code{covariates}/
#'   \code{event_covariates} are not accepted together with this, since
#'   TDA's console summary carries no covariate columns to show.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return \code{tda_seqmd} returns a data frame: \code{id}, \code{time},
#'   \code{event} (the binary indicator), one \code{period} dummy column
#'   per distinct time point, and any \code{covariates} -- or, with
#'   \code{summary = TRUE}, \code{time}, \code{riskset}, \code{events}
#'   instead. \code{tda_seqev}
#'   returns a data frame, one row per transition that occurs
#'   in the data: \code{from}, \code{to}, \code{count}.
#'   \code{tda_seqevd} returns a data frame: \code{time}, \code{cases},
#'   one match-count column per event type found, and \code{total}.
#' @family sequence analysis
#' @examples
#' # three individuals, six time points (0-5), states 1-3; the target
#' # event is a 1 -> 2 transition, with two time-independent covariates
#' s <- data.frame(Y0 = c(1, 1, 2), Y1 = c(1, 2, 1), Y2 = c(2, 1, 1),
#'                 Y3 = c(1, 1, 2), Y4 = c(1, 1, 2), Y5 = c(3, 3, 3))
#' cov <- data.frame(v1 = c(1, 2, 2), v2 = c(-10, -15, -20))
#'
#' tda_seqev(s)    # which transitions actually occur, and how often
#' tda_seqevd(s)   # event counts over time
#'
#' d <- tda_seqmd(s, event = c(1, 2), covariates = cov)
#' d
#' # fit it: a period dummy for every distinct time point, the way
#' # factor(time) would in an ordinary logistic regression
#' pd <- grep("^period", names(d), value = TRUE)
#' tda_qreg(stats::as.formula(paste("event ~", paste(c(pd, "v1", "v2"),
#'                                                    collapse = "+"))),
#'         d, model = "logit", intercept = FALSE)
#'
#' # the same data, but with an ID column and the states/covariates all
#' # together in one data frame, the more realistic shape -- sequences=
#' # and covariates= just name columns in data= directly
#' full <- cbind(id = 1:3, s, cov)
#' tda_seqmd(paste0("Y", 0:5), event = c(1, 2), data = full,
#'          covariates = c("v1", "v2"))
#' @export
tda_seqmd <- function(sequences, event, data = NULL, covariates = NULL,
                      event_covariates = NULL, select = NULL, tp = NULL,
                      summary = FALSE, options = list(),
                      dir = tempfile("tda"), ...) {
    d <- .seq_recode(sequences, data)
    nm <- names(d)
    if (length(event) != 2L)
        stop("`event` must be c(from, to)")

    opts <- list(ev = sprintf("[1,%d,%d]", event[1L], event[2L]))
    if (!is.null(select)) {
        # select= is a raw string, never translated the way definitions=
        # elsewhere in this package is -- so it has to be rewritten here
        # from the real column names given in sequences= to the internal
        # Y1, Y2, ... TDA was actually given, the same problem tda_frml's
        # own Surv() arguments had and for the identical reason: a
        # TDA syntax error otherwise, checked, not
        # merely a naming preference. R's comparison operators
        # (Y0 == 1, Y0 = 1) are also accepted here, translated to TDA's
        # own eq(Y0,1) form before that column-name rewrite -- see
        # .select_translate_ops().
        select <- .select_translate_ops(select)
        orig <- if (is.character(sequences)) sequences
               else names(as.data.frame(sequences))
        for (i in seq_along(orig))
            if (orig[i] != nm[i])
                select <- gsub(sprintf("\\b%s\\b", orig[i]), nm[i], select)
        opts$sel <- select
    }
    if (!is.null(tp))
        opts$tp <- if (length(tp) > 1L) paste(tp, collapse = ",") else tp

    # summary = TRUE is a separate mode of TDA's seqmd,
    # not a different reading of the same output: called without a
    # right-hand-side file at all, it prints a Time/RiskSet/Events table
    # to the console instead of writing the full per-case, per-period
    # records seqmd otherwise writes -- checked against a
    # real run (with no rhs, no file gets written, and this table is all
    # that appears). covariates=/event_covariates= have nothing to
    # contribute to that table (TDA's console summary does not carry
    # them), so they are not accepted here.
    if (isTRUE(summary)) {
        if (!is.null(covariates) || !is.null(event_covariates))
            stop("summary = TRUE does not support covariates= or ",
                 "event_covariates= -- TDA's console summary has no ",
                 "covariate columns to show")
        opts <- c(opts, .tda_extra(options))
        res <- tda_run(c(tda_nvar(d), sprintf("seqdef = %s;",
                                              paste(nm, collapse = ",")),
                         do.call(tda_block, c(list(name = "seqmd"), opts))),
                       data = d, dir = dir, ...)
        err <- grep("^Error", res$output, value = TRUE)
        if (length(err))
            stop("TDA could not run seqmd: ", err[1L], call. = FALSE)
        i <- grep("^Time\\s+Risk Set\\s+Events", res$output)
        if (!length(i))
            return(NULL)
        body <- res$output[(i[1L] + 2L):length(res$output)]
        # The table ends at the next blank line or the "---" rule TDA
        # prints after every command's output -- reading to the end
        # of the console otherwise picks up trailing lines ("Current
        # memory: ..." and the like) as spurious extra, malformed rows.
        end <- which(!nzchar(trimws(body)) | grepl("^-+$", trimws(body)))
        if (length(end))
            body <- body[seq_len(end[1L] - 1L)]
        body <- body[nzchar(trimws(body))]
        toks <- strsplit(trimws(body), "\\s+")
        m <- do.call(rbind, toks)
        smy <- data.frame(time = as.numeric(m[, 1L]),
                          riskset = as.numeric(m[, 2L]),
                          events = as.numeric(m[, 3L]))
        smy <- .overlay_num(smy, res$exports[["seqmd.risk"]])
        attr(smy, "run") <- res
        return(smy)
    }

    covnm <- NULL
    if (!is.null(covariates)) {
        # covariates= can also just name columns in data=, the same way
        # sequences= can, rather than requiring a separately-built data
        # frame for what is usually already sitting right there.
        cov <- if (!is.null(data) && is.character(covariates))
            data[covariates] else as.data.frame(covariates)
        orig_covnames <- names(cov)
        covnm <- .tda_names(orig_covnames)
        names(cov) <- covnm
        d <- cbind(d, cov)
        opts$v <- paste(covnm, collapse = ",")
    }
    # event_covariates= (xe=): unlike v=, a time-varying
    # covariate -- a separate value per time point, not one fixed value
    # per case -- given as a named list, each element the per-time-point
    # column names for one such covariate (in the same order as the
    # sequence states themselves), checked against a real
    # run: xe= takes one bracketed group per covariate ([S0,S1,...,S5]),
    # not a plain comma list the way v= does, and reports back "Number of
    # event variables: 1" for one such group, not one per column
    # inside it.
    if (!is.null(event_covariates)) {
        if (is.null(names(event_covariates)) || any(!nzchar(names(event_covariates))))
            stop("`event_covariates` must be a named list, one name per ",
                 "time-varying covariate")
        xecols <- unlist(event_covariates, use.names = FALSE)
        if (is.null(data))
            stop("`event_covariates` needs `data` to read its columns from")
        xenm <- .tda_names(xecols)
        # a column already in the frame (a covariate named in both v=
        # and xe=, the manual's xe=[V1],[V2] of 6.18.3) is not added twice
        new <- !(xenm %in% names(d))
        if (any(new)) {
            xed <- data[xecols[new]]
            names(xed) <- xenm[new]
            d <- cbind(d, xed)
        }
        pos <- 1L
        groups <- vector("character", length(event_covariates))
        for (k in seq_along(event_covariates)) {
            n_k <- length(event_covariates[[k]])
            groups[k] <- paste0("[", paste(xenm[pos:(pos + n_k - 1L)],
                                           collapse = ","), "]")
            pos <- pos + n_k
        }
        opts$xe <- paste(groups, collapse = ",")
    }
    opts$dtda <- "desc.tda"
    opts <- c(opts, .tda_extra(options))

    res <- tda_run(c(tda_nvar(d), sprintf("seqdef = %s;",
                                          paste(nm, collapse = ",")),
                     do.call(tda_block, c(list(name = "seqmd"), opts,
                                          list(rhs = "out.d")))),
                   data = d, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run seqmd: ", err[1L], call. = FALSE)
    # export first, so out.d is not read when seqmd.table covers it
    etab0 <- if (.use_exports()) res$exports[["seqmd.table"]]
    tab <- if (is.matrix(etab0)) .export_frame(etab0)
           else tda_file(res, "out.d")
    if (is.null(tab))
        stop("seqmd produced no output")
    # dtda= is TDA's, reliable description of the output's real
    # column structure (id/time/target-event/period dummies, each with
    # its own "# comment" TDA itself writes) -- read back here rather
    # than assumed from column position and count, which is exactly
    # the kind of assumption that silently breaks the moment xe= adds
    # columns in a shape not otherwise verified. v='s columns keep
    # their real names directly (no comment to read for those); xe='s
    # own columns are named from event_covariates instead of TDA's
    # auto-generated "S0_D"-style name.
    # The export names every column, so nothing is appended after it;
    # the file names only the columns TDA commented, and the v=/xe=
    # ones are added below.
    desc_all <- .dtda_names_from_export(res, ncol(tab))
    desc <- desc_all %||% .read_dtda_names(file.path(res$dir, "desc.tda"))
    if (!is.null(desc) && length(desc) >= 3L) {
        # .read_dtda_names() returns TDA's variable names (ID, TIME,
        # Z, P1, P2, ...), not the comment text next to them -- mapped
        # here to this function's own, more readable equivalents
        # directly, since that structure (id/time/target-event/period
        # dummies, always in that order) is itself the reliable part
        # confirmed by dtda=, not something to re-derive generically.
        renamed <- ifelse(desc == "ID", "id",
                   ifelse(desc == "TIME", "time",
                   ifelse(desc == "Z", "event",
                   ifelse(grepl("^P[0-9]+$", desc),
                         paste0("period", sub("^P", "", desc)), desc))))
        names(tab)[seq_along(desc)] <- renamed
        if (is.null(desc_all)) {
            rest <- seq.int(length(desc) + 1L, ncol(tab))
            rest_names <- character(0)
            if (length(covnm)) rest_names <- c(rest_names, orig_covnames)
            if (!is.null(event_covariates))
                rest_names <- c(rest_names, names(event_covariates))
            if (length(rest) == length(rest_names))
                names(tab)[rest] <- rest_names
        }
        else {
            # The export names every column, but with TDA's
            # spellings: a v= covariate keeps the mangled name it was
            # given, and an xe= column gets an auto-generated one
            # ("V1_D").  The caller asked for its own names, and the
            # file path substitutes them, so this one does too --
            # v= by name, xe= positionally at the end, which is the
            # order the file path appends them in.
            if (length(covnm)) {
                i <- match(covnm, names(tab))
                names(tab)[i[!is.na(i)]] <- orig_covnames[!is.na(i)]
            }
            if (!is.null(event_covariates)) {
                # the caller's names go on the columns TDA made for the
                # xe= groups it named (the "_D" ones), in order; any
                # further xe= output (the event counters E_, ED_, E0_..
                # of the counting form) keeps TDA's names
                k <- length(event_covariates)
                xd <- which(grepl("_D$", desc_all))
                if (length(xd) >= k)
                    names(tab)[xd[seq_len(k)]] <- names(event_covariates)
            }
        }
    } else {
        # dtda= itself failed to produce a readable description -- fall
        # back to the same position-based naming used before dtda= was
        # read at all, rather than leave the columns unnamed.
        nper <- ncol(tab) - 3L - length(covnm) -
            if (!is.null(event_covariates)) length(event_covariates) else 0L
        names(tab) <- c("id", "time", "event", paste0("period", seq_len(nper)),
                        if (length(covnm)) orig_covnames,
                        if (!is.null(event_covariates)) names(event_covariates))
    }
    tab <- .overlay_num(tab, res$exports[["seqmd.table"]])
    attr(tab, "run") <- res
    tab
}

# Shared with tda_seqm(): recode a sequence data frame's columns to
# integer states named Y1, Y2, ... -- the same conversion, kept in one
# place so both stay consistent.
.seq_recode <- function(sequences, data = NULL) {
    # sequences can name columns within a separate data= directly
    # (sequences = paste0("Y", 0:5), data = the_full_data_frame) rather
    # than requiring the caller to subset a real data frame to just the
    # state columns first -- checked that this matters: real
    # sequence data typically also carries an ID column and covariates
    # alongside the states in the same data frame (TDA's
    # examples/exam/seq.d4, e.g.), and requiring a pre-subsetted
    # data frame with no way to say "these columns, from this data" was
    # exactly what made it unclear how to use these functions at all.
    if (!is.null(data) && is.character(sequences)) {
        miss <- setdiff(sequences, names(data))
        if (length(miss))
            stop("`sequence` names columns not in `data`: ",
                 paste(miss, collapse = ", "))
        sequences <- data[sequences]
    }
    d <- as.data.frame(sequences)
    if (!ncol(d))
        stop("no sequence columns")
    chr <- vapply(d, function(z) is.character(z) || is.factor(z), logical(1))
    if (any(chr)) {
        lev <- sort(unique(unlist(lapply(d[chr], as.character))))
        for (j in which(chr))
            d[[j]] <- as.integer(factor(as.character(d[[j]]), levels = lev))
    }
    for (j in seq_along(d))
        d[[j]] <- as.integer(d[[j]])
    names(d) <- paste0("Y", seq_len(ncol(d)))
    d
}

# TDA's seqdef() can be called several times, sn=1, sn=2, ... to
# hold several independent sequence data structures at once, built from
# entirely different columns of the same underlying data -- confirmed
# directly against a real manual example (seq.d4's Y0..Y5 states as
# sn=1 and S0..S5 as a second, separate sequence, not a
# covariate, as sn=2). A single spec (a data frame, matrix, or
# character vector of column names) keeps working exactly as before,
# implicitly sn=1; sequences= given as a list of several such specs
# instead builds one seqdef() per element, sn=1, sn=2, ..., in order.
# Also returns orig_names, the real column name behind each internal
# SQ<i>_<j> one, per structure -- needed to translate select='s
# text from the user's real column names (confirmed a real,
# separate bug: select= is a raw string, never translated the way
# tda_frml's definitions= is, so a real column name in it was a
# TDA syntax error the moment this internal SQ<i>_<j> naming
# was introduced for multi-structure support).
.seq_recode_multi <- function(sequences, data = NULL, type = 1) {
    multi <- is.list(sequences) && !is.data.frame(sequences)
    if (!multi)
        sequences <- list(sequences)
    orig_names <- lapply(sequences, function(s) {
        if (is.character(s)) s else names(as.data.frame(s))
    })
    specs <- lapply(sequences, .seq_recode, data = data)
    n <- length(specs)
    for (i in seq_len(n))
        names(specs[[i]]) <- paste0("SQ", i, "_", seq_along(specs[[i]]))
    d <- do.call(cbind, specs)
    # type 2 is seqdef's m=2: the columns are state/time pairs, and
    # the time axis comes from the times rather than the column count
    # type may be given per structure: the manual's seq3.cf declares
    # a type 1 structure over Y0..Y7 and a type 2 over three state/time
    # pairs, and its output shows both in one table (3.4.2, Box 10).
    ty <- rep_len(as.integer(type), n)
    defs <- vapply(seq_len(n), function(i)
        sprintf(if (ty[i] == 2L) "seqdef(sn=%d,m=2) = %s;"
                else "seqdef(sn=%d) = %s;",
                i, paste(names(specs[[i]]), collapse = ",")),
        character(1))
    list(data = d, defs = defs, n = n, orig_names = orig_names,
        int_names = lapply(specs, names))
}

#' Information about currently defined sequence data structures
#'
#' Runs TDA's \code{seq;} against the sequence structure(s) given --
#' the same table TDA itself prints after defining them, listing each
#' structure's number of state variables, its time axis, and the
#' distinct states it actually contains.
#'
#' @inheritParams tda_seqmd
#' @param type the sequence data type, TDA's \code{seqdef} \code{m=}:
#'   \code{1} (default) one column per time point, \code{2} columns read
#'   as state/time pairs, one pair per spell.
#' @return A data frame, one row per sequence structure: \code{sn},
#'   \code{type}, \code{variables} (number of time points),
#'   \code{tmin}/\code{tmax} (the time axis), \code{nstates}, and
#'   \code{states} (a list column, the actual state values, since their
#'   count varies by structure).
#' @family sequence analysis
#' @examples
#' s <- data.frame(Y0 = c(1, 1, 2), Y1 = c(1, 2, 1), Y2 = c(2, 1, 1))
#' tda_seq_info(s)
#'
#' # several structures at once, the same way tda_seqev's sn= does
#' d <- data.frame(Y0 = c(1, 1, 2), Y1 = c(1, 2, 1), S0 = c(1, 3, 1),
#'                 S1 = c(3, 1, 3))
#' tda_seq_info(list(c("Y0", "Y1"), c("S0", "S1")), data = d)
#' @export

tda_seq_info <- function(sequences, data = NULL, options = list(),
                         type = 1, dir = tempfile("tda"), ...) {
    ms <- .seq_recode_multi(sequences, data, type = type)
    # seq -- with no arguments at all -- is a bare command, not a
    # tda_block()-style name(...) call (checked: "seq()"
    # itself is a "Unknown command" error); every options=
    # entry still has to reach TDA somehow, so anything given there is
    # simply ignored rather than silently miscompiled into invalid
    # syntax, since seq takes none.
    if (length(options))
        warning("`options` is ignored: seq takes no options of its own")
    res <- tda_run(c(tda_nvar(ms$data), ms$defs, "seq;"),
                   data = ms$data, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run seq: ", err[1L], call. = FALSE)
    # seqdef() itself auto-prints this same table after every call, not
    # only when seq is asked for explicitly -- checked, and
    # each one only reflects the structures defined up to that point, so
    # the *last* such table in the console output (after every seqdef()
    # and the explicit seq at the end) is the complete, correct one, not
    # the first, which the generic "find a matching header" approach
    # elsewhere in this file would otherwise pick.
    # built from seq.defs (one row per structure) and
    # seq.states (one block per structure, its state values), so
    # the console table below is only the tdaR.use_exports = FALSE path.
    e <- .seq_info_from_exports(res)
    if (!is.null(e))
        return(e)
    i <- grep("^Structure Type", res$output)
    if (!length(i))
        return(NULL)
    i <- i[length(i)]
    body <- res$output[(i + 2L):length(res$output)]
    end <- which(!nzchar(trimws(body)) | grepl("^-+$", trimws(body)))
    if (length(end))
        body <- body[seq_len(end[1L] - 1L)]
    body <- body[nzchar(trimws(body))]
    toks <- strsplit(trimws(body), "\\s+")
    fixed <- do.call(rbind, lapply(toks, function(t) t[1:6]))
    states <- lapply(toks, function(t) as.integer(t[-(1:6)]))
    out <- data.frame(sn = as.integer(fixed[, 1L]),
                      type = as.integer(fixed[, 2L]),
                      variables = as.integer(fixed[, 3L]),
                      tmin = as.numeric(fixed[, 4L]),
                      tmax = as.numeric(fixed[, 5L]),
                      nstates = as.integer(fixed[, 6L]),
                      states = I(states))
    # Every seqdef() prints one of these tables, so the exports hold one
    # block per table; the parser reads the last table, so the last
    # block is the matching one.  The state lists are flushed one per
    # structure, which is why the last nrow(out) of them are taken
    # rather than the last one.
    if (.use_exports()) {
        db <- .exports_blocks(res$exports, "seq.defs")
        if (length(db)) {
            m <- db[[length(db)]]
            if (is.matrix(m) && nrow(m) == nrow(out) && ncol(m) == 6L)
                out <- .overlay_cols(out, m,
                                     c("sn", "type", "variables", "tmin",
                                       "tmax", "nstates"))
        }
        sb <- .exports_blocks(res$exports, "seq.states")
        if (length(sb) >= nrow(out)) {
            sb <- sb[seq.int(length(sb) - nrow(out) + 1L, length(sb))]
            if (all(vapply(seq_along(sb), function(k)
                    is.matrix(sb[[k]]) && length(sb[[k]]) ==
                        length(out$states[[k]]), NA)))
                out$states <- I(lapply(seq_along(sb), function(k)
                    .overlay_col(out$states[[k]], as.vector(sb[[k]]))))
        }
    }
    out
}

# seq.defs holds one row per structure (sn, type, variables, tmin,
# tmax, nstates) and seq.states one block per structure carrying its
# state values.  Every seqdef() call re-prints this table, so the
# exports hold one block per printed table too: the LAST nrow(defs)
# state blocks belong to the last, complete table -- the same one the
# console parser deliberately takes.
.seq_info_from_exports <- function(res) {
    if (!.use_exports())
        return(NULL)
    db <- .exports_blocks(res$exports, "seq.defs")
    if (!length(db))
        return(NULL)
    m <- db[[length(db)]]
    if (!is.matrix(m) || ncol(m) != 6L)
        return(NULL)
    sb <- .exports_blocks(res$exports, "seq.states")
    if (length(sb) < nrow(m))
        return(NULL)
    sb <- sb[seq.int(length(sb) - nrow(m) + 1L, length(sb))]
    states <- lapply(sb, function(x)
        if (is.matrix(x)) as.integer(as.vector(x)) else NULL)
    if (any(vapply(states, is.null, NA)))
        return(NULL)
    if (!identical(vapply(states, length, 0L), as.integer(m[, 6L])))
        return(NULL)
    data.frame(sn = as.integer(m[, 1L]), type = as.integer(m[, 2L]),
               variables = as.integer(m[, 3L]),
               tmin = as.numeric(m[, 4L]), tmax = as.numeric(m[, 5L]),
               nstates = as.integer(m[, 6L]),
               states = I(states))
}

# Rewrites select='s text from the real column names given for
# structure sn to the internal SQ<sn>_<j> names TDA was actually given
# -- the same kind of translation tda_frml's definitions= needed
# for Surv()'s arguments, and for the identical reason: a raw string
# handed to TDA verbatim, with no other mechanism to make the user's
# own column names reach it.
# Lets select= be written with R's comparison operators (Y0 == 1,
# or Y0 = 1) instead of requiring TDA's eq(Y0,1) function-call form
# up front. A plain text substitution on purpose, not R's parser --
# select= has to keep accepting a bare "=" too (not valid R syntax to
# pass unquoted as an argument at all), and a real TDA expression
# (an existing eq(...)/ge(...) call, say) needs to pass through
# untouched rather than be reparsed and rewritten. Longer/composite
# operators are matched before shorter ones (== before the bare =, >=
# before >, ...) so a real == is never mistaken for two consecutive =
# matches.
.select_translate_ops <- function(select) {
    if (is.null(select))
        return(select)
    ident <- "[A-Za-z_][A-Za-z0-9_.]*"
    val <- "-?[A-Za-z0-9_.]+"
    ops <- list(c("==", "eq"), c("!=", "ne"), c(">=", "ge"), c("<=", "le"),
               c(">", "gt"), c("<", "lt"), c("=", "eq"))
    for (op in ops) {
        pat <- sprintf("(%s)\\s*%s\\s*(%s)", ident, op[1L], val)
        select <- gsub(pat, sprintf("%s(\\1,\\2)", op[2L]), select,
                       perl = TRUE)
    }
    select
}

.seq_translate_select <- function(select, ms, sn) {
    if (is.null(select))
        return(select)
    select <- .select_translate_ops(select)
    on <- ms$orig_names[[sn]]
    inm <- ms$int_names[[sn]]
    for (i in seq_along(on))
        if (on[i] != inm[i])
            select <- gsub(sprintf("\\b%s\\b", on[i]), inm[i], select)
    select
}

#' @rdname tda_seqmd
#' @export
tda_seqev <- function(sequences, data = NULL, sn = 1L, select = NULL,
                      options = list(), dir = tempfile("tda"), ...) {
    ms <- .seq_recode_multi(sequences, data)
    if (sn < 1L || sn > ms$n)
        stop("`sn` must be between 1 and ", ms$n,
             " (the number of sequence structures given)")
    opts <- .tda_extra(options)
    if (!is.null(select))
        opts$sel <- .seq_translate_select(select, ms, sn)
    if (ms$n > 1L)
        opts$sn <- sn
    # seqev's syntax requires a right-hand-side file name (confirmed
    # directly against its help text and a real run: without one, it is
    # a TDA syntax error, not a fallback to console output) --
    # results are only ever written to that file, never printed to the
    # console at all, so there was never any "Event type ... Number"
    # text for the old console-parsing version of this function to find.
    res <- tda_run(c(tda_nvar(ms$data), ms$defs,
                     do.call(tda_block, c(list(name = "seqev"), opts,
                                          list(rhs = "out.d")))),
                   data = ms$data, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run seqev: ", err[1L], call. = FALSE)
    etab0 <- if (.use_exports()) res$exports[["seqev.table"]]
    tab <- if (is.matrix(etab0)) .export_frame(etab0)
           else tda_file(res, "out.d")
    if (is.null(tab))
        stop("seqev produced no output")
    names(tab) <- c("from", "to", "count")[seq_len(ncol(tab))]
    tab <- .overlay_num(tab, res$exports[["seqev.table"]])
    attr(tab, "run") <- res
    tab
}

#' @rdname tda_seqmd
#' @export
tda_seqevd <- function(sequences, data = NULL, sn = 1L, select = NULL,
                       options = list(), dir = tempfile("tda"), ...) {
    ms <- .seq_recode_multi(sequences, data)
    if (sn < 1L || sn > ms$n)
        stop("`sn` must be between 1 and ", ms$n,
             " (the number of sequence structures given)")
    opts <- .tda_extra(options)
    if (!is.null(select))
        opts$sel <- .seq_translate_select(select, ms, sn)
    if (ms$n > 1L)
        opts$sn <- sn
    opts$dtda <- "desc.tda"
    res <- tda_run(c(tda_nvar(ms$data), ms$defs,
                     do.call(tda_block, c(list(name = "seqevd"), opts,
                                          list(rhs = "out.d")))),
                   data = ms$data, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run seqevd: ", err[1L], call. = FALSE)
    # export first, so out.d is not read when seqevd.table covers it
    etab0 <- if (.use_exports()) res$exports[["seqevd.table"]]
    tab <- if (is.matrix(etab0)) .export_frame(etab0)
           else tda_file(res, "out.d")
    if (is.null(tab))
        stop("seqevd produced no output")
    # dtda= names each column by the transition it actually counts
    # (EV1_2, EV2_1, ...), not a generic position -- checked
    # (TDA's comment: "number of events (1,2)"), and this is real
    # information a generic "events1"/"events2" throws away: which
    # transition each column is, not just how many there were.
    desc <- .dtda_names_from_export(res, ncol(tab)) %||%
        .read_dtda_names(file.path(res$dir, "desc.tda"))
    if (!is.null(desc) && length(desc) == ncol(tab)) {
        names(tab) <- ifelse(desc == "TIME", "time",
                      ifelse(desc == "NCAS", "cases",
                      ifelse(desc == "NEV", "total", tolower(desc))))
    } else {
        # dtda= itself failed to produce a readable description -- fall
        # back to the same generic, position-based naming used before
        # dtda= was read at all, rather than leave the columns unnamed.
        ntype <- ncol(tab) - 3L
        names(tab) <- c("time", "cases",
                        paste0("events", seq_len(ntype)), "total")
    }
    tab <- .overlay_num(tab, res$exports[["seqevd.table"]])
    attr(tab, "run") <- res
    tab
}

# Several sequence-description commands (seqgc, seqlg, seqsd, seqsi) name
# their output columns dynamically -- seqgc's depend on how many
# distinct states the data actually has, for instance -- so rather than
# guess or hardcode a column count, this reads them back from the dtda
# description file every one of these commands can write, taking the
# name TDA itself gave each column from its own "# comment" annotation.
.read_dtda_names <- function(path) {
    if (!file.exists(path))
        return(NULL)
    lines <- readLines(path, warn = FALSE)
    m <- regmatches(lines, regexec("^\\s*(\\S+)\\s*<.*#\\s*(.*)$", lines))
    nm <- vapply(m, function(x) if (length(x) >= 2L) x[2L] else NA_character_,
                 character(1))
    nm[!is.na(nm)]
}

# The column names TDA wrote into its dtda description, as a typed
# string vector straight off the export channel -- the same channel
# coeff.names and rate.est.names use.  No text is reconstructed and no
# regex runs here: the producer recognises its writer's fixed
# "NAME <w>[fmt] = cK" shape at the point of writing, and hands over
# the names themselves.
.dtda_names_from_export <- function(res, ncol = NULL) {
    if (is.null(res) || !.use_exports())
        return(NULL)
    nm <- res$exports[["dtda.names"]]
    if (!is.character(nm) || !length(nm))
        return(NULL)
    # The export names EVERY column, including the v=/covariate ones
    # whose description lines carry no "# comment" for the file
    # reader's regex to match.  A caller passes its column count
    # and gets names only when they cover the table exactly -- so a
    # caller using these never has to append anything, and a shape it
    # does not recognise falls back to the file rather than naming
    # columns wrongly.
    if (!is.null(ncol) && length(nm) != ncol)
        return(NULL)
    nm
}

#' Sequence characteristics, lengths, gaps and state distributions
#'
#' Descriptive summaries of a sequence data structure: \code{tda_seqgc}
#' basic characteristics per case (length, number of states visited,
#' number of transitions, time spent in and number of episodes in each
#' state), \code{tda_seqlg} sequence lengths and gaps, \code{tda_seqsd}
#' the distribution of states over time (how many cases are in each state
#' at each time point), \code{tda_seqsi} a state indicator matrix (one
#' row per case-time combination, one column per state, 1 where that case
#' is in that state at that time), and \code{tda_seqen} entropy of the
#' state distribution at each time point. Column names come from TDA's
#' own description of each output, not a fixed layout guessed in advance
#' -- \code{tda_seqgc}'s columns, for instance, depend on how many
#' distinct states the data actually has.
#'
#' @param sequences a data frame or matrix of one column per time point,
#'   in order, the same shape \code{\link{tda_seqm}} takes -- optionally
#'   with further, non-sequence columns (an id, say) alongside them; use
#'   \code{variables}/\code{id} to say which is which rather than
#'   pre-subsetting \code{sequences} yourself.
#' @param tp for \code{tda_seqsi}, the time points to build the indicator
#'   matrix at -- required, as a vector or a single TDA range expression
#'   such as \code{"0(1)10"}.
#' @param select optional expression (\code{sel=}) restricting which
#'   sequences are used -- written with R's comparison operators
#'   (\code{Y0 == 1}, or \code{Y0 = 1}) or TDA's function form
#'   (\code{eq(Y0,1)}) interchangeably; \code{&}/\code{|} combine several
#'   conditions the same way in both.
#' @param variables which columns of \code{sequences} are the actual
#'   time-point columns, in order -- TDA's \code{seqdef=Y0,,Y7}, which
#'   likewise names the sequence out of a data matrix holding other
#'   variables too rather than requiring them to be dropped first. A
#'   character vector of names in \code{sequences}; defaults to every
#'   column (so passing an already-subsetted \code{sequences} still works
#'   exactly as before).
#' @param id further columns of \code{sequences} to carry through to the
#'   output as-is, one value per case -- \code{seqgc}/\code{seqlg}/
#'   \code{seqpm}/\code{seqsi}'s \code{v=}, which writes the case's
#'   raw value onto its output row, no more.
#'   A character vector of names in \code{sequences}. \code{seqsd}
#'   and \code{seqen} have no per-case row to attach a value to (one row
#'   per time point, over all cases) and error if given one --
#'   TDA accepts \code{v=} there without
#'   complaint but never actually writes the column.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @param type the sequence data type, TDA's \code{seqdef} \code{m=}:
#'   \code{1} (default) one column per time point, \code{2} columns read
#'   as state/time pairs, one pair per spell, with the time axis taken
#'   from the times rather than from the column count.
#' @return A data frame, one row per case (\code{tda_seqgc},
#'   \code{tda_seqlg}) or per time point (\code{tda_seqsd},
#'   \code{tda_seqen}), or per case-time combination (\code{tda_seqsi}).
#' @family sequence analysis
#' @examples
#' s <- data.frame(Y0 = c(1, 1, 2), Y1 = c(1, 2, 1), Y2 = c(2, 1, 1),
#'                 Y3 = c(1, 1, 2), Y4 = c(1, 1, 2), Y5 = c(3, 3, 3))
#' tda_seqgc(s)   # length, states visited, transitions, time in each state
#' tda_seqlg(s)   # lengths and gaps
#' tda_seqsd(s)   # how many cases are in each state, at each time point
#' tda_seqsi(s, tp = "0(1)5")   # state indicator matrix
#' tda_seqen(s)   # entropy of the state distribution at each time point
#'
#' # an id column alongside the sequence, positively selected rather
#' # than requiring s2[, -1]
#' s2 <- cbind(id = 101:103, s)
#' tda_seqgc(s2, variables = names(s), id = "id")
#' @export
tda_seqgc <- function(sequences, select = NULL, variables = NULL,
                      id = NULL, options = list(), type = 1, dir = tempfile("tda"),
                      ...)
    .seq_desc("seqgc", sequences, variables, id, list(sel = select),
              options, dir, allow_id = TRUE, type = type, ...)

#' @rdname tda_seqgc
#' @export
tda_seqlg <- function(sequences, select = NULL, variables = NULL,
                      id = NULL, options = list(), type = 1, dir = tempfile("tda"),
                      ...)
    .seq_desc("seqlg", sequences, variables, id, list(sel = select),
              options, dir, allow_id = TRUE, ...)

#' @rdname tda_seqgc
#' @export
tda_seqsd <- function(sequences, select = NULL, variables = NULL,
                      options = list(), dir = tempfile("tda"), ...)
    .seq_desc("seqsd", sequences, variables, NULL, list(sel = select),
              options, dir, allow_id = FALSE, ...)

#' Write sequence data back out
#'
#' \code{seqpd} prints the sequence data in one of four layouts, chosen
#' by \code{layout} (TDA's \code{m=}):
#' \describe{
#'   \item{1}{every sequence structure side by side, one row per case}
#'   \item{2}{one row per case and structure}
#'   \item{3}{one row per case, structure and position}
#'   \item{4}{one row per spell: its origin and destination state, and
#'     the times it spans}
#' }
#'
#' @inheritParams tda_seqgc
#' @param layout 1 to 4, the shape of the result, as listed above.
#' @param type the sequence data type, TDA's \code{seqdef} \code{m=}:
#'   \code{1} (default) one column per time point, \code{2} columns read
#'   as state/time pairs, one pair per spell, with the time axis taken
#'   from the times rather than from the column count.
#' @param second names of the columns of a second sequence structure,
#'   declared beside the first as the manual's \code{seq5.cf} does; the
#'   tables then carry one column per structure.
#' @param second_type the sequence data type of the second structure,
#'   as \code{type}; default 2, state/time pairs.
#' @return A data frame in the requested layout.
#' @examples
#' # four cases over three states, one sequence structure
#' s <- data.frame(ID = 1:4,
#'                 Y1 = c(1, 2, 3, 1), Y2 = c(1, 2, 3, 3),
#'                 Y3 = c(2, 2, 1, 3), Y4 = c(2, 3, 1, 1))
#' tda_seqpd(s, id = "ID")
#'
#' # one row per spell instead, with the times each spans
#' tda_seqpd(s, layout = 4, id = "ID")
#' @export
tda_seqpd <- function(sequences, layout = 1, select = NULL,
                      variables = NULL, id = NULL, options = list(),
                      type = 1, second = NULL, second_type = 2,
                      dir = tempfile("tda"), ...) {
    if (!layout %in% 1:4)
        stop("`layout` is 1, 2, 3 or 4 -- see ?tda_seqpd", call. = FALSE)
    .seq_desc("seqpd", sequences, variables, id,
              list(m = layout, sel = select), options, dir,
              allow_id = TRUE, type = type, second = second,
              second_type = second_type, ...)
}

#' @rdname tda_seqgc
#' @export
tda_seqsi <- function(sequences, tp, select = NULL, variables = NULL,
                      id = NULL, options = list(), dir = tempfile("tda"),
                      ...) {
    if (missing(tp))
        stop("`tp` is required: the time points to build the indicator ",
             "matrix at")
    .seq_desc("seqsi", sequences, variables, id,
             list(tp = if (length(tp) > 1L) paste(tp, collapse = ",") else tp,
                  sel = select),
             options, dir, allow_id = TRUE, ...)
}

#' @rdname tda_seqgc
#' @export
tda_seqen <- function(sequences, select = NULL, variables = NULL,
                      options = list(), dir = tempfile("tda"), ...)
    .seq_desc("seqen", sequences, variables, NULL, list(sel = select),
              options, dir, allow_id = FALSE, ...)

# variables= positively selects, from `sequences`, the columns that make
# up the actual time-point sequence -- TDA's seqdef=Y0,,Y7 does the
# same thing against its data matrix, naming the sequence rather than
# requiring every other variable to be absent. Defaulting to every
# column of `sequences` keeps an already-subsetted caller working
# unchanged. id= is the separate, v= passthrough: further named
# columns whose raw per-case value TDA writes onto the output row as-is
# (t_seq.c's "additional variables" loops, e.g. seqgc's), not usable on
# seqsd/seqen, which have no per-case row at all -- checked:
# TDA accepts v= there without error but the dtda file it writes claims
# a column its output never contains, silently.
.seq_select <- function(sequences, variables, id, allow_id, cmd) {
    sdf <- as.data.frame(sequences)
    if (!is.null(id)) {
        if (!allow_id)
            stop(cmd, " has no `id`: its output is one row per time ",
                 "point, over all cases, not one row per case -- there ",
                 "is nothing to attach a per-case value to (TDA accepts ",
                 "`v=` here without error, but never actually writes it)")
        if (!is.character(id))
            stop("`id` must be a character vector of column names in ",
                 "`sequences`")
        miss <- setdiff(id, names(sdf))
        if (length(miss))
            stop("`id` not found in `sequences`: ",
                 paste(miss, collapse = ", "))
    }
    vars <- if (is.null(variables)) setdiff(names(sdf), id) else variables
    if (!is.character(vars))
        stop("`variables` must be a character vector of column names in ",
             "`sequences`")
    miss <- setdiff(vars, names(sdf))
    if (length(miss))
        stop("`variables` not found in `sequences`: ",
             paste(miss, collapse = ", "))
    dup <- intersect(vars, id)
    if (length(dup))
        stop("named in both `variables` and `id`: ",
             paste(dup, collapse = ", "))
    list(sequences = sdf[vars], id = if (!is.null(id)) sdf[id] else NULL)
}

.seq_desc <- function(cmd, sequences, variables, id, fixed_opts, options,
                      dir, allow_id, type = 1, second = NULL,
                      second_type = 2, ...) {
    picked <- .seq_select(sequences, variables, id, allow_id, cmd)
    d <- .seq_recode(picked$sequences)
    orig_n <- ncol(d)
    # select= is a raw string, never translated the way definitions=
    # elsewhere in this package is -- rewritten here the same way as
    # tda_seqev/tda_seqevd/tda_seqmd: R's comparison operators
    # (Y0 == 1, Y0 = 1) to TDA's eq(Y0,1) form, and the real column
    # names given in `sequences` to the internal Y1, Y2, ... TDA was
    # actually given (the same problem tda_frml's Surv() arguments
    # had, for the identical reason -- a TDA syntax error
    # otherwise, checked, not merely a naming preference).
    if (!is.null(fixed_opts$sel)) {
        sel <- .select_translate_ops(fixed_opts$sel)
        on <- names(picked$sequences)
        for (i in seq_along(on))
            if (on[i] != names(d)[i])
                sel <- gsub(sprintf("\\b%s\\b", on[i]), names(d)[i], sel)
        fixed_opts$sel <- sel
    }
    # v= writes further, non-sequence columns (an id, say) alongside the
    # command's output -- checked: by
    # comparing the output file with and without it: the requested
    # column appears as an extra field on every row. These are

    # not part of the sequence itself, so they are cbind()ed onto the
    # already-recoded sequence data rather than threaded through
    # .seq_recode(), which only knows about sequence columns. The
    # original column count is captured before the cbind() (not derived
    # from sequences itself, which is not reliably the same shape as d
    # for every caller) so the seqdef block only ever names the real
    # sequence columns, not the appended extras.
    vnm <- NULL
    if (!is.null(picked$id)) {
        vnm <- .tda_names(names(picked$id))
        d <- cbind(d, stats::setNames(picked$id, vnm))
    }
    opts <- fixed_opts[!vapply(fixed_opts, is.null, logical(1))]
    if (!is.null(vnm))
        opts$v <- paste(vnm, collapse = ",")
    opts <- c(opts, .tda_extra(options), list(dtda = "desc.t"))
    # TDA has two sequence data types.  Type 1 is one column per
    # time point, which every caller here has always sent.  Type 2
    # -- seqdef's m=2 -- reads the columns as state/time pairs,
    # one pair per spell, and builds the time axis from the times
    # rather than from the column count: the manual's section 3.4.2
    # Box 6 shows the same six columns read as three variables over
    # 0 to 8 instead of six over 0 to 5.
    seqdef_cmd <- if (identical(as.integer(type), 2L))
                      "seqdef(m = 2) = %s;"
                  else "seqdef = %s;"
    # A run may declare more than one structure, and the commands that
    # print or compare sequences then work across all of them: the
    # manual's seqpd example (3.4.3) defines a type 1 structure over
    # Y0..Y7 and a type 2 structure over three state/time pairs, and its
    # tables carry a column per structure. second= names the columns of
    # the second one; they are carried into the data unrecoded, since
    # the pairs are read as they stand.
    second_cols <- NULL
    if (!is.null(second)) {
        sec <- if (is.character(second)) sequences[, second, drop = FALSE]
               else as.data.frame(second)
        second_cols <- .tda_names(paste0("S2_", names(sec)))
        names(sec) <- second_cols
        d <- cbind(d, sec)
    }
    defs <- sprintf(seqdef_cmd, paste(names(d)[seq_len(orig_n)],
                                     collapse = ","))
    if (!is.null(second_cols))
        defs <- c(defs, sprintf("seqdef(sn = 2, m = %d) = %s;",
                                as.integer(second_type),
                                paste(second_cols, collapse = ",")))
    res <- tda_run(c(tda_nvar(d), defs,
                     do.call(tda_block, c(list(name = cmd), opts,
                                          list(rhs = "out.d")))),
                   data = d, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run ", cmd, ": ", err[1L], call. = FALSE)
    # the export is tried FIRST, so out.d is not read at all
    # when the producer covers this command -- the reader used to parse
    # the file and then overlay the exported values onto it, which meant
    # deleting the run directory broke a table whose every value was
    # already in hand.
    etab <- if (.use_exports())
        res$exports[[paste0(cmd, ".table")]]
    tab <- if (is.matrix(etab)) as.data.frame(etab) else tda_file(res, "out.d")
    if (is.null(tab))
        stop(cmd, " produced no output")
    # The export already names v='s appended extras; the FILE does not,
    # because TDA writes those lines without the "# comment" the file
    # reader matches on, so they are appended by hand in that path.
    nm_all <- .dtda_names_from_export(res, ncol(tab))
    nm <- nm_all %||% .read_dtda_names(file.path(res$dir, "desc.t"))
    if (is.null(nm_all) && !is.null(vnm) &&
        length(nm) == ncol(tab) - length(vnm))
        nm <- c(nm, names(picked$id))
    else if (!is.null(nm_all) && !is.null(vnm)) {
        i <- match(vnm, nm)
        nm[i[!is.na(i)]] <- names(picked$id)[!is.na(i)]
    }
    if (length(nm) == ncol(tab))
        names(tab) <- tolower(nm)
    # These commands write their table to a file at the print format,
    # so the numbers come back at that precision; where a producer
    # exists the exported doubles replace them.  The run itself is
    # attached the way the other file-backed returns carry it, so the
    # output is reachable from a result that is otherwise a bare frame.
    if (!is.matrix(etab))
        tab <- .overlay_num(tab, res$exports[[paste0(cmd, ".table")]])
    attr(tab, "run") <- res
    tab
}

#' Derive a sequence data structure from episode data
#'
#' \code{seqpe}: builds a sequence -- a state at every time point -- from
#' episode data (a start time, an end time, and origin/destination states
#' for each spell). The result can be read back with a fresh call to
#' \code{\link{tda_seqm}} or any of the other sequence functions once
#' written out and re-read; this function only performs the conversion,
#' since TDA's output is a plain data file rather than something kept
#' in memory across calls the way \code{\link{tda_seqm}} manages
#' internally.
#'
#' @param data a data frame of episodes, one row per spell.
#' @param id name of the case identifier column in \code{data}.
#' @param origin,destination names of the origin- and destination-state
#'   columns.
#' @param start,end names of the start- and end-time columns.
#' @param tp definition of the time points the sequence is built over, as
#'   a vector or a single TDA range expression such as \code{"0(1)10"}.
#' @param missing value to use for a gap at the beginning
#'   (\code{missing_start}), at the end (\code{missing_end}), or anywhere
#'   else (\code{missing}); TDA's default is \code{-1} for all three.
#' @param missing_start,missing_end see \code{missing}.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A data frame: the case id, then one column per time point in
#'   \code{tp}, each holding the state the case was in at that time.
#' @family sequence analysis
#' @examples
#' # two cases, two spells each
#' e <- data.frame(id = c(1, 1, 2, 2), org = c(1, 2, 1, 3),
#'                 des = c(2, 3, 3, 1), ts = c(0, 3, 0, 2),
#'                 tf = c(3, 5, 2, 5))
#' tda_seqpe(e, id = "id", origin = "org", destination = "des",
#'          start = "ts", end = "tf", tp = "0(1)5")
#' @export
tda_seqpe <- function(data, id, origin, destination, start, end, tp,
                      missing = NULL, missing_start = NULL,
                      missing_end = NULL, options = list(),
                      dir = tempfile("tda"), ...) {
    d <- as.data.frame(data)
    lab <- names(d)
    names(d) <- .tda_names(lab)
    rename <- function(x) names(d)[match(x, lab)]

    opts <- list(id = rename(id), org = rename(origin),
                des = rename(destination), ts = rename(start),
                tf = rename(end),
                tp = if (length(tp) > 1L) paste(tp, collapse = ",") else tp)
    if (!is.null(missing)) opts$m <- missing
    if (!is.null(missing_start)) opts$m1 <- missing_start
    if (!is.null(missing_end)) opts$m2 <- missing_end
    opts <- c(opts, .tda_extra(options))

    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block, c(list(name = "seqpe"), opts,
                                          list(rhs = "out.d")))),
                   data = d, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run seqpe: ", err[1L], call. = FALSE)
    etab0 <- if (.use_exports()) res$exports[["seqpe.table"]]
    tab <- if (is.matrix(etab0)) .export_frame(etab0)
           else tda_file(res, "out.d")
    if (is.null(tab))
        stop("seqpe produced no output")
    names(tab) <- c("id", paste0("t", seq_len(ncol(tab) - 1L) - 1L))
    tab
}

#' Maximum likelihood for a user-defined transition rate model
#'
#' \code{frml}: like \code{\link{tda_fml}}, but for episode (survival)
#' data rather than plain cross-sectional rows -- \code{fn} is the
#' log-likelihood contribution of one \emph{episode}, and can refer to
#' \code{time} (the position within it) the way \code{\link{tda_rate}}'s
#' own \code{define} can. This is not \code{fml} under another name:
#' the same definitions that work with \code{fml} on plain rows fail
#' here with \dQuote{no episode data defined} until an actual episode
#' structure (\code{survival::Surv()}) is given.
#'
#' \code{fn}'s expression cannot reference the raw duration/status columns
#' from \code{data} directly -- \code{Surv()}'s columns are consumed
#' by \code{edef()} into TDA's internal episode representation, and are
#' not separately visible by name afterward. An
#' \code{fn} written against the original column names fails with those
#' names listed as \emph{parameters to estimate} (TDA's convention for any
#' identifier it does not otherwise recognise), not as a lookup error.
#' Use TDA's built-in accessors for an episode instead: \code{ts} and
#' \code{tf} (start and finish time), \code{org} and \code{des} (origin
#' and destination state), \code{time} (\code{tf - ts}, and the same
#' built-in \code{\link{tda_rate}}'s \code{define} can use). A
#' covariate, by contrast, is referenced by its own name from
#' \code{formula}, translated automatically the same way
#' \code{\link{tda_fml}} does it.
#'
#' @param formula a \code{Surv(...)} formula, the same as
#'   \code{\link{tda_rate}}'s.
#' @param data a data frame.
#' @param definitions a character vector of assignments, evaluated in
#'   order inside TDA, the last of which must define \code{fn} -- see
#'   \code{\link{tda_fml}} for the same convention on plain data, and the
#'   description above for the episode-specific variables available here. Or
#'   an unevaluated \code{\{ \}} block of plain R assignments instead,
#'   translated automatically the same way \code{tda_fml}'s
#'   \code{definitions} is.
#'   Names on the right that are neither a covariate from \code{formula}
#'   nor one of TDA's built-ins, and are not defined earlier, are
#'   parameters to be estimated. \code{formula}'s \code{Surv(...)}
#'   arguments can also be referenced directly by their real column
#'   name when they are a simple column reference (\code{Surv(TFP,
#'   DES)}'s \code{TFP}/\code{DES}, say, work the same as TDA's
#'   built-ins \code{tf - ts}/\code{des}) -- matching how TDA's
#'   manual writes this (\code{examples/exam/frml1.cf}, referencing its
#'   own \code{DES}/\code{DUR} columns directly, not through any
#'   built-in alias). A computed \code{Surv()} argument
#'   (\code{Surv(TFin - TStart + 1, DES)}) has no one column to
#'   translate to, so only TDA's built-ins reach it.
#' @param start optional named list of starting values, in the order the
#'   definitions introduce the parameters.
#' @param id,spell columns identifying the case and numbering its spells,
#'   for multi-episode data -- the same pair \code{\link{tda_rate}}
#'   takes; a case's spells then enter the likelihood together.
#'   \code{strata()} in the formula is refused rather than accepted:
#'   \code{fml} has no stratification mechanism and would silently
#'   ignore the term (fit per stratum instead).
#' @param constraints optional linear constraints on the parameters --
#'   \code{frml}'s \code{con=}, the identical mechanism and
#'   \code{bN}-by-position convention as \code{\link{tda_fml}}'s
#'   \code{constraints} (see there for the full explanation and a worked
#'   example of the naming trap); it works the same way here.
#' @param residuals ask TDA to also compute, per episode, its
#'   contribution to \code{fn} at the converged parameters -- \code{frml}'s
#'   own \code{pres=}, the identical mechanism as \code{\link{tda_fml}}'s
#'   own \code{residuals} (see there for what it actually is: not a
#'   classical observed-minus-fitted residual), computed per episode
#'   here.
#' @param residual_vars with \code{residuals = TRUE}, extra columns to
#'   write alongside \code{fn}'s value, one per episode -- \code{frml}'s
#'   own \code{v=}, the same mechanism as \code{\link{tda_fml}}'s
#'   \code{residual_vars}. Only reliably supports a covariate from
#'   \code{formula}; TDA's episode built-ins (\code{ts}, \code{tf},
#'   \code{org}, \code{des}, \code{time}) cannot be named here, despite
#'   being available inside
#'   \code{fn}: they are expression-only aliases, not real data-matrix
#'   variables -- \code{v = ts,tf} is a
#'   TDA syntax error (their real, differently-named
#'   underlying variables, e.g. \code{Tstart} for \code{ts}, do work,
#'   but are internal, auto-generated by \code{edef()} and not
#'   something to build an interface around safely).
#' @param protocol ask TDA to also write its iteration-by-iteration
#'   diagnostic log -- \code{frml}'s \code{prot=}, the same
#'   mechanism as \code{\link{tda_fml}}'s \code{protocol}.
#' @param control optional \code{\link{tda_control}}.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return An object of class \code{tda_fit}. With \code{residuals = TRUE},
#'   also carries \code{$residuals}.
#' @family sequence analysis
#' @examples
#' set.seed(1)
#' n <- 200
#' d <- data.frame(x = rnorm(n))
#' haz <- exp(0.5 + 0.8 * d$x)
#' d$dur <- rexp(n, haz)
#' d$status <- 1
#' # the exponential hazard model, written out by hand -- matches
#' # tda_rate(Surv(dur, status) ~ x, d, model = "exponential"). des and
#' # (tf - ts) are TDA's built-ins, not d's dur/status columns -- see
#' # Details.
#' fit <- tda_frml(Surv(dur, status) ~ x, d,
#'                 c("xb = b0 + x * b1", "gamma = exp(xb)",
#'                   "fn = des * xb - gamma * (tf - ts)"))
#' coef(fit)
#' coef(tda_rate(Surv(dur, status) ~ x, d, model = "exponential"))
#'
#' # the same model, written in R syntax instead
#' fit2 <- tda_frml(Surv(dur, status) ~ x, d, {
#'     xb = b0 + x * b1
#'     gamma = exp(xb)
#'     fn = des * xb - gamma * (tf - ts)
#' })
#' coef(fit2)
#' @export
tda_frml <- function(formula, data, definitions, start = NULL,
                     id = NULL, spell = NULL,
                     control = NULL, constraints = NULL, residuals = FALSE,
                     residual_vars = NULL, protocol = FALSE,
                     options = list(), dir = tempfile("tda"),
                     ...) {
    definitions <- .fml_check_definitions(definitions, substitute(definitions),
                                          "the log-likelihood contribution of one episode")

    # strata() would be silently DROPPED here, not used: a run with a
    # strata() term completes and returns the identical coefficient as
    # the same call without it, which is the worst of the options --
    # refuse it instead so a stratified analysis is never quietly
    # unstratified.
    if (any(grepl("strata\\(", attr(stats::terms(formula, data = data),
                                     "term.labels"))))
        stop("tda_frml does not stratify; fit per stratum instead ",
             "(fml has no strata mechanism -- a strata() term would be ",
             "silently ignored)", call. = FALSE)
    d <- .tda_design(formula, data, id = id, spell = spell)
    # Failing loudly here instead of silently, since a plausible-looking
    # but silently-unstratified fit is a worse outcome than an error.
    if (!is.null(d$sname))
        stop("tda_frml does not support strata() -- it is silently ",
             "ignored by TDA's frml() (checked: a ",
             "stratified and an unstratified call give an identical ",
             "fit), not merely unimplemented here")
    # The definitions are written against the data frame's column
    # names, so they have to be rewritten to the names TDA was given --
    # xlab/xname is exactly that mapping, already built by .tda_design(),
    # but it only covers the formula's covariates. Surv()'s
    # arguments need the identical treatment: checked, a real
    # run referencing Surv(TFP, DES)'s TFP/DES by name inside
    # definitions -- exactly the way TDA's manual example
    # (examples/exam/frml1.cf) references its DES/DUR columns
    # directly -- was a TDA syntax error ("probably wrong
    # reference to a variable") without this, since those columns are
    # actually named Tfin/Des internally (.tda_design()'s naming),
    # and nothing translated the user's Surv() argument text to
    # match. Only covers a simple column reference (Surv(TFP, DES)), not
    # a computed one (Surv(TFin - TStart + 1, DES)) -- there is no
    # single column name to translate to for the latter, so it is left
    # alone; TDA's lowercase built-ins (ts/tf/org/des) still work
    # for that case, as they always did.
    surv_args <- as.list(formula[[2L]])[-1L]
    surv_names <- switch(as.character(length(surv_args)),
                        "2" = c("Tfin", "Des"),
                        "3" = c("Tstart", "Tfin", "Des"),
                        "4" = c("Tstart", "Tfin", "Org", "Des"))
    for (i in seq_along(surv_args))
        if (is.symbol(surv_args[[i]])) {
            lab <- as.character(surv_args[[i]])
            if (lab != surv_names[i])
                definitions <- gsub(sprintf("\\b%s\\b", lab), surv_names[i],
                                    definitions)
        }
    for (i in seq_along(d$xlab))
        if (d$xlab[i] != d$xname[i])
            definitions <- gsub(sprintf("\\b%s\\b", d$xlab[i]), d$xname[i],
                                definitions)
    # residual_vars only reliably supports covariates from the formula
    # -- checked: that TDA's episode
    # built-ins (ts, tf, org, des, time) are expression-only aliases,
    # not real v=-usable variable names (v = ts,tf is a TDA
    # syntax error) -- so an unrecognised name is rejected here, with a
    # clear reason, rather than silently forwarded into one.
    rv_tda <- if (!is.null(residual_vars)) {
        miss <- setdiff(residual_vars, d$xlab)
        if (length(miss))
            stop("`residual_vars` only supports covariates from `formula` ",
                 "-- not TDA's episode built-ins (ts, tf, org, des, ",
                 "time), which are not real v=-usable variable names: ",
                 paste(miss, collapse = ", "))
        d$xname[match(residual_vars, d$xlab)]
    }

    opts <- .fml_opts(start, constraints, control, options)
    opts <- .fml_residuals_opts(opts, residuals, rv_tda)
    opts <- .fml_protocol_opts(opts, protocol)
    cmd <- do.call(tda_block, c(list(name = "frml"), opts,
                                list(rhs = paste(definitions,
                                                 collapse = ",\n    "))))
    res <- .tda_fit(d, cmd, dir, ...)
    .tda_result(res, list(n = d$n, xlab = d$xlab, xname = d$xname),
               match.call(), "tda_frml",
               list(residuals = .fml_residuals_read(res, residuals, residual_vars),
                    protocol = .fml_protocol_read(res, protocol)),
               data = data)
}

#' Evaluate an inclusion function over an interval
#'
#' \code{evalfi}: given an expression and an interval domain for each of
#' its arguments, returns an interval guaranteed to contain the
#' expression's true range over that domain -- interval arithmetic, not a
#' numeric evaluation at a point. The bound is conservative rather than
#' tight where a variable appears more than once (the standard
#' "dependency problem" of interval arithmetic: each occurrence is
#' widened independently, as if they were unrelated variables that
#' happened to share a domain) -- \code{evalfi(x = c(0, 3), "x*x-2*x+1")}
#' returns \code{[-5, 10]}, not the true range \code{[0, 4]} that
#' \code{(x-1)^2} takes over \code{[0, 3]}. \code{evalfi1} is the same but
#' for the first derivative.
#'
#' \code{^} is not a valid operator in an interval expression -- write
#' a power as repeated multiplication (\code{x*x}, not \code{x^2}).
#'
#' @param expr the expression, as a string, in TDA's language.
#' @param ... one interval per argument the expression uses, named for
#'   the argument, each as \code{c(lower, upper)} --
#'   \code{evalfi(x = c(0, 3), expr = "...")}.
#' @param derivative if \code{TRUE}, use \code{evalfi1} (the first
#'   derivative's inclusion function) instead of \code{evalfi}.
#' @param options a named list of further TDA options, passed through.
#' @param dir working directory.
#' @return A named numeric vector, \code{c(lower, upper)}.
#' @family sequence analysis
#' @examples
#' # true range of (x-1)^2 over [0,3] is [0,4] -- this is wider, since x
#' # appears twice and interval arithmetic cannot see they are the same x
#' tda_evalfi("x*x-2*x+1", x = c(0, 3))
#' tda_evalfi("x*x-2*x+1", x = c(0, 3), derivative = TRUE)  # 2x-2, exact
#'                                                          # here: [-2,4]
#' @export
tda_evalfi <- function(expr, ..., derivative = FALSE, options = list(),
                       dir = tempfile("tda")) {
    args <- list(...)
    if (!length(args) || is.null(names(args)) || any(!nzchar(names(args))))
        stop("give a named interval for every argument, e.g. x = c(0, 3)")
    opts <- stats::setNames(lapply(args, function(a) paste(a, collapse = ",")),
                            names(args))
    opts <- c(opts, list(fmt = "18.12"), .tda_extra(options))
    cmd <- if (isTRUE(derivative)) "evalfi1" else "evalfi"
    res <- tda_run(do.call(tda_block, c(list(name = cmd), opts, list(rhs = expr))),
                   dir = dir)
    err <- grep("^Error|Can't evaluate", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not evaluate this: ", err[1L], call. = FALSE)
    # The echoed argument domain ("Arguments: x=[0,3]") is bracketed the
    # same way a result is, and comes first -- matching on "function
    # value"'s/"gradient"'s line, not just the first bracket found, is
    # what tells them apart. evalfi1 prints both the function value and
    # the gradient; derivative = TRUE wants the second, not the first.
    # export first: the interval is exported at the print site, so no
    # bracket-matching against the console is needed.  evalfi1 prints
    # the function value and the gradient separately and each has its
    # own key, which is what told them apart in the parse too.
    ekey <- if (isTRUE(derivative)) "evalfi.gradient" else "evalfi.value"
    em <- if (.use_exports()) res$exports[[ekey]]
    if (is.matrix(em) && ncol(em) == 2L && nrow(em) >= 1L)
        return(stats::setNames(as.vector(em[1L, ]), c("lower", "upper")))
    pat <- if (isTRUE(derivative)) "gradient" else "function value"
    i <- grep(pat, res$output)
    if (!length(i))
        stop("evalfi produced no readable result")
    m <- regmatches(res$output[i[1L]],
                    regexpr("\\[\\s*[0-9.eE+-]+\\s*,\\s*[0-9.eE+-]+\\s*\\]",
                            res$output[i[1L]]))
    if (!length(m) || !nzchar(m))
        stop("evalfi produced no readable result")
    bounds <- as.numeric(strsplit(gsub("[][ ]", "", m), ",")[[1]])
    stats::setNames(bounds, c("lower", "upper"))
}
