# Assembling R objects from the direct exports (CONTRIBUTING.md).
#
# Two rules, both from review:
#
# 1. Never coerce mixed content into a matrix.  A numeric export plus
#    its label vector becomes a DATA FRAME with one typed column each --
#    cbind()ing a character vector onto a numeric matrix silently turns
#    everything into a character matrix, which is exactly the accident
#    the exports exist to avoid.
#
# 2. Two printed tables stay two objects.  Binding across blocks is
#    only defined when the blocks conform -- same number of columns,
#    same column names in the same order.  TDA can legitimately print
#    non-conforming siblings (a life table's per-destination columns
#    follow each block's destination set), and an unconditional
#    do.call(rbind, ...) would glue columns that do not belong
#    together.  .exports_bind() refuses that loudly instead.

# Join a numeric export with its label export into a proper data frame.
# labels_key defaults to "<name>.names"; a missing label export gives a
# frame of the numeric columns alone.
.exports_frame <- function(exports, name,
                           labels_key = paste0(name, ".names"),
                           label_col = "Variable",
                           col_names = NULL) {
    m <- exports[[name]]
    if (is.null(m))
        return(NULL)
    d <- as.data.frame(m)
    if (!is.null(col_names) && length(col_names) == ncol(d))
        names(d) <- col_names
    lb <- exports[[labels_key]]
    if (!is.null(lb) && length(lb) == nrow(d)) {
        d <- cbind(stats::setNames(data.frame(as.character(lb),
                                              stringsAsFactors = FALSE),
                                   label_col),
                   d)
    }
    d
}

# rbind a list of export matrices or data frames, but only when they
# conform; otherwise stop with a message that names the offender.  The
# deliberate contrast with do.call(rbind, ...): shape agreement is
# checked.
.exports_bind <- function(blocks) {
    if (!length(blocks))
        return(NULL)
    if (length(blocks) == 1L)
        return(blocks[[1L]])
    ref <- blocks[[1L]]
    ref_nc <- ncol(ref)
    ref_nm <- colnames(ref)
    for (k in seq_along(blocks)[-1L]) {
        b <- blocks[[k]]
        if (ncol(b) != ref_nc ||
            !identical(colnames(b), ref_nm))
            stop("blocks do not conform (block ", k, " has ",
                 ncol(b), " columns", if (!is.null(colnames(b)))
                     paste0(": ", paste(colnames(b), collapse = ", ")),
                 "; block 1 has ", ref_nc,
                 if (!is.null(ref_nm))
                     paste0(": ", paste(ref_nm, collapse = ", ")),
                 ") -- they are separate tables and stay separate",
                 call. = FALSE)
    }
    do.call(rbind, blocks)
}

# The master switch: exports are preferred wherever a reader has been
# switched, and options(tdaR.use_exports = FALSE) forces every such
# reader back onto the text/file parser, which remains as the fallback.
.use_exports <- function() {
    isTRUE(getOption("tdaR.use_exports", TRUE))
}

# The full-precision square matrix for cov/corr/rcorr fits: the export,
# dimnamed with the caller's labels, or NULL so the caller falls back to
# the text parser (flag off, absent export, or a shape mismatch).
.exports_sqmat <- function(res, key, d) {
    if (!.use_exports())
        return(NULL)
    m <- res$exports[[key]]
    if (!is.matrix(m) || nrow(m) != ncol(m) || nrow(m) != length(d$xlab))
        return(NULL)
    dimnames(m) <- list(d$xlab, d$xlab)
    m
}

# Overlay the export's full-precision numeric columns onto a parsed
# estimates table when shapes line up.  key: coeff / rate.est /
# qreg.est.  The export's NaN (the printed ---) becomes NA to match the
# parser.  Returns the table unchanged when the flag is off, the export
# is absent, or the rows do not align.
.est_overlay <- function(est, res, key, which = NULL) {
    if (!.use_exports() || !is.data.frame(est))
        return(est)
    key <- key[key %in% names(res$exports)][1L]
    if (is.na(key))
        return(est)
    # A command can print several coefficient tables in one run (zreg
    # prints two), so the export arrives as key, key.2, ...  The parser
    # keeps the LAST table it finds, so the last conforming block is the
    # matching one -- verified against zreg, where the first block holds
    # a different fit entirely and overlaying it would have replaced the
    # estimates with the wrong numbers rather than merely failing.
    blocks <- .exports_blocks(res$exports, key)
    blocks <- Filter(function(b) is.matrix(b) && nrow(b) == nrow(est),
                     blocks)
    if (!length(blocks))
        return(est)
    # `which` names the printed table this frame came from, for callers
    # that pick a block themselves (.reg_run's block = "first"/"last");
    # without it the last conforming block is used.
    i <- if (is.null(which)) length(blocks)
         else max(1L, min(as.integer(which), length(blocks)))
    m <- blocks[[i]]
    v <- if (key == "coeff") m else m[, -1L, drop = FALSE]
    if (ncol(v) != 4L)
        return(est)
    v[is.nan(v)] <- NA_real_
    for (j in seq_along(cn <- c("Coeff", "Error", "C/Error", "Signif")))
        if (cn[j] %in% names(est))
            est[[cn[j]]] <- v[, j]
    est
}

# Put an exported column back into a parsed frame without changing the
# column's type.  Every export is a double, so assigning one
# straight onto a column the parser had read as integer silently
# widens it -- which broke identical() against an integer state vector
# and would change the type of every count column in the sequence
# tables.  An integer column stays integer whenever the exported
# values really are whole numbers in range; anything else becomes the
# double it is.
.overlay_col <- function(old, v) {
    v[is.nan(v)] <- NA_real_
    if (is.integer(old)) {
        fin <- v[!is.na(v)]
        if (!length(fin) ||
            (all(fin == trunc(fin)) && all(abs(fin) <= .Machine$integer.max)))
            return(as.integer(v))
    }
    v
}

# Replace a parsed data frame's numeric content with the export's
# full-precision values when shapes agree (column names and any
# non-numeric columns stay the parser's). NaN in the export becomes NA.
.overlay_num <- function(df, m) {
    if (!.use_exports() || !is.data.frame(df) || !is.matrix(m))
        return(df)
    if (nrow(m) != nrow(df) || ncol(m) != ncol(df))
        return(df)
    for (j in seq_len(ncol(df)))
        if (is.numeric(df[[j]]))
            df[[j]] <- .overlay_col(df[[j]], m[, j])
    df
}

# Overlay named columns of a parsed frame from an export whose columns
# are a subset, in order: the parser's non-numeric columns (a
# variable name, a class label) have no counterpart in the export, so
# the column counts .overlay_num insists on never line up.  cols names
# the parsed columns the export's columns correspond to, left to right.
.overlay_cols <- function(df, m, cols) {
    if (!.use_exports() || !is.data.frame(df) || !is.matrix(m))
        return(df)
    if (nrow(m) != nrow(df) || ncol(m) != length(cols) ||
        !all(cols %in% names(df)))
        return(df)
    for (j in seq_along(cols))
        df[[cols[j]]] <- .overlay_col(df[[cols[j]]], m[, j])
    df
}

# All blocks a repeated-table producer wrote, in production order:
# "key", "key.2", "key.3", ... .  Returns them stacked when they
# conform, which is how a reader that parsed one concatenated frame
# out of several printed blocks gets its counterpart.
.exports_blocks <- function(exports, key) {
    ks <- c(key, grep(paste0("^", key, "\\.[0-9]+$"), names(exports),
                      value = TRUE))
    ks <- ks[ks %in% names(exports)]
    if (!length(ks))
        return(list())
    ks <- ks[order(c(1L, as.integer(sub(".*\\.", "", ks[-1L]))))]
    unname(exports[ks])
}

.exports_stack <- function(exports, key) {
    b <- .exports_blocks(exports, key)
    if (!length(b))
        return(NULL)
    .exports_bind(b)
}

# The generic console tap (print.values) stages every number printf1()
# formats as (line, value).  This pulls the values off one output line,
# which makes any scalar TDA prints reachable at full precision without
# a producer of its own -- the fallback for the long tail of "Label:
# value" lines that have no table behind them.
# Find a tapped value by the LABEL TDA printed in front of it, rather
# than by grepping res$output for its line -- the only form of the
# lookup that works with no printed text at all.
#
# The stored label is the literal text of the format string in front of
# the conversion, which can be only the TAIL of what the reader sees:
# TDA sometimes splits one printed line across two printf1() calls, so
# "Rank of least squares data matrix:" is stored as "least squares data
# matrix:".  Matching is therefore on the tail, after dropping the
# trailing colon from both sides.
.tap_by_label <- function(res, label, which = 1L) {
    if (!.use_exports())
        return(NULL)
    pv <- res$exports[["print.values"]]
    lb <- res$exports[["print.labels"]]
    if (!is.matrix(pv) || !is.character(lb) || nrow(pv) != length(lb))
        return(NULL)
    norm <- function(x) sub("[[:space:]:]+$", "", trimws(x))
    want <- norm(label)
    have <- norm(lb)
    i <- which(nzchar(have) &
               (have == want | endsWith(want, have) | startsWith(have, want)))
    if (!length(i))
        return(NULL)
    v <- pv[i[1L], 2L]
    if (length(v)) v else NULL
}

# `which` picks one of the values on that line, in print order -- a
# line often prints more than one number ("Best minimal function value:
# X   best lower bound: Y"), so requiring exactly one silently fell back
# to the text for every such line.
# The doubles behind output line k, from the console tap; NULL when the
# tap is off or has nothing for that line.  For a parser that has found
# its line by text and split it into numbers, .tap_line() gives the same
# numbers unrounded when the counts agree -- the text stays in charge of
# WHICH line and how many fields, the tap supplies the values.
.tap_line <- function(res, k) {
    if (!.use_exports())
        return(NULL)
    pv <- res$exports[["print.values"]]
    if (!is.matrix(pv))
        return(NULL)
    v <- pv[pv[, 1L] == k, 2L]
    if (length(v)) v else NULL
}

.tap_or <- function(res, k, parsed) {
    v <- .tap_line(res, k)
    if (!is.null(v) && length(v) == length(parsed) &&
        all(is.na(parsed) | abs(v - parsed) <= 5.1e-5 * pmax(1, abs(v))))
        v
    else parsed
}

# .tap_values reads the tap behind the FIRST line matching `pattern`.
# Where a command prints the same line once per case -- subm's Distance:,
# one per comparison -- the values wanted are one per matching line, so
# this walks every match instead and returns them in printed order.
.tap_values_each <- function(res, pattern, which = 1L) {
    if (!.use_exports())
        return(NULL)
    pv <- res$exports[["print.values"]]
    if (!is.matrix(pv) || is.null(res$output))
        return(NULL)
    i <- grep(pattern, res$output)
    if (!length(i))
        return(NULL)
    out <- vapply(i, function(k) {
        v <- pv[pv[, 1L] == k, 2L]
        if (length(v) >= which) v[which] else NA_real_
    }, numeric(1L))
    if (anyNA(out)) NULL else out
}

.tap_values <- function(res, pattern, which = NULL) {
    if (!.use_exports())
        return(NULL)
    pv <- res$exports[["print.values"]]
    if (!is.matrix(pv) || is.null(res$output))
        return(NULL)
    i <- grep(pattern, res$output)
    if (!length(i))
        return(NULL)
    v <- pv[pv[, 1L] == i[1L], 2L]
    if (!length(v))
        return(NULL)
    if (!is.null(which))
        return(if (which <= length(v)) v[which] else NULL)
    v
}

# A frame built from an export, typed the way read.table() would have
# typed the file it replaces: a column whose values are all whole
# numbers in range comes back integer, not double.
#
# Without this the two paths are not interchangeable -- switching a
# reader to the export silently changed integer columns to numeric, and
# identical() against them started failing in callers that had nothing
# to do with the switch.
.export_frame <- function(m) {
    d <- as.data.frame(m)
    for (j in seq_along(d)) {
        v <- d[[j]]
        if (!is.numeric(v))
            next
        fin <- v[is.finite(v)]
        if (length(fin) && all(fin == trunc(fin)) &&
            all(abs(fin) <= .Machine$integer.max))
            d[[j]] <- as.integer(v)
    }
    d
}

# The counterpart to .overlay_num(): BUILD a frame from an
# export instead of replacing the numeric content of a parsed one.
# `cols` names the export's columns, left to right; a reader that can
# call this no longer needs the printed text at all, which is the
# property test-phase3.R checks by emptying res$output.
#
# Returns NULL -- sending the caller back to its parser -- when the flag
# is off, the export is missing, or the column count does not match, so
# a producer whose shape drifts degrades to the old path rather than
# building a mislabelled table.
.frame_from_export <- function(res, key, cols, row_names = NULL) {
    if (!.use_exports())
        return(NULL)
    m <- res$exports[[key]]
    if (!is.matrix(m) || ncol(m) != length(cols))
        return(NULL)
    m[is.nan(m)] <- NA_real_
    d <- as.data.frame(m)
    names(d) <- cols
    if (!is.null(row_names) && length(row_names) == nrow(d))
        rownames(d) <- row_names
    d
}

# For the estimate tables: BUILD the frame from the numeric
# export and its label export, rather than overlaying onto the parsed
# console table.
#
# Each family's label export carries the row's whole identity, "|"
# separated, which is what makes this possible without the header line:
#   coeff.names     "<wave>|<variable>"          (wave "-" for the intercept)
#   rate.est.names  "<sn>|<org>|<des>|<MT>|<name>"  (plus P<t> for MEXP2)
#   qreg.est.names  "<cat>|<term>|<variable>"
# The numeric export is the four printed columns, with rate/qreg
# carrying their Idx in front.
#
# Returns NULL when anything does not line up, sending the caller back
# to the parser -- a label export whose field count drifts must not
# produce a confidently mislabelled table.
#
# The index/identifier columns are numeric, not integer, because that is
# what the parser produces: switching the type would be a user-visible
# change with nothing to gain, and identical() between the two paths is
# what the phase-4 tests check.
.est_one_block <- function(key, m, lb, unmangle) {
    if (!is.matrix(m) || !is.character(lb) || nrow(m) != length(lb))
        return(NULL)
    parts <- strsplit(lb, "|", fixed = TRUE)
    nf <- unique(lengths(parts))
    if (length(nf) != 1L)
        return(NULL)
    # coeff and fml.est are the whole four columns; rate.est and qreg.est
    # carry their Idx in front of them
    v <- if (key == "coeff" || key == "fml.est") m
         else m[, -1L, drop = FALSE]
    if (ncol(v) != 4L)
        return(NULL)
    v[is.nan(v)] <- NA_real_
    fld <- function(i) vapply(parts, `[`, "", i)
    # TDA prints "-" where a column has no value for this row -- the
    # ordered models have no category, the intercept has no wave -- and
    # as.numeric() on it is an NA plus a "NAs introduced by coercion"
    # warning on every single fit.  The NA is right; the warning is
    # noise, so the conversion is explicit about what "-" means instead
    # of letting the coercion complain.
    num <- function(i) {
        x <- trimws(fld(i))
        x[x == "-" | x == ""] <- NA_character_
        suppressWarnings(as.numeric(x))
    }
    head_cols <- switch(
        key,
        coeff = list(Idx = as.numeric(seq_along(lb)),
                     Wave = num(1L),
                     Variable = unmangle(fld(2L))),
        rate.est = if (nf == 5L)
            list(Idx = as.numeric(m[, 1L]), SN = num(1L),
                 Org = num(2L), Des = num(3L),
                 MT = fld(4L), Variable = unmangle(fld(5L)))
        else NULL,
        # fml/freg/frml print "Parameter", "Value" and "Value/E"; by the
        # time tda_estimates() is done the parsed frame calls them
        # Parameter, Coeff and C/Error, and the overlay only applies when
        # the two frames agree on names -- so these match that, not the
        # raw header.
        fml.est = list(Idx = as.numeric(seq_along(lb)),
                       Parameter = unmangle(lb)),
        qreg.est = if (nf == 3L)
            list(Idx = as.numeric(m[, 1L]), Cat = num(1L),
                 Term = fld(2L), Variable = unmangle(fld(3L)))
        else NULL,
        NULL)
    if (is.null(head_cols))
        return(NULL)
    est <- as.data.frame(head_cols, stringsAsFactors = FALSE,
                         check.names = FALSE)
    est[c("Coeff", "Error", "C/Error", "Signif")] <- as.data.frame(v)
    est
}

.est_from_export <- function(res, key, unmangle = identity) {
    if (!.use_exports())
        return(NULL)
    key <- key[key %in% names(res$exports)][1L]
    if (is.na(key))
        return(NULL)
    # A command can print several estimate tables in one run (zreg does),
    # so the exports arrive as key, key.2, ... with their labels
    # alongside.  One frame is built per block and the LIST is returned
    # when there is more than one, exactly as tda_estimates() does --
    # returning only the first silently handed back the wrong fit's
    # numbers, which is how this was caught.
    mats <- .exports_blocks(res$exports, key)
    labs <- .exports_blocks(res$exports, paste0(key, ".names"))
    if (length(mats) > 1L) {
        if (length(labs) != length(mats))
            return(NULL)
        out <- lapply(seq_along(mats), function(i)
            .est_one_block(key, mats[[i]], labs[[i]], unmangle))
        if (any(vapply(out, is.null, NA)))
            return(NULL)
        return(out)
    }
    .est_one_block(key, res$exports[[key]],
                   res$exports[[paste0(key, ".names")]], unmangle)
}

# The label/value pairs prn_sfmt prints (contingency measures, glm's
# deviance block, and the scattered scalars beside them) come back as
# two parallel exports; this turns them into one named numeric vector,
# in print order, so a reader can pick by the label TDA itself used.
.exports_sfmt <- function(res) {
    if (!.use_exports())
        return(NULL)
    v <- res$exports[["sfmt.values"]]
    l <- res$exports[["sfmt.labels"]]
    if (!is.matrix(v) || !is.character(l) || nrow(v) != length(l))
        return(NULL)
    stats::setNames(as.vector(v), trimws(l))
}

# Overlay a list of parsed blocks with suffixed exports (key, key.2,
# ...), in order; blocks whose shape does not match stay parsed.
.overlay_blocks <- function(blocks, exports, key) {
    if (!.use_exports() || !length(blocks))
        return(blocks)
    ks <- c(key, paste0(key, ".", seq_len(max(0L, length(blocks) - 1L)) + 1L))
    for (k in seq_along(blocks))
        if (!is.null(exports[[ks[k]]]))
            blocks[[k]] <- .overlay_num(blocks[[k]], exports[[ks[k]]])
    blocks
}
