#!/usr/bin/env Rscript
# Argument-name consistency across the exported API.
#
# Two things go wrong when wrappers are written one at a time:
#
#  1. ONE NAME, TWO MEANINGS.  `opt` used to be a formal on tda_dma
#     (preprocessing) and tda_pdatd (the distance measure) -- TDA's own
#     spelling carried into R, where it says nothing. Both are renamed
#     (preprocess, measure). This check lists names shared by several
#     functions so the meanings can be compared; `select` (a case
#     selection expression, TDA's sel=) is the model to follow.
#
#  2. PARTIAL-MATCHING TRAPS.  R matches a supplied name against formals
#     BEFORE `...`, so a formal named `options` swallows a call written
#     `opt = 2` -- it never reaches the dots and `options` stops being a
#     list. Any formal that is a strict prefix of another formal in the
#     same function, or of `options`, is flagged.
#
#     Rscript tools/check_api_names.R
suppressMessages(library(tdaR))
ns <- asNamespace("tdaR")
ex <- sort(grep("^tda_", getNamespaceExports("tdaR"), value = TRUE))
fl <- lapply(ex, function(f) names(formals(get(f, envir = ns))))
names(fl) <- ex

# 1. shared names, for eyeballing meaning
all_nm <- table(unlist(lapply(fl, function(x) setdiff(x, "..."))))
shared <- sort(all_nm[all_nm > 1L], decreasing = TRUE)
cat("argument names used by more than one function:", length(shared), "\n")

# 2. partial-matching traps.
#
# The trap is NOT two formals where one prefixes the other -- an exact
# name always wins over a partial match, so tda_seqmd(event = x) binds
# to `event`, not `event_covariates`.  My first version of this check
# flagged 24 of those and every one was a false positive.
#
# The real trap is a formal BEFORE `...` that prefixes a name the user
# would reasonably type but which is NOT a formal: a TDA option going to
# the command through the dots.  `options` swallowing `opt` is the case
# that bit.
hlp <- readLines(system.file("extdata", "tda.hlp", package = "tdaR"),
                 warn = FALSE)
heads <- grep("^##", hlp)
hn <- trimws(sub("^##", "", hlp[heads]))
cmd_opts <- function(cmd) {
    i <- match(cmd, hn)
    if (is.na(i)) return(character(0))
    to <- if (i < length(heads)) heads[i + 1L] - 1L else length(hlp)
    sec <- hlp[seq.int(heads[i] + 1L, to)]
    m <- regmatches(sec, regexpr("^\\s{2,}([a-zA-Z][a-zA-Z0-9_]*)\\s*=", sec))
    unique(trimws(sub("=$", "", trimws(m))))
}
bad <- character(0)
for (f in ex) {
    fo <- fl[[f]]
    if (!("..." %in% fo)) next
    before <- fo[seq_len(match("...", fo) - 1L)]
    b <- paste(deparse(get(f, envir = ns)), collapse = " ")
    m <- regmatches(b, regexpr('name = "[a-z0-9_]+"', b))
    if (!length(m)) next
    opts <- cmd_opts(sub('name = "', "", sub('"$', "", m[1L])))
    opts <- setdiff(opts, fo)          # a formal of its own is safe
    # DIRECTION MATTERS, and I got it wrong twice.  R matches a SUPPLIED
    # name that is a PREFIX OF a formal: supplying opt= where the formal
    # is `options` binds to options.  The reverse is safe -- supplying
    # gt= where the formal is `g` goes to the dots, because gt is longer
    # than g.  Checking startsWith(option, formal) flagged `g` swallowing
    # gt= on three graph wrappers and `s` swallowing sc= on tda_map, and
    # all four were false positives.
    # Only `options` counts as harmful.  Most swallows are the API
    # working as intended: v= is swallowed by `variables`, sel= by
    # `select`, k= by `kernel` -- in each case the formal IS the
    # descriptive name for that option and is what the user should
    # write.  `options` is different: it is a generic container, so a
    # setting it swallows is simply lost, and the user gets "options
    # must be a named list" or, worse, silence.
    # Wrappers where the setting already HAS a descriptive formal, so
    # being unable to type TDA's own spelling is the intended outcome
    # rather than a loss:
    #   tda_dma      opt -> preprocess
    #   tda_pdatd    opt -> measure
    #   tda_g, tda_dmet, tda_ptree
    #                opt is gdd's input format, and is set upstream by
    #                tda_graph(form = "edges"/"dissimilarity") -- the
    #                user never passes it to these at all
    #   tda_npreg    opt -> method ("mean"/"quantile"/"frequency").  TDA
    #                documents 4=lowess and 5=midmeans too, but opt=4 is
    #                a genuine syntax error in the shipped source --
    #                established earlier, see the comment in tda_npreg
    #                itself -- so method covers everything that works.
    #   tda_sddf     opt -> self_consistent (opt=2)
    #   tda_qreg     opt is "model-specific options"; the models that use
    #                it have their own formals (parameterization, nq,
    #                waves, ...)
    #   tda_locate_line
    #                the first block in the body is gdd, whose opt=1 is
    #                the edge-list input format, set internally; the
    #                gloc block that carries `options` has no opt=
    covered <- c("tda_dma", "tda_pdatd", "tda_g", "tda_dmet", "tda_ptree",
                 "tda_npreg", "tda_sddf", "tda_qreg", "tda_locate_line")
    if (f %in% covered) next
    for (q in opts)
        if (startsWith("options", q) && q != "options" &&
            "options" %in% before)
            bad <- c(bad, sprintf("%s: %s= is swallowed by `options`", f, q))
}
cat("partial-matching traps:", length(bad), "\n")
if (length(bad)) {
    cat(paste(unique(bad), collapse = "\n"), "\n")
    quit(status = 1)
}
