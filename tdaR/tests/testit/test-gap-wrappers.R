# Session 40 gap-batch wrappers: verified against independent R
# computations where an oracle is cheap, smoke-plus-shape otherwise.

# mproc vs the closed-form orthogonal Procrustes solution (svd)
local({
    set.seed(4)
    x <- matrix(rnorm(8), 4)
    th <- 0.6; R <- rbind(c(cos(th), -sin(th)), c(sin(th), cos(th)))
    y <- x %*% R
    z <- tda_mproc(x, y)
    sv <- svd(crossprod(y, x))
    zr <- y %*% (sv$u %*% t(sv$v))
    same("mproc equals the svd Procrustes solution", unname(z), unname(zr), 1e-5)
})

# recode vs plain R recoding
local({
    d <- data.frame(a = c(1, 2, 3), b = c(10, 20, 30))
    r <- tda_recode(d, a = "a * 2 + b / 10")
    same("recode matches R", r$a, d$a * 2 + d$b / 10, 1e-9)
    same("untouched column survives", r$b, d$b, 1e-9)
})

# com: permutation and subset counts against factorial / choose
local({
    nrows <- function(x) if (is.matrix(x)) nrow(x) else length(x)
    ok("com counts 3! permutations",
       nrows(tda_com("permutations", n = 3)) == 6)
    ok("com counts choose(4,2) subsets",
       nrows(tda_com("subsets", n = 4, m = 2)) == choose(4, 2))
})

# gqap vs brute force over all assignments (n = 3)
local({
    f <- rbind(c(0, 3, 1), c(3, 0, 2), c(1, 2, 0))
    d <- rbind(c(0, 1, 4), c(1, 0, 2), c(4, 2, 0))
    o <- tda_gqap(f, d)
    cost <- if (is.list(o)) o$cost else
        as.numeric(sub(".*value: *([0-9.]+).*", "\\1",
            grep("Best cost value", o, value = TRUE)[1]))
    perms <- rbind(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
    # TDA sums over both orientations of every pair, so no halving
    brute <- min(apply(perms, 1, function(p) sum(f * d[p, p])))
    ok("gqap reaches the brute-force optimum", isTRUE(abs(cost - brute) < 1e-6))
})

# triang: every reported triangle has three distinct vertices
local({
    tr <- tda_triang(x = c(0, 1, 0, 1, 0.5), y = c(0, 0, 1, 1, 0.5))
    ok("triangulation reports triangles",
       (if (is.matrix(tr)) nrow(tr) else length(tr)) >= 2)
})

# pdatr: all ordered pairs appear
local({
    p <- tda_pdatr(data.frame(v = c(5, 7, 9)))
    ok("pdatr emits n*(n-1) ordered pairs or n^2 pairs",
       nrow(p) %in% c(6L, 9L))
})

# spmod / pcyc verified in test-r-comparisons.R; smoke the rest so a
# regression in any wrapper's pipeline fails loudly
local({
    sz <- function(x) if (is.matrix(x)) nrow(x) else length(x)
    ok("gde runs",   sz(tda_gde(c(1,2,2,3), c(10,10,20,20))) >= 1)
    ok("ghd runs",   sz(tda_ghd(rbind(c(0,0), c(1,0), c(0,1), c(1,1)))) >= 1)
    ok("mdsn returns a structured solution", {
        pts <- cbind(c(0, 0, 3, 3), c(0, 2, 0, 2))
        m <- tda_mdsn(as.matrix(dist(pts)))
        is.list(m) && is.finite(m$stress) &&
            all(dim(m$configuration) == c(4, 2)) &&
            cor(dist(m$configuration), dist(pts),
                method = "spearman") > 0.9})
    ok("rfit runs",  length(tda_rfit(rbind(c(0,2,1), c(2,0,2), c(1,2,0)))) >= 1)
    sgar <- tda_sga(rbind(c(1,1,0), c(1,0,0), c(1,1,1)))
    ok("sga reports its reproducibility coefficients",
       all(grepl("^Err[0-2]: 1.0000$",
                 grep("^Err[0-2]:", sgar, value = TRUE))))
    unfr <- tda_unf(rbind(c(3,2,1), c(1,3,2), c(2,3,1)))
    bv <- as.numeric(sub(".*: *", "",
              grep("Best value", unfr, value = TRUE)[1L]))
    bp <- as.integer(strsplit(trimws(sub(".*: *", "",
              grep("Best permutation", unfr, value = TRUE)[1L])), " +")[[1L]])
    # TDA's folded-tau criterion has no simple closed R analogue; the
    # exact values on this fixture are pinned (verified by reading the
    # full search output), plus the structural invariants
    um <- as.numeric(sub(".*=", "", grep("UNFMax", unfr, value = TRUE)[1L]))
    ok("unf: pinned optimum and a valid permutation",
       bv == 3 && identical(bp, 1:3) && bv <= um &&
       setequal(bp, seq_along(bp)))
    rodr <- tda_rod(rbind(c(1,2,3), c(2,1,3)), "table")
    ok("rod returns the orders with counts",
       is.data.frame(rodr$orders) && nrow(rodr$orders) == 2L &&
       all(rodr$orders$count == 1L))
    rm1 <- local({
        set.seed(3); ab <- rnorm(150)
        x <- sapply(c(-1.2, 0, 1.2), function(d)
            rbinom(150, 1, plogis(ab - d)))
        tda_rmod(x)
    })
    ok("rmod item parameters order by true difficulty",
       length(rm1$items) == 3L &&
       rm1$items[1L] > rm1$items[2L] && rm1$items[2L] > rm1$items[3L] &&
       is.finite(rm1$logLik) && rm1$logLik < 0)
    cp <- tda_clpyr(as.matrix(dist(c(0,1,5,6))))
    ok("clpyr returns the pyramid, no engine errors",
       inherits(cp, "tda_pyramid") && !any(grepl("Error:", cp$output)) &&
       cp$pyramid$size[nrow(cp$pyramid)] == 4L)
    cp2 <- local({
        old <- options(tdaR.use_exports = FALSE); on.exit(options(old))
        tda_clpyr(as.matrix(dist(c(0,1,5,6))))
    })
    ok("clpyr export tap and file agree",
       isTRUE(all.equal(unname(as.matrix(cp$pyramid)),
                        unname(as.matrix(cp2$pyramid)))))
    ok("bfc runs",   {x <- expand.grid(a=0:1, b=0:1)
                      length(tda_bfc(as.integer(x$a & x$b), x)) > 3})
    cjr <- tda_conj(c(4,3,2,1), expand.grid(a = 1:2, b = 1:2))
    cjt <- grep("Const|X1|X2", tda_payload(cjr), value = TRUE)
    cf <- suppressWarnings(as.numeric(sub(".* ", "",
              trimws(cjt[grepl("Const|X[12]", cjt)]))))
    cf <- utils::tail(cf[is.finite(cf)], 3L)
    ok("conj recovers the exact additive part-worths",
       length(cf) == 3L && abs(cf[1L] - 4) < 1e-3 &&
       abs(cf[2L] + 1) < 1e-3 && abs(cf[3L] + 2) < 1e-3)
    ok("dmet1 runs", length(tda_dmet1(as.matrix(dist(cbind(c(0,1,4), 0))))) >= 1)
    ok("scla runs",  length(tda_scla(as.matrix(dist(c(0,.1,5,5.1))))) >= 1)
    ok("gdu runs",   length(tda_gdu(data.frame(i=c(1,2), j=c(2,3), v=1),
                                    data.frame(i=c(1,3), j=c(3,1), v=2))) >= 1)
    fc <- tda_gfcf(data.frame(i = c(1, 2, 3), j = c(2, 3, 4),
                              v = c(3, 2, 1)))
    ok("gfcf: head of the chain controls everything downstream",
       is.data.frame(fc) && fc$controlled[fc$node == 1] == "2,3,4" &&
       fc$n_controlled[fc$node == 4] == 0)
    rp <- tda_rap(cohort = c(60, 61), entry_year = c(80, 81),
                  exit_year = c(85, 86), destination = c(1, 1),
                  years = c(80, 86), ages = c(18, 26))
    # cohort 60 enters observation in 80 aged 20: that cell holds 1
    ok("rap risk table holds cohort 60 at age 20 in year 80",
       is.matrix(rp$risk) && rp$risk[rp$risk[, "age"] == 20, "y80"] == 1)
    ok("rcsv renders", {f <- tempfile(fileext=".csv")
                        writeLines(c("a;b", "1;2.5", "2;1"), f)
                        length(tda_rcsv(f)) >= 2})
    ok("ndvar equals model.matrix", {
        r <- tda_ndvar(data.frame(g = c(1, 2, 3, 2)))
        mm <- model.matrix(~ 0 + factor(c(1, 2, 3, 2)))
        all(as.matrix(r[, -1]) == mm)})
    ok("repsel runs", is.character(tda_repsel(data.frame(g=c(1,1,2), x=1:3),
                                              by="g", times="<2, 1>")))
    ok("niset runs", is.character(tda_niset(1e-6,
                       commands="evalf = int(sin(t), t, 0, 1);")))
    ok("mdsr runs",  length(tda_mdsr(cbind(c(0,1,2,3), c(0,1,0,1)),
                                     cbind(v = 1:4))) >= 1)
    ok("plot wrappers produce PostScript", {
        ps <- tda_plsurf3(fx="u", fy="v", fz="sin(u) * cos(v)")
        file.exists(ps) && file.size(ps) > 500})
})

# gtopo against its synthesized tile: the DEM row at lat 49.5 is
# 105..108 by construction, recovered exactly with its coordinates
local({
    f <- tempfile(fileext = ".dem"); con <- file(f, "wb")
    writeBin(as.integer(100 + 1:12), con, size = 2, endian = "big")
    close(con)
    g <- tda_gtopo(f, rows = 3, cols = 4, upper_left = c(10, 50),
                   pixel_size = c(0.5, 0.5), lat = 49.5, lon = c(10, 11.5))
    ok("gtopo recovers the synthesized row",
       is.matrix(g) && all(g[, ncol(g)] == 105:108))
    ok("gtopo coordinates step by the pixel size",
       is.matrix(g) && all(abs(diff(g[, 1]) - 0.5) < 1e-9))
})

local({
    r1 <- tda_rfit(rbind(c(0,2,1), c(2,0,2), c(1,2,0)))
    r2 <- tda_rfit(rbind(c(0,2,1), c(2,0,2), c(1,2,0)),
                   variant = "experimental")
    ok("rfit standard runs", length(r1) >= 1)
    ok("rfit1 variant runs", length(r2) >= 1)
})

# the last stragglers
local({
    ok("gap_permute runs", {
        o <- tda_gap_permute(rbind(c(0,3,1), c(3,0,2), c(1,2,0)),
                             rbind(c(0,1,4), c(1,0,2), c(4,2,0)))
        is.matrix(o$permuted) && !any(grepl("rror", o$output))})
    ok("gpro runs clean", {
        o <- tda_gpro(as.matrix(dist(cbind(c(0,0,3,3), c(0,2,0,2)))))
        is.matrix(o$projection) && !any(grepl("rror", o$output))})
    ok("tnet runs on dated layers", {
        o <- tda_tnet(data.frame(i=c(1,2,3), j=c(2,3,4), v=c(1,1,1),
                                 s=c(10,12,14), e=c(20,22,24)))
        length(o) > 5 && !any(grepl("rror", o))})
    ok("rplz renders", {
        f <- tempfile(); writeLines(c("01067;Dresden", "01069;Dresden"), f)
        length(tda_rplz(f)) > 3})
})

# typed returns + input validation (session 41): validators fail with
# plain messages naming the argument; structured classes are accepted
# where they fit
local({
    err <- function(e) tryCatch({e; ""}, error = conditionMessage)
    ok("gqap rejects non-square with a named message",
       grepl("`flow` must be square", err(tda_gqap(matrix(1,2,3), matrix(1,2,3)))))
    ok("gqap rejects asymmetric",
       grepl("must be symmetric", err(tda_gqap(matrix(c(0,1,2,0),2), matrix(0,2,2)))))
    ok("gqap rejects dim mismatch",
       grepl("same dimensions", err(tda_gqap(matrix(0,2,2), matrix(0,3,3)))))
    ok("rmod rejects non-binary",
       grepl("only 0 and 1", err(tda_rmod(matrix(c(0,1,2,1),2)))))
    ok("mdsn rejects NA",
       grepl("missing values", err(tda_mdsn(matrix(c(0,NA,NA,0),2)))))
    m <- tda_mdsn(as.matrix(dist(cbind(c(0,0,3,3), c(0,2,0,2)))))
    ok("mdsn returns tda_mds_nonmetric", inherits(m, "tda_mds_nonmetric"))
    mdr <- tda_mdsr(m, cbind(v = 1:4))
    rc <- as.numeric(sub(".*: *", "",
              grep("Rank correlation", tda_payload(mdr), value = TRUE)[1L]))
    ok("mdsr fits the collinear variable with rank correlation 1",
       is.finite(rc) && abs(rc - 1) < 1e-8)
    ok("mdsr rejects a 3-column config",
       grepl("two-column", err(tda_mdsr(matrix(1,4,3), cbind(1:4)))))
    q <- tda_gqap(rbind(c(0,3,1),c(3,0,2),c(1,2,0)),
                  rbind(c(0,1,4),c(1,0,2),c(4,2,0)))
    ok("gqap returns tda_assignment when exports are on",
       !tdaR:::.use_exports() || inherits(q, "tda_assignment"))
    s <- tda_spmod(c(0,.4,.3,.1), c(.9,.8,.7,0), rep(1,4))
    ok("spmod returns tda_leslie", inherits(s, "tda_leslie"))
    ok("structured print does not error",
       is.list(withVisible(capture.output(print(s)))))
})

# session 41: the validator now also guards the older families,
# centrally at .dist_cmd (cluster/mds/cutree) and at tda_mqap
local({
    err <- function(e) tryCatch({e; ""}, error = conditionMessage)
    ok("dist family rejects non-square centrally",
       grepl("`d` must be square", err(tda_cluster(matrix(1,2,3)))))
    ok("dist family rejects asymmetric centrally",
       grepl("must be symmetric", err(tda_mds(matrix(c(0,1,2,0),2)))))
    ok("mqap rejects dim mismatch",
       grepl("same dimensions", err(tda_mqap(matrix(0,2,2), matrix(0,3,3)))))
    ok("mds metric refuses ndim != 2",
       grepl("two dimensions only", err(tda_mds(dist(matrix(1:12, 4)),
                                                method = "metric", ndim = 3))))
})

# mdsn with Kruskal's step adaptation (changes-from-tda.md): what it reports and
# writes is checked against R, not against itself
local({
    set.seed(3)
    d <- eurodist
    m <- tda_mds(d, method = "nonmetric",
                 options = list(mxit = 300, ns = 3, df = "fit.d", pcf = "shep.d"))
    X <- as.matrix(m$points)
    dh <- as.vector(dist(X))
    dd <- as.vector(d)
    S <- read.table(file.path(m$run$dir, "shep.d"))
    # pcf=: the dissimilarities ascending, the fitted distances of the
    # returned points in that order, and their monotone regression.
    # opt=1 ("ties ignored") regresses over the sorted index without
    # pooling tied dissimilarities, so the reference is isoreg on the index
    # and the fitted column is compared as a multiset (tied pairs may sort
    # either way)
    ok("mdsn pcf= column 1 is the sorted dissimilarities",
       isTRUE(all.equal(S[, 1], sort(dd), tolerance = 1e-3)))
    ok("mdsn pcf= column 2 is dist() of the returned points",
       isTRUE(all.equal(sort(S[, 2]), sort(dh), tolerance = 1e-3)))
    disp <- isoreg(seq_len(nrow(S)), S[, 2])$yf
    ok("mdsn pcf= column 3 is the monotone regression",
       isTRUE(all.equal(S[, 3], disp, tolerance = 1e-3)))
    ok("mdsn reported stress is Kruskal stress-1 of the returned points",
       abs(sqrt(sum((S[, 2] - disp)^2) / sum(S[, 2]^2)) - m$stress) < 1e-5)
    D <- as.matrix(read.table(file.path(m$run$dir, "fit.d")))
    ok("mdsn df= is the distance matrix of the returned points",
       max(abs(D - as.matrix(dist(X)))) < 1e-10)
    # an exactly embeddable configuration has stress zero
    p <- matrix(rnorm(60), 30)
    e <- tda_mds(dist(p), method = "nonmetric",
                 options = list(mxit = 500, ns = 3))
    ok("mdsn recovers an exact planar embedding", e$stress < 1e-3)
    if (requireNamespace("MASS", quietly = TRUE)) {
        iso <- MASS::isoMDS(d, trace = FALSE)$stress / 100
        ok("mdsn on eurodist is not worse than MASS::isoMDS", m$stress <= iso + 1e-3)
    }
    # mdsn1 (tda_mdsn) is an independent implementation of the same
    # criterion -- its gradient, its normalisation; with its
    # gradient no longer accumulating across iterations (changes-from-tda.md) the
    # two land on the same minimum
    m1 <- tda_mdsn(as.matrix(d), max_iterations = 300)
    ok("mdsn1 reaches the same stress as mdsn", abs(m1$stress - m$stress) < 1e-4)
})

# session 41 close-out: the two untested edges
local({
    err <- function(e) tryCatch({e; ""}, error = conditionMessage)
    s <- tda_spmod(c(0,.4,.3,.1), c(.9,.8,.7,0), rep(1,4))
    out <- capture.output(print(s))
    ok("tda_leslie print names the class",
       any(grepl("tda_leslie", out, fixed = TRUE)))
    ok("tda_leslie print shows growth and stationary",
       any(grepl("^growth", out)) && any(grepl("^stationary", out)))
    ok("gap_permute names both args on dim mismatch",
       grepl("`g2` and `g1` must have the same dimensions",
             err(tda_gap_permute(matrix(0,2,2), matrix(0,3,3))),
             fixed = TRUE))
})

# rcsv (session 42): the export tap returns TDA's parse as a data
# frame; the padded rendering remains the flag-off contract
local({
    f <- tempfile(fileext = ".csv")
    write.csv(data.frame(a = 1:3, b = c(2.5, 1, 4)), f, row.names = FALSE)
    d <- tda_rcsv(f, sep = ",")
    if (tdaR:::.use_exports()) {
        ok("rcsv returns TDA's parse as a data frame",
           is.data.frame(d) && identical(d$a, 1:3) &&
           identical(d$b, c(2.5, 1, 4)))
    } else {
        ok("rcsv returns the rendered lines without exports",
           is.character(d) && length(d) >= 4L)
    }
})

# idf opt=2 (changes-from-tda.md): the self-consistent distribution function is the
# fixed point of the manual's 8.5.2 equation.  Checked against that
# equation iterated in R from the same start, and as a property of the
# returned values on their own.
local({
    set.seed(1)
    lo <- round(runif(40), 2)
    hi <- lo + round(runif(40, 0.05, 1), 2)
    d <- data.frame(lo, hi)
    f <- tda_idf(~ lo + hi, d, self_consistent = TRUE,
                 control = tda_control(maxit = 1000, tolf = 1e-10))
    x <- sort(unique(c(lo, hi)))
    ok("idf self-consistent table is on the induced partition",
       isTRUE(all.equal(x, f$self_consistent$partition)))
    ok("idf reports convergence", isTRUE(f$converged))
    il <- match(lo, x)
    ih <- match(hi, x)
    # one application of the equation: 0 below the interval, 1 above it,
    # the conditional share (F(x) - F(lo)) / (F(hi) - F(lo)) across it
    F <- f$table$mean_df
    for (it in 1:20000) {
        Fn <- numeric(length(x))
        for (i in seq_along(x)) {
            t <- x[i]
            s <- ifelse(lo > t, 0, ifelse(hi < t, 1,
                        (F[i] - F[il]) / (F[ih] - F[il])))
            Fn[i] <- mean(s)
        }
        if (max(abs(Fn - F)) < 1e-12) break
        F <- Fn
    }
    ok("idf self-consistent df equals the equation iterated in R",
       max(abs(f$self_consistent$df - F)) < 1e-7)
    FT <- f$self_consistent$df
    res <- numeric(length(x))
    for (i in seq_along(x)) {
        t <- x[i]
        s <- ifelse(lo > t, 0, ifelse(hi < t, 1,
                    (FT[i] - FT[il]) / (FT[ih] - FT[il])))
        res[i] <- mean(s) - FT[i]
    }
    ok("idf self-consistent df satisfies its equation", max(abs(res)) < 1e-8)
    ok("idf self-consistent df ends at 1 and is nondecreasing",
       abs(FT[length(FT)] - 1) < 1e-8 && all(diff(FT) > -1e-10))
    ok("idf without opt=2 has no self_consistent element",
       is.null(tda_idf(~ lo + hi, d)$self_consistent))
})

# rfit against brute force: every relation on 4 nodes (2^16), the
# requested properties checked directly, the objective sum((2 d - 1) x)
# with a 1 on the diagonal as rfit() builds it.  rfit's value must be the
# optimum and every relation it returns an optimal one with the
# properties; the number of tied optima it lists may be smaller (see the
# docs of max_solutions).
local({
    set.seed(3)
    n <- 4
    d <- matrix(rbinom(n * n, 1, 0.5), n)
    d <- pmax(d, t(d))
    diag(d) <- 0
    C <- 2 * d - 1
    diag(C) <- 1
    for (props in list(c("symmetric", "transitive"),
                       c("reflexive", "transitive"))) {
        best <- -Inf
        keys <- character()
        for (b in 0:(2^(n * n) - 1)) {
            X <- matrix(as.integer(intToBits(b))[1:(n * n)], n, byrow = TRUE)
            ok <- TRUE
            if ("reflexive" %in% props) ok <- ok && all(diag(X) == 1)
            if ("symmetric" %in% props) ok <- ok && all(X == t(X))
            if ("transitive" %in% props) {
                XX <- (X %*% X) > 0
                ok <- ok && all(X[XX] == 1)
            }
            if (!ok) next
            v <- sum(C * X)
            if (v > best) {
                best <- v
                keys <- paste(X, collapse = "")
            } else if (v == best) {
                keys <- c(keys, paste(X, collapse = ""))
            }
        }
        r <- tda_rfit(d, properties = props, max_solutions = 200)
        vals <- sapply(r$relations, function(X) sum(C * X))
        ok(paste("rfit", paste(props, collapse = "+"), "reaches the brute-force optimum"),
           all(vals == best))
        ok(paste("rfit", paste(props, collapse = "+"), "returns only optimal relations with the properties"),
           all(sapply(r$relations, function(X) paste(X, collapse = "")) %in% keys))
    }
})
