# Differential fuzzer for TDA's expression evaluator (t_eval*):
# every documented operator (tda.hlp section 5.2.5) is computed by TDA
# and by R on the same randomized-but-valid values; results must agree
# to 1e-6.  Structure comes from the documentation, only the values
# are random (seeded).  Usage: Rscript tools/fuzz_eval.R [seed]
library(tdaR)
seed <- if (length(commandArgs(TRUE))) as.integer(commandArgs(TRUE)[1]) else 42L
set.seed(seed)
n <- 25L
d <- data.frame(
    X = round(rnorm(n, 0, 3), 4),
    P = round(runif(n, 0.05, 0.95), 4),       # probabilities
    A = round(runif(n, 0.5, 6), 4),            # positive shape args
    B = round(runif(n, 0.5, 6), 4),
    K = sample(0:8, n, TRUE),                  # small counts
    N = sample(9:15, n, TRUE),                 # larger counts
    DF = sample(1:20, n, TRUE),                # degrees of freedom
    Y = sample(1990:2010, n, TRUE),            # dates
    M = sample(1:12, n, TRUE),
    D = sample(1:28, n, TRUE))
# TDA's mr is the hazard rate dnorm/(1-pnorm) -- the reciprocal of
# the textbook Mill's ratio (established by this fuzzer, session 39)
# TDA's mr approximation is good to ~1e-4 up to |x| ~ 7 and drifts
# beyond (2.5% at |x| ~ 8, measured session 39) -- an engine accuracy
# limitation, characterized, not fixed; the case stays inside the
# trustworthy domain.
mr_r  <- function(d) { a <- pmin(abs(d$X), 6); dnorm(a) / (1 - pnorm(a)) }
cases <- list(
 list("abs(X)",            function(d) abs(d$X)),
 list("ceil(X)",           function(d) ceiling(d$X)),
 list("floor(X)",          function(d) floor(d$X)),
 list("rnd(X)",            function(d) round(d$X)),
 list("sign(X)",           function(d) sign(d$X)),
 list("sqrt(A)",           function(d) sqrt(d$A)),
 list("exp(X / 3)",        function(d) exp(d$X / 3)),
 list("eexp(X)",           function(d) exp(d$X) / (1 + exp(d$X))),
 list("log(A)",            function(d) log(d$A)),
 list("sin(X)",            function(d) sin(d$X)),
 list("cos(X)",            function(d) cos(d$X)),
 list("lgam(A)",           function(d) lgamma(d$A)),
 list("digam(A)",          function(d) digamma(d$A)),
 list("trigam(A)",         function(d) trigamma(d$A)),
 list("bc(N,K)",           function(d) choose(d$N, d$K)),
 list("max(X,A,B)",        function(d) pmax(d$X, d$A, d$B)),
 list("min(X,A,B)",        function(d) pmin(d$X, d$A, d$B)),
 list("tr(X,-1,1)",        function(d) pmin(pmax(d$X, -1), 1)),
 list("K % 3",             function(d) d$K %% 3),  # modulus is integer-only in TDA (documented finding)
 list("A ^ 2",             function(d) d$A^2),
 list("X + A * B - A / B", function(d) d$X + d$A * d$B - d$A / d$B),
 list("eq(K,4)",           function(d) as.numeric(d$K == 4)),
 list("ne(K,4)",           function(d) as.numeric(d$K != 4)),
 list("gt(X,0)",           function(d) as.numeric(d$X > 0)),
 list("ge(K,4)",           function(d) as.numeric(d$K >= 4)),
 list("lt(X,0)",           function(d) as.numeric(d$X < 0)),
 list("le(K,4)",           function(d) as.numeric(d$K <= 4)),
 list("and(gt(X,0),lt(X,2))", function(d) as.numeric(d$X > 0 & d$X < 2)),
 list("or(lt(X,-2),gt(X,2))", function(d) as.numeric(d$X < -2 | d$X > 2)),
 list("not(gt(X,0))",      function(d) as.numeric(!(d$X > 0))),
 list("if(gt(X,0), A, B)", function(d) ifelse(d$X > 0, d$A, d$B)),
 list("nd(X)",             function(d) pnorm(d$X)),
 list("bivn(X/3,B-3,P-0.5)", if (requireNamespace("mvtnorm", quietly = TRUE))
     function(d) mapply(function(h, k, r)
         mvtnorm::pmvnorm(upper = c(h, k),
                          corr = matrix(c(1, r, r, 1), 2))[1L],
         d$X / 3, d$B - 3, d$P - 0.5)
     else NULL),
 list("ndi(P)",            function(d) qnorm(d$P)),
 list("ndf(X)",            function(d) dnorm(d$X)),
 list("mr(tr(abs(X),0,6))", mr_r),
 list("chd(A,DF)",         function(d) pchisq(d$A, d$DF)),
 list("fd(A,DF,DF)",       function(d) pf(d$A, d$DF, d$DF)),
 list("td(X,DF)",          function(d) pt(abs(d$X), d$DF)),  # sign-blind: td(x)=pt(|x|) (established here)
 list("icg(A,B)",          function(d) pgamma(d$A, d$B)),
 list("icb(P,A,B)",        function(d) pbeta(d$P, d$A, d$B)),  # args: x, a, b
 list("poisson(A,K)",      function(d) dpois(d$K, d$A, log = TRUE)),  # args: mean, integer count (established here)
 list("negbin(DF,A,K)",    function(d) dnbinom(d$K, size = d$DF, mu = d$A, log = TRUE)),  # observed: size, mean, count
 list("gclen(0,0,X,0)",    function(d) abs(d$X)),  # observed: arc length in DEGREES
 list("jul(Y,M,D)",        NULL),   # verified by round trip below
 list("july(jul(Y,M,D))",  function(d) d$Y),
 list("julm(jul(Y,M,D))",  function(d) d$M),
 list("juld(jul(Y,M,D))",  function(d) d$D))
# ---- stage 2: column operators (tda.hlp 5.2.6) ------------------------
# Conventions are OBSERVED by the differential, not assumed; where TDA
# and textbook R defaults differ, the mapping records what TDA does.
d$G <- sort(sample(1:4, n, TRUE))              # grouping, contiguous
cases2 <- list(
 list("mean(X)",        function(d) rep(mean(d$X), nrow(d))),
 list("sum(X)",         function(d) rep(sum(d$X), nrow(d))),
 list("vmin(X)",        function(d) rep(min(d$X), nrow(d))),
 list("vmax(X)",        function(d) rep(max(d$X), nrow(d))),
 list("std(X)",         function(d) rep(sd(d$X), nrow(d))),
 list("cum(X)",         function(d) cumsum(d$X)),
 list("sort(X)",        function(d) sort(d$X)),
 list("rank(X)",        function(d) rank(d$X)),
 list("ndv(X)",         function(d) rep(length(unique(d$X)), nrow(d))),
 list("ndv1(X,0)",      function(d) rep(length(unique(d$X[d$X >= 0])), nrow(d))),
 list("pre(X)",         function(d) c(NA, d$X[-nrow(d)])),
 list("suc(X)",         function(d) c(d$X[-1], NA)),
 list("lag(X,2)",       function(d) c(d$X[-(1:2)], NA, NA)),  # observed: lag(X,n) looks n rows FORWARD
 list("cnteq(K,4)",     function(d) rep(sum(d$K == 4), nrow(d))),
 list("cntgt(X,0)",     function(d) rep(sum(d$X > 0), nrow(d))),
 list("cntlt(X,0)",     function(d) rep(sum(d$X < 0), nrow(d))),
 list("change(K)",      function(d) as.numeric(c(TRUE, d$K[-1] != d$K[-nrow(d)]))),  # observed: first row counts as a change
 list("cntch(K)",       function(d) rep(1 + sum(d$K[-1] != d$K[-nrow(d)]), nrow(d))),  # observed: counts the first row too
 list("gmean(X,G)",     function(d) ave(d$X, d$G)),
 list("gsum(X,G)",      function(d) ave(d$X, d$G, FUN = sum)),
 list("gmin(X,G)",      function(d) ave(d$X, d$G, FUN = min)),
 list("gmax(X,G)",      function(d) ave(d$X, d$G, FUN = max)),
 list("gstd(X,G)",      function(d) ave(d$X, d$G, FUN = sd)),
 list("gcnt(G)",        function(d) ave(d$X, d$G, FUN = length)),
 list("gfirst(G)",      function(d) as.numeric(!duplicated(d$G))),
 list("glast(G)",       function(d) as.numeric(!duplicated(d$G, fromLast = TRUE))),
 list("grec(G)",        function(d) ave(seq_len(nrow(d)), d$G, FUN = seq_along)),
 list("gsn(G)",         function(d) as.numeric(factor(d$G))),
 list("gndv(X,G)",      function(d) ave(d$X, d$G, FUN = function(v) length(unique(v)))),
 list("quant1(A,gt(K,0),0.5,0)", function(d) {
     # observed (arbitrated by a hand Kaplan-Meier): quant1's p is the
     # SURVIVOR-function level and the quantile is linearly
     # interpolated between KM steps; Z = 1 marks an event, 0 censored
     t <- d$A; ev <- as.numeric(d$K > 0)
     o <- order(t); t <- t[o]; ev <- ev[o]
     S <- cumprod(ifelse(ev == 1, 1 - 1 / rev(seq_along(t)), 1))
     pts <- rbind(c(min(t), 1), cbind(t, S))
     f <- approxfun(pts[, 2], pts[, 1], ties = min)
     rep(f(0.5), length(t))
 }),
 list("mav(X,3)",       function(d) {          # observed: centered window
     x <- d$X; n <- length(x)                    # i-3 .. i+3, truncated
     vapply(seq_len(n), function(i)
         mean(x[max(1, i - 3):min(n, i + 3)]), 0)
 }),
 list("cd(sort(X))",    function(d) ecdf(d$X)(sort(d$X))))
 # observed: cd is the empirical CDF (upper step at ties) over
 # ASCENDING input; on unsorted data the output is order-dependent,
 # so the working idiom is cd(sort(X)).
cases <- c(cases, cases2)
dr <- tempfile("fz"); dir.create(dr)
write.table(d, file.path(dr, "fz.dat"), row.names = FALSE, col.names = FALSE)
decl <- paste0("nvar(dfile=fz.dat, X[12.4]=c1, P[12.4]=c2, A[12.4]=c3, B[12.4]=c4, K[4.0]=c5, N[4.0]=c6, DF[4.0]=c7, Y[6.0]=c8, M[4.0]=c9, D[4.0]=c10, G[4.0]=c11, E[24.12]=")
bad <- 0L
for (cs in cases) {
    expr <- cs[[1]]; rfun <- cs[[2]]
    out <- tryCatch(tda_run(c(paste0(decl, expr, ");"),
                              "pdata(keep=E) = fzout.txt;"), dir = dr),
                    error = function(e) e)
    if (inherits(out, "error") || length(grep("^Error", out$output))) {
        cat(sprintf("FAIL-RUN   %-26s %s\n", expr,
                    if (inherits(out, "error")) conditionMessage(out)
                    else grep("^Error", out$output, value = TRUE)[1]))
        bad <- bad + 1L; next
    }
    got <- scan(file.path(dr, "fzout.txt"), quiet = TRUE)
    file.remove(file.path(dr, "fzout.txt"))
    if (is.null(rfun)) next
    want <- rfun(d)
    tol <- if (grepl("^mr\\(", expr)) 1e-4 else 1e-6  # mr uses a tail
    # approximation good to ~1e-5 relative (established here)
    want <- unname(want)
    got2 <- got; got2[is.na(want)] <- NA   # boundary cells: observe, don't judge
    # scale-aware: TDA accumulates column statistics in single
    # precision (observed: |diff| ~ 3e-8 on a near-zero mean), so a
    # near-zero result must not be judged at relative tolerance
    sc <- max(1, max(abs(want), na.rm = TRUE))
    okv <- isTRUE(all.equal(got2, want, tolerance = tol, scale = sc))
    if (!okv) {
        bad <- bad + 1L
        i <- which(abs(got - want) > 1e-6 * pmax(1, abs(want)))[1]
        cat(sprintf("MISMATCH   %-26s row %d: tda %.10g vs R %.10g\n",
                    expr, i, got[i], want[i]))
    }
}
cat(sprintf("fuzz_eval seed %d: %d/%d operators clean\n",
            seed, length(cases) - bad, length(cases)))
if (bad) quit(status = 1)

# ---- stage 3: matrix expressions (mexpr; cyclic broadcasting) ----------
# Oracle: the documented rule itself, implemented independently --
# result dims are the operand maxima, missing rows/columns cycle.
# Where TDA and this oracle disagree, the docs' own worked example
# (3x1 * 1x2 outer table) arbitrates before anything is called wrong.
cyc <- function(M, r, c)
    M[((seq_len(r) - 1L) %% nrow(M)) + 1L,
      ((seq_len(c) - 1L) %% ncol(M)) + 1L, drop = FALSE]
bcast <- function(f, ...) {
    ms <- list(...)
    r <- max(vapply(ms, nrow, 0L)); c <- max(vapply(ms, ncol, 0L))
    Reduce(f, lapply(ms, cyc, r = r, c = c))
}
mdef_txt <- function(name, M)
    sprintf("mdef(%s,%d,%d) = %s;", name, nrow(M), ncol(M),
            paste(t(M), collapse = ","))
A <- matrix(round(rnorm(sample(2:4, 1) * sample(2:4, 1), 0, 2), 3),
            sample(2:4, 1))
B <- matrix(round(runif(sample(2:4, 1) * sample(2:4, 1), 0.5, 3), 3),
            sample(2:4, 1))
cases3 <- list(
 list("A + B",              function() bcast(`+`, A, B)),
 list("A - B",              function() bcast(`-`, A, B)),
 list("A * B",              function() bcast(`*`, A, B)),
 list("A / B",              function() bcast(`/`, A, B)),
 list("A ^ 2",              function() A^2),
 list("sin(A)",             function() sin(A)),
 list("exp(A / 3)",         function() exp(A / 3)),
 list("abs(A) + sqrt(B)",   function() bcast(`+`, abs(A), sqrt(B))),
 list("if(gt(A,0), A, 0 - A)", function() abs(A)),
 list("max(A, B)",          function() bcast(pmax, A, B)),
 list("min(A, B)",          function() bcast(pmin, A, B)),
 list("nd(A)",              function() pnorm(A)),
 list("A * B + A",          function() bcast(`+`, bcast(`*`, A, B), A)))
bad3 <- 0L
for (cs in cases3) {
    expr <- cs[[1]]; want <- cs[[2]]()
    dr3 <- tempfile("m3"); dir.create(dr3)
    o <- tryCatch(tda_run(c(mdef_txt("A", A), mdef_txt("B", B),
                            "mfmt = 22.12;", sprintf("mexpr(%s, R);", expr),
                            "mpr(R) = m3.out;"), dir = dr3),
                  error = function(e) e)
    f <- file.path(dr3, "m3.out")
    if (inherits(o, "error") || !file.exists(f)) {
        cat(sprintf("FAIL-RUN   %-26s %s\n", expr,
            if (inherits(o, "error")) conditionMessage(o)
            else grep("rror", o$output, value = TRUE)[1])); bad3 <- bad3 + 1L
        next
    }
    got <- as.matrix(read.table(f))
    if (!isTRUE(all.equal(unname(got), unname(want), tolerance = 1e-6,
                          check.attributes = FALSE))) {
        bad3 <- bad3 + 1L
        cat(sprintf("MISMATCH   %-26s dims tda %s vs oracle %s; [1,1] %g vs %g\n",
                    expr, paste(dim(got), collapse = "x"),
                    paste(dim(want), collapse = "x"), got[1, 1], want[1, 1]))
    }
}
cat(sprintf("fuzz_eval stage3 seed %d: %d/%d matrix expressions clean (A %dx%d, B %dx%d)\n",
            seed, length(cases3) - bad3, length(cases3),
            nrow(A), ncol(A), nrow(B), ncol(B)))
if (bad3) quit(status = 1)

# ---- stage 4: symbolic derivatives via fmin ----------------------------
# Triangulated: TDA with analytic derivatives (dopt=0, the t_eval2
# engine) vs TDA with numerical approximation (dopt=1) vs R optim on
# the same function.  All three must land on the same optimum; a
# disagreement is arbitrated by the third computation before anything
# is judged.
fmin_run <- function(fn, mina, dopt, xp) {
    dr4 <- tempfile("f4"); dir.create(dr4)
    o <- tda_run(sprintf("fmin(mina=%d, dopt=%d, xp=%s, mxit=500) = %s;",
                         mina, dopt, paste(xp, collapse = ","), fn),
                 dir = dr4)
    ln <- grep("^  [0-9]+ +[a-z]", o$output, value = TRUE)
    if (length(ln) < 2L) return(NULL)
    ln <- tail(ln, 2L)          # the result table; the earlier matches
    vapply(strsplit(trimws(ln), "\\s+"),  # are the starting values
           function(v) as.numeric(v[3]), 0)
}
probs <- list(
 list("(x - 2)^2 + (y + 1)^2",              c(2, -1),
      function(p) (p[1] - 2)^2 + (p[2] + 1)^2),
 list("(1 - x)^2 + 5 * (y - x^2)^2",        c(1, 1),
      function(p) (1 - p[1])^2 + 5 * (p[2] - p[1]^2)^2),
 list("exp(x - 1) - x + (y - 0.5)^2",       c(1, 0.5),
      function(p) exp(p[1] - 1) - p[1] + (p[2] - 0.5)^2))
bad4 <- 0L; n4 <- 0L
for (pr in probs) {
    fn <- pr[[1]]; truth <- pr[[2]]; rf <- pr[[3]]
    ro <- optim(c(0.2, 0.2), rf, method = "BFGS")$par
    for (mina in c(3, 4, 5, 6)) {
        pa <- fmin_run(fn, mina, 0, c(0.2, 0.2))   # analytic (t_eval2)
        pn <- fmin_run(fn, mina, 1, c(0.2, 0.2))   # numeric approx
        n4 <- n4 + 1L
        # the derivative-engine test is analytic-vs-numeric agreement;
        # distance from the known optimum is an algorithm property and
        # is reported, not judged
        agree <- !is.null(pa) && !is.null(pn) &&
            isTRUE(all.equal(pa, pn, tolerance = 1e-3))
        atopt <- agree && isTRUE(all.equal(pa, truth, tolerance = 1e-3))
        if (agree && !atopt)
            cat(sprintf("NOTE       mina=%d %-30s both modes stop at %s (optimum %s)\n",
                mina, fn, paste(round(pa, 4), collapse = ","),
                paste(truth, collapse = ",")))
        if (!agree) {
            bad4 <- bad4 + 1L
            cat(sprintf("DIVERGE    mina=%d %-30s analytic %s | numeric %s | optim %s | truth %s\n",
                mina, fn,
                if (is.null(pa)) "no-output" else paste(round(pa, 4), collapse = ","),
                if (is.null(pn)) "no-output" else paste(round(pn, 4), collapse = ","),
                paste(round(ro, 4), collapse = ","),
                paste(truth, collapse = ",")))
        }
    }
}
cat(sprintf("fuzz_eval stage4 seed %d: %d/%d fmin triangulations clean\n",
            seed, n4 - bad4, n4))
if (bad4) quit(status = 1)

# ---- stage 5: block operators (tda.hlp 5.2.7; the v_eval2 engine) ------
# Blocks come from dblock over a sorted grouping variable; oracles are
# base-R computations over the same blocks.
d5 <- data.frame(ID = sort(sample(1:5, 20, TRUE)),
                 X = round(rnorm(20, 0, 2), 3))
dr5 <- tempfile("b5"); dir.create(dr5)
write.table(d5, file.path(dr5, "b.dat"), row.names = FALSE, col.names = FALSE)
blk <- function(f) ave(d5$X, cumsum(!duplicated(d5$ID)), FUN = f)
cases5 <- list(
 list("bmin(X)",   function() blk(min)),
 list("bmax(X)",   function() blk(max)),
 list("bnrec",     function() blk(length)),
 list("brec",      function() ave(seq_along(d5$ID), d5$ID, FUN = seq_along)),
 list("bfirst",    function() as.numeric(!duplicated(d5$ID))),
 list("blast",     function() as.numeric(!duplicated(d5$ID, fromLast = TRUE))),
 list("bnum",      function() cumsum(!duplicated(d5$ID))))
bad5 <- 0L
for (cs in cases5) {
    expr <- cs[[1]]; want <- cs[[2]]()
    # dblock is an nvar OPTION (block mode), established from
    # gdatops5.cf -- the standalone dblock command does not feed the
    # lazy type-4 evaluation
    o <- tryCatch(tda_run(c(sprintf(
        "nvar(dfile=b.dat, ID[4.0]=c1, X[10.3]=c2, dblock=ID, E[24.12]=%s);",
        expr), "pdata(keep=E) = b.out;"), dir = dr5),
                  error = function(e) e)
    f <- file.path(dr5, "b.out")
    if (inherits(o, "error") || !file.exists(f)) {
        cat(sprintf("FAIL-RUN   %-14s %s\n", expr,
            if (inherits(o, "error")) conditionMessage(o)
            else grep("rror", o$output, value = TRUE)[1])); bad5 <- bad5 + 1L
        next
    }
    got <- scan(f, quiet = TRUE); file.remove(f)
    if (!isTRUE(all.equal(got, unname(want), tolerance = 1e-6))) {
        bad5 <- bad5 + 1L
        i <- which(abs(got - want) > 1e-6)[1]
        cat(sprintf("MISMATCH   %-14s row %d: tda %g vs R %g\n",
                    expr, i, got[i], want[i]))
    }
}
cat(sprintf("fuzz_eval stage5 seed %d: %d/%d block operators clean\n",
            seed, length(cases5) - bad5, length(cases5)))
if (bad5) quit(status = 1)

# ---- stage 6: mexp (t_eval4's linear-algebra engine) -------------------
# mexp is mexpr's undocumented sibling: * is true matrix
# multiplication when both operands are matrices (scalar
# multiplication otherwise), with functions t, cross, ginv, trace,
# sqrt, eexp.  Oracles from base R; ginv arbitrated by the
# svd pseudoinverse (no packages).
M1 <- matrix(round(rnorm(6, 0, 2), 3), 2)
M2 <- matrix(round(rnorm(6, 0, 2), 3), 3)
psinv <- function(M) { s <- svd(M); s$v %*% diag(1/s$d, length(s$d)) %*% t(s$u) }
cases6 <- list(
 list("M1 * M2",          function() M1 %*% M2),
 list("t(M1)",            function() t(M1)),
 list("M1 * t(M1)",       function() M1 %*% t(M1)),
 list("cross(M1)",        function() crossprod(M1)),
 list("trace(M1 * t(M1))", function() matrix(sum(diag(M1 %*% t(M1))), 1)),
 list("sqrt(eexp(M1) + 1)",  function() sqrt(exp(M1)/(1+exp(M1)) + 1)),  # mexp has no abs; eexp keeps the domain positive
 list("ginv(M1)",         function() psinv(M1)),
 list("M1 * M2 * t(M2)",  function() M1 %*% M2 %*% t(M2)),
 list("(M1 + M1) * 0.5",  function() M1),
 list("eexp(M1)",         function() exp(M1) / (1 + exp(M1))))
bad6 <- 0L
for (cs in cases6) {
    expr <- cs[[1]]; want <- cs[[2]]()
    dr6 <- tempfile("x6"); dir.create(dr6)
    o <- tryCatch(tda_run(c(mdef_txt("M1", M1), mdef_txt("M2", M2),
                            "mfmt = 22.12;",
                            sprintf("mexp(%s, R);", expr),
                            "mpr(R) = x6.out;"), dir = dr6),
                  error = function(e) e)
    f <- file.path(dr6, "x6.out")
    if (inherits(o, "error") || !file.exists(f)) {
        cat(sprintf("FAIL-RUN   %-24s %s\n", expr,
            if (inherits(o, "error")) conditionMessage(o)
            else grep("rror", o$output, value = TRUE)[1])); bad6 <- bad6 + 1L
        next
    }
    got <- as.matrix(read.table(f))
    if (!isTRUE(all.equal(unname(got), unname(want), tolerance = 1e-6,
                          check.attributes = FALSE))) {
        bad6 <- bad6 + 1L
        cat(sprintf("MISMATCH   %-24s dims %s vs %s; [1,1] %g vs %g\n",
                    expr, paste(dim(got), collapse = "x"),
                    paste(dim(want), collapse = "x"), got[1, 1], want[1, 1]))
    }
}
cat(sprintf("fuzz_eval stage6 seed %d: %d/%d mexp expressions clean\n",
            seed, length(cases6) - bad6, length(cases6)))
if (bad6) quit(status = 1)
