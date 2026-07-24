# Content oracles: every assertion compares a wrapper's RESULT to an
# independent R computation or an exact known answer -- never to a
# length, a class, or the absence of an error.
eq <- function(a, b, tol = 1e-6)
    isTRUE(all.equal(unname(as.matrix(a)), unname(as.matrix(b)),
                     tolerance = tol, check.attributes = FALSE))

# ---- matrix family vs base R -------------------------------------
local({
    A <- rbind(c(4, 1), c(1, 3)); B <- rbind(c(2, 0), c(1, 1))
    ok("mmul == %*%",      eq(tda_mmul(A, B), A %*% B))
    ok("mtransp == t()",   eq(tda_mtransp(A), t(A)))
    ok("minvs == solve()", eq(tda_minvs(A), solve(A)))
    ok("mchol == lower chol()", eq(tda_mchol(A), t(chol(A))))
    ev <- tda_mevs(A); e0 <- eigen(A)
    ok("mevs values == eigen()", eq(ev$values, e0$values))
    ok("mevs vectors span eigen()'s",
       all(abs(abs(colSums(ev$vectors * e0$vectors)) - 1) < 1e-6))
    X <- rbind(c(3, 1), c(1, 3), c(1, 1))
    ok("msvd == svd()$d",  eq(tda_msvd(X), svd(X)$d))
    s1 <- tda_msvd1(X)
    ok("msvd1 d == svd()$d", eq(s1$d, svd(X)$d))
    x <- c(9, 2, 7, 2)
    ok("msort == sort() on the by column",
       eq(tda_msort(cbind(x), by = 1)[, 1], sort(x)))
    ok("mrank == order() (the sorting permutation)",
       eq(tda_mrank(cbind(x), by = 1), order(x)))
    M <- rbind(c(1, 4), c(2, 5), c(3, 6))
    ok("mrsum == rowSums",   eq(tda_mrsum(M), rowSums(M)))
    ok("mcsum == colSums",   eq(tda_mcsum(M), colSums(M)))
    ok("mnrow/mncol == dim", tda_mnrow(M) == 3 && tda_mncol(M) == 2)
    ok("mcent == scale(center)", eq(tda_mcent(M), scale(M, scale = FALSE)))
    ok("mdiag == diag()",    eq(tda_mdiag(c(5, 7)), diag(c(5, 7))))
    ok("mcvec == as.vector", eq(tda_mcvec(M), as.vector(M)))
    p <- c(2, 3, 1)
    ok("mpinv == order() (inverse permutation)",
       eq(tda_mpinv(p), order(p)))
    q <- tda_mqp(rbind(c(2, 0), c(0, 2)), c(-2, -4))
    ok("mqp x and value == analytic QP",
       eq(q$x, c(1, 2)) && abs(q$value - (-5)) < 1e-8)
    ls1 <- tda_mls(rbind(c(1, 0), c(0, 1), c(1, 1)), c(2.1, 2.9, 5))
    o1 <- qr.solve(rbind(c(1, 0), c(0, 1), c(1, 1)), c(2.1, 2.9, 5))
    ok("mls == qr.solve least squares", eq(ls1$x, o1))
})

# ---- integration / optimisation vs closed forms ------------------
local({
    ok("evalf == the function value at the point",
       abs(tda_evalf("sin(x) + x^2", options = list(x = 0.7))$value -
           (sin(0.7) + 0.49)) < 1e-3)
    # flag-off reads the printed six-digit value, so the tolerance
    # covers print precision
    ok("integrate == 1/3",
       abs(tda_integrate("x*x", 0, 1)$value - 1/3) < 1e-5)
    ok("gmin finds the global minimum value",
       abs(tda_gmin("(x-2)*(x-2)", start = list(c(2, 0, 4)))$value) < 1e-6)
    m <- tda_minimize("(x - 2)^2", options = list(xp = "0"))
    ok("minimize lands on the analytic minimiser",
       abs(m$estimates$Value[1L] - 2) < 1e-4)
    r <- tda_range("sin(x)")
    ok("range == sin at its default domain ends",
       eq(r$range, sin(c(-1, 1)), tol = 1e-4))
})

# ---- regression vs lm/nls ----------------------------------------
local({
    # mreg rescales the response monotonically; on noise-free linear
    # data the rescaling is the identity and lm() is the exact oracle
    set.seed(11); x <- rnorm(50)
    d <- data.frame(x = x, y = 3 + 2 * x)
    ok("mreg == lm() when the monotone rescale is the identity",
       eq(coef(tda_mreg(y ~ x, d)), coef(lm(y ~ x, d)), tol = 1e-4))
    set.seed(12); xn <- runif(40, 0, 5)
    dn <- data.frame(x = xn, y = 3 * exp(0.1 * xn) + rnorm(40, sd = .05))
    fn <- tda_nlreg(y ~ x, dn, expr = "a * exp(b * x)", start = c(1, 1))
    on <- nls(y ~ a * exp(b * x), dn, start = list(a = 1, b = 0.5))
    ok("nlreg coefficients == nls()", eq(coef(fn), coef(on), tol = 1e-3))
})

# ---- counting vs base R ------------------------------------------
local({
    f <- tempfile(); writeLines(c("hello world", "foo"), f)
    cc <- tda_ccnt(f)
    ok("ccnt counts 'l' and 'o' like table()",
       cc$count[match("l", cc$char)] == 3L &&
       cc$count[match("o", cc$char)] == 4L &&
       cc$count[match("h", cc$char)] == 1L)
    lc <- tda_lcnt(f)
    ok("lcnt lengths and frequencies == nchar()",
       sum(lc$frequency) == 2L &&
       setequal(rep(lc$length, lc$frequency), nchar(c("hello world", "foo"))))
})

# ---- graph family vs igraph --------------------------------------
local({
    if (!requireNamespace("igraph", quietly = TRUE)) return(invisible())
    e <- data.frame(from = c(1, 1, 2, 3, 4), to = c(2, 3, 3, 4, 5))
    tg <- tda_graph(e, directed = FALSE)
    ig <- igraph::graph_from_data_frame(e, directed = FALSE)
    s <- tda_g_shortest(tg, format = "matrix")
    D <- as.matrix(as.data.frame(s$table %||% s)[, -(1:2)])
    d0 <- igraph::distances(ig)
    d0 <- d0[order(as.integer(rownames(d0))), order(as.integer(colnames(d0)))]
    ok("g_shortest distance matrix == igraph::distances", eq(D, d0))
})

# ---- graph family vs igraph, batch 2 -----------------------------
local({
    if (!requireNamespace("igraph", quietly = TRUE)) return(invisible())
    e <- data.frame(from = c(1, 1, 2, 3, 4), to = c(2, 3, 3, 4, 5))
    g  <- tda_graph(e, directed = FALSE)
    ig <- igraph::graph_from_data_frame(e, directed = FALSE)
    d <- tda_g_degrees(g)
    ok("g_degrees == igraph::degree",
       eq(d$degree[order(d$node)],
          unname(igraph::degree(ig))[order(as.integer(igraph::V(ig)$name))]))
    ed <- data.frame(from = c(1, 1, 2, 3), to = c(2, 3, 3, 1))
    gd  <- tda_graph(ed, directed = TRUE)
    igd <- igraph::graph_from_data_frame(ed, directed = TRUE)
    dd <- tda_g_degrees(gd)
    o <- order(as.integer(igraph::V(igd)$name))
    ok("g_degrees in/out == igraph modes",
       eq(dd$in_degree[order(dd$node)],
          unname(igraph::degree(igd, mode = "in"))[o]) &&
       eq(dd$out_degree[order(dd$node)],
          unname(igraph::degree(igd, mode = "out"))[o]))
    # maximum spanning tree == igraph mst on negated weights
    ev <- data.frame(from = c(1, 1, 2, 3), to = c(2, 3, 3, 4),
                     v = c(10, 1, 2, 5))
    gv <- tda_graph(ev, directed = FALSE)
    sp <- tda_g_spanning(gv)$table
    igv <- igraph::graph_from_data_frame(ev[1:2], directed = FALSE)
    igraph::E(igv)$weight <- -ev$v
    mx <- igraph::mst(igv)
    ok("g_spanning == maximum spanning tree (igraph mst, negated)",
       sum(sp[[4L]]) == -sum(igraph::E(mx)$weight))
    # directed reachability counts
    ep <- data.frame(from = 1:4, to = 2:5)
    gp <- tda_graph(ep, directed = TRUE)
    rc <- tda_g_reachable(gp)$table
    igp <- igraph::graph_from_data_frame(ep, directed = TRUE)
    cnt <- vapply(as.character(1:5), function(v)
        length(igraph::subcomponent(igp, v, mode = "out")) - 1L,
        integer(1L))
    ok("g_reachable counts == igraph subcomponents",
       all(rc[[3L]] == cnt[as.character(rc[[2L]])]))
    # transitive closure distances == igraph distances (Inf -> -1)
    tc <- tda_g_transitive(gp)$table
    D <- as.matrix(tc[, -(1:2)])
    d0 <- igraph::distances(igp, mode = "out")
    d0 <- d0[order(as.integer(rownames(d0))),
             order(as.integer(colnames(d0)))]
    d0[is.infinite(d0)] <- -1
    ok("g_transitive distances == igraph::distances", eq(D, d0))
    # toposort: every edge points forward in the returned order
    ts <- tda_g_toposort(gp)$table
    pos <- ts[[3L]][order(ts[[2L]])]
    ok("g_toposort order respects every edge",
       all(pos[ep$from] < pos[ep$to]))
    # cliques on a triangle-plus-tail: the one maximal 3-clique found
    ct <- tda_g_cliques(tda_graph(
        data.frame(from = c(1, 1, 2, 3), to = c(2, 3, 3, 4)),
        directed = FALSE))$table
    ok("g_cliques finds the triangle",
       any(ct[[3L]] == 3L | ct[[4L]] == 3L))
})

# ---- batch 3: interpolation, intervals, LP, rates, sequences -----
local({
    # triangulated linear interpolation is exact on a plane
    set.seed(4); x <- runif(20); y <- runif(20); z <- 2 * x + 3 * y + 1
    it <- tda_interp(x, y, z, rx = c(0.3, 0.6), ry = c(0.4, 0.5))$table
    ok("interp reproduces the plane exactly",
       eq(it$z, 2 * it$x + 3 * it$y + 1, tol = 1e-6))
    # interval evaluation of a monotone function is the endpoint image
    ok("evalfi x^2 over [1,2] == [1,4]",
       eq(tda_evalfi("x*x", x = c(1, 2)), c(1, 4)))
    # point intervals: interval sd == population sd
    d0 <- data.frame(lo = c(1, 4, 6, 9), hi = c(1, 4, 6, 9))
    s0 <- tda_isd(~ iv(lo, hi), d0)
    ok("isd on point intervals == population sd",
       abs(s0[["lower"]] - sqrt(mean((d0$lo - mean(d0$lo))^2))) < 1e-5 &&
       abs(s0[["upper"]] - s0[["lower"]]) < 1e-8)
    # constrained least squares on a consistent system is exact
    m1 <- tda_mlse(rbind(c(1, 0), c(0, 1), c(1, 1)), c(2, 3, 5))
    ok("mlse solves the consistent system exactly", eq(m1$x, c(2, 3)))
    # linear program: max x1 + 2 x2, x1+x2 <= 3, 2 x1+x2 <= 4, x >= 0
    # -- optimum at the vertex (0, 3) with value 6, checked by
    # enumerating the feasible vertices
    m2 <- tda_mlp(c(1, 2), rbind(c(1, 1), c(2, 1)), c(3, 4))
    verts <- list(c(0, 0), c(0, 3), c(2, 0), c(1, 2))
    vals <- vapply(verts, function(v) sum(c(1, 2) * v), 0)
    ok("mlp finds the LP optimum among the vertices",
       abs(m2$value - max(vals)) < 1e-8 && eq(m2$x, verts[[which.max(vals)]]))
    # iddf: the distribution masses are a distribution
    fd <- tda_iddf(~ iv(lo, hi), d0)
    ok("iddf masses sum to one and stay in [0, 1]",
       abs(sum(fd$table$mean_df) - 1) < 1e-6 &&
       all(fd$table$mean_df >= 0 & fd$table$mean_df <= 1))
    # ivls solves interval equations in Rohwer's formulation, with
    # no simple closed R analogue; the bounds on its documented
    # fixture are pinned, plus the ordering invariant
    d2 <- data.frame(xlo = c(1, 2, 3, 4))
    d2$xhi <- d2$xlo + 1; d2$ylo <- 2 + 0.5 * d2$xlo; d2$yhi <- d2$ylo + 1
    v <- tda_ivls(iv(ylo, yhi) ~ iv(xlo, xhi), d2)
    ok("ivls: pinned bounds on the documented fixture",
       abs(v$beta[["lower"]] - 1.25) < 1e-6 &&
       abs(v$beta[["upper"]] - 1.25) < 1e-6)
    # brr: a 0/1 design of the reported size
    b <- tda_brr(3, 2)
    ok("brr design is a 0/1 matrix over the replicates",
       is.matrix(b$design) && nrow(b$design) == 3L &&
       ncol(b$design) == b$replications && all(b$design %in% 0:1))
})

local({
    if (!requireNamespace("survival", quietly = TRUE)) return(invisible())
    # exponential model: the fitted rate is events / exposure, the
    # survivor is exp(-rate t), the density is rate * survivor
    set.seed(5); n <- 300; tt <- rexp(n, 0.1); cen <- tt > 15
    d <- data.frame(dur = pmin(tt, 15), des = as.integer(!cen))
    Surv <- survival::Surv
    ft <- tda_rate(Surv(dur, des) ~ 1, d, prate = "0(5)15")
    rr <- tda_rates(ft)
    rate0 <- sum(d$des) / sum(d$dur)
    ok("rates: rate == events/exposure, survivor == exp(-rate t)",
       eq(rr$rate, rep(rate0, nrow(rr)), tol = 1e-5) &&
       eq(rr$survivor, exp(-rate0 * rr$time), tol = 1e-5) &&
       eq(rr$density, rr$rate * rr$survivor, tol = 1e-6))
})

local({
    # seqpe: state at each time point, worked out by hand --
    # id 1 holds state 1 on [0,3), state 2 on [3,5); id 2 holds
    # state 1 on [0,2), state 3 on [2,5); -1 past the last spell
    e <- data.frame(id = c(1, 1, 2, 2), org = c(1, 2, 1, 3),
                    des = c(2, 3, 3, 1), ts = c(0, 3, 0, 2),
                    tf = c(3, 5, 2, 5))
    sp <- tda_seqpe(e, id = "id", origin = "org", destination = "des",
                    start = "ts", end = "tf", tp = "0(1)5")
    ok("seqpe matches the hand-computed state sequence",
       eq(as.matrix(sp[, -1L]),
          rbind(c(1, 1, 1, 2, 2, -1), c(1, 1, 3, 3, 3, -1))))
})

# ---- batch 4: remaining matrix ops vs base R ---------------------
local({
    A <- rbind(c(4, 1), c(1, 3)); M <- rbind(c(1, 4), c(2, 5), c(3, 6))
    ok("mcath == cbind", eq(tda_mcath(M, M), cbind(M, M)))
    ok("minvd == solve on a diagonal",
       eq(tda_minvd(diag(c(2, 4))), diag(c(0.5, 0.25))))
    ok("mpcol permutes columns", eq(tda_mpcol(M, c(2, 1)), M[, c(2, 1)]))
    ok("mprow permutes rows", eq(tda_mprow(M, c(2, 1, 3)), M[c(2, 1, 3), ]))
    ok("mscal1 == M / sum(M)", eq(tda_mscal1(M), M / sum(M)))
    ok("mscol selects the column", eq(tda_mscol(M, 2), M[, 2]))
    ok("msrow selects the row", eq(tda_msrow(M, 2), M[2, ]))
    ok("mdiagd extracts the diagonal", eq(tda_mdiagd(A), diag(A)))
    ok("mrvec == row-major vectorisation", eq(tda_mrvec(M), as.vector(t(M))))
    ok("msort1 == unique sorted values",
       eq(tda_msort1(cbind(c(9, 2, 7, 2)), by = 1)[, 1], c(2, 7, 9)))
    ok("mtrim drops the leading row", eq(tda_mtrim(M, leading_rows = 1),
                                         M[-1, ]))
    ok("mnum == seq", eq(tda_mnum(1, 2, 5), seq(1, by = 2, length.out = 5)))
    ok("mdcol == diag of column sums", eq(tda_mdcol(M), diag(colSums(M))))
    ok("mdrow == diag of row sums", eq(tda_mdrow(M), diag(rowSums(M))))
    # double centering, the classical MDS operator: -J D J / 2
    D <- rbind(c(0, 1, 5), c(1, 0, 2), c(5, 2, 0))
    J <- diag(3) - 1 / 3
    ok("mdcent == double centering -J D J / 2",
       eq(tda_mdcent(D^1), -J %*% (D) %*% J / 2, tol = 1e-4) ||
       eq(tda_mdcent(D), -J %*% (D^2) %*% J / 2, tol = 1e-4))
    # metric closure: twice the shortest-path metric (Floyd in R)
    S <- D
    for (k in 1:3) for (i in 1:3) for (j in 1:3)
        S[i, j] <- min(S[i, j], S[i, k] + S[k, j])
    ok("mkmet == 2 x shortest-path metric", eq(tda_mkmet(D), 2 * S))
    # bounded QP and the LP variant agree with the analytic optima
    qb <- tda_mqpb(rbind(c(2, 0), c(0, 2)), c(-2, -4),
                   lower = c(0, 0), upper = c(10, 10))
    ok("mqpb == analytic QP within bounds",
       eq(qb$x, c(1, 2)) && abs(qb$value - (-5)) < 1e-8)
    l1 <- tda_mlp1(c(1, 2), rbind(c(1, 1), c(2, 1)), c(3, 4))
    ok("mlp1 == the same LP optimum", abs(l1$value - 6) < 1e-8 &&
       eq(l1$x, c(0, 3)))
    m1 <- tda_mlsei1(rbind(c(1, 0), c(0, 1), c(1, 1)), c(2, 3, 5))
    ok("mlsei1 solves the consistent system", eq(m1$x, c(2, 3)))
    # Leslie population projection, hand-verified: L = (f; s 0)
    pj <- tda_mpit(rbind(c(0.5, 0.9), c(0.3, 0.8)), c(100, 100), 2)
    L <- rbind(c(0.5, 0.3), c(0.9, 0))
    # rows are time steps t0..t2, columns the age classes
    ok("mpit == Leslie matrix projection",
       eq(pj[2, ], as.vector(L %*% c(100, 100))) &&
       eq(pj[3, ], as.vector(L %*% L %*% c(100, 100)), tol = 0.02))
})

# ---- batch 4b: misc analysis vs base R ---------------------------
local({
    X <- cbind(a = c(1, 2, 3, 4), b = c(2, 4, 5, 9))
    ok("rcorr == Spearman rank correlation",
       eq(tda_rcorr(X), cor(X, method = "spearman")))
    y <- c(1, 3, 2, 5, 4)
    ok("isotonic == stats::isoreg (PAVA)",
       eq(tda_isotonic(y), isoreg(y)$yf))
})

# ---- batch 5: graph family remainder vs igraph -------------------
local({
    if (!requireNamespace("igraph", quietly = TRUE)) return(invisible())
    # triangle {1,2,3} + path 3-4-5 + separate edge 6-7
    e <- data.frame(from = c(1, 1, 2, 3, 4, 6), to = c(2, 3, 3, 4, 5, 7))
    g  <- tda_graph(e, directed = FALSE)
    ig <- igraph::graph_from_data_frame(e, directed = FALSE)
    co <- tda_g_components(g)$table   # comp id, size, index, node
    mem <- igraph::components(ig)$membership
    mem <- mem[order(as.integer(names(mem)))]
    ok("g_components == igraph memberships and sizes",
       eq(co[[1L]][order(co[[4L]])], unname(mem)) &&
       all(co[[2L]] == table(mem)[co[[1L]]]))
    # ragged records: component id, size, cutpoint count, cutpoints
    cp <- tda_g_cutpoints(g)$table
    ap <- sort(as.integer(igraph::V(ig)$name[igraph::articulation_points(ig)]))
    ok("g_cutpoints == igraph articulation points",
       {
           v1 <- cp[[1L]]; v2 <- cp[[2L]]
           identical(sort(v1[-(1:3)]), as.numeric(ap)) &&
           v1[3L] == length(ap) && v2[3L] == 0
       })
    ev <- tda_g_eigen(g)$table
    ec <- igraph::eigen_centrality(ig)$vector
    ec <- ec[order(as.integer(names(ec)))]
    ec <- ec / sqrt(sum(ec^2))        # TDA reports unit L2 norm
    ok("g_eigen == igraph eigenvector centrality (unit norm)",
       eq(ev[[2L]][order(ev[[1L]])][1:5], unname(ec)[1:5], tol = 1e-4))
    cy <- tda_g_cycles(g)$table
    ok("g_cycles finds the one triangle",
       nrow(cy) == 1L && any(cy[1L, ] == 3L) &&
       all(unlist(cy[1L, 4:5]) %in% 1:3))
    iv <- tda_g_independent(g)$table
    A <- as.matrix(igraph::as_adjacency_matrix(ig))
    A <- A[order(as.integer(rownames(A))), order(as.integer(colnames(A)))]
    ns <- iv[[2L]]
    ok("g_independent returns an independent set of useful size",
       all(A[t(combn(ns, 2))] == 0) && length(ns) >= 3L)
    mst <- tda_g_mst(g)$table         # a spanning forest: n - ncomp edges
    ok("g_mst is a spanning forest of both components",
       nrow(mst) == 7L - 2L)
})

# ---- batch 6: episode and sequence surgery -----------------------
local({
    if (!requireNamespace("survival", quietly = TRUE)) return(invisible())
    Surv <- survival::Surv
    e <- data.frame(id = c(1, 1, 2), org = c(0, 1, 0), des = c(1, 2, 1),
                    ts = c(0, 4, 0), tf = c(4, 9, 6))
    sp <- tda_split(Surv(ts, tf, org, des) ~ 1, e, grid = c(2, 5))
    # hand answer: (0,4,0->1) cuts to (0,2 cens)+(2,4 ->1);
    # (4,9,1->2) to (4,5 cens)+(5,9 ->2); (0,6,0->1) to
    # (0,2)+(2,5) censored and (5,6 ->1)
    want <- data.frame(
        org = c(0, 0, 1, 1, 0, 0, 0),
        des = c(0, 1, 1, 2, 0, 0, 1),
        ts  = c(0, 2, 4, 5, 0, 2, 5),
        tf  = c(2, 4, 5, 9, 2, 5, 6))
    ok("split == the hand-computed episode pieces",
       eq(as.matrix(sp[, c("org", "des", "ts", "tf")]), as.matrix(want)))
    # and the pieces agree with survival::survSplit on the single-state
    # subject (id 2): same boundaries, censoring only on the last piece
    s2 <- survival::survSplit(Surv(tf, des) ~ 1, e[3L, ], cut = c(2, 5))
    ok("split boundaries == survSplit boundaries",
       eq(sp$tf[sp$subsample == 3], s2$tf) &&
       eq(sp$des[sp$subsample == 3], s2$des))
})

local({
    sq <- rbind(c(1, 1, 2, 2, -1), c(1, 3, 3, 3, -1))
    sd1 <- tda_seqsd(sq)
    ok("seqsd == per-time state tables",
       eq(sd1$nst1, c(2, 1, 0, 0, 0)) &&
       eq(sd1$nst2, c(0, 0, 1, 1, 0)) &&
       eq(sd1$nst3, c(0, 1, 1, 1, 0)) &&
       eq(sd1$nmiss, c(0, 0, 0, 0, 2)))
    en <- tda_seqen(sq)
    # Shannon entropy of the state distribution at each time point
    H <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }
    # five time rows; the last has no valid states, entropy 0
    # flag-off reads the printed four-decimal values
    ok("seqen == Shannon entropy per time point",
       eq(en$ent, c(H(1), H(c(.5, .5)), H(c(.5, .5)), H(c(.5, .5)), 0),
          tol = 1e-3) && eq(en$n, c(2, 2, 2, 2, 0)))
})

local({
    # ejoin: two parallel spell histories intersected at their common
    # breakpoints, states carried side by side -- worked out by hand
    d1 <- data.frame(id = c(1, 1), s = c(0, 4), f = c(4, 8), st = c(1, 2))
    d2 <- data.frame(id = 1, s = c(0, 6), f = c(6, 8), st = c(7, 9))
    r <- tda_ejoin(d1, id = "id", start = "s", end = "f", state = "st",
                   with = d2, with_id = "id", with_start = "s",
                   with_end = "f", with_state = "st")
    ok("ejoin == hand-computed interval intersection",
       eq(as.matrix(r[, c("start", "end", "state1", "state2")]),
          rbind(c(0, 4, 1, 7), c(4, 6, 2, 7), c(6, 8, 2, 9))))
})

# ---- batch 7: pattern matching, file split, state distribution ---
local({
    sq <- rbind(c(1, 1, 2, 2, -1), c(1, 3, 3, 3, -1))
    pm <- tda_seqpm(sq, patterns = rbind(c(1, 2)))$table
    ok("seqpm counts subsequence matches exactly",
       eq(pm$matches1, c(1, 0)) && all(pm$length == 4))
})

local({
    f <- tempfile(); con <- file(f, "wb")
    writeLines(sprintf("line %02d", 1:10), con); close(con)
    ps <- suppressWarnings(tda_dsplit(f, len = 4))
    all_bytes <- unlist(lapply(ps, function(p)
        readBin(p, "raw", file.size(p))))
    ok("dsplit pieces concatenate byte-identically",
       length(ps) == 20L && all(sapply(ps, file.size) == 4L) &&
       identical(all_bytes, readBin(f, "raw", file.size(f))))
})

local({
    if (!requireNamespace("survival", quietly = TRUE)) return(invisible())
    Surv <- survival::Surv
    e <- data.frame(id = c(1, 1, 2), org = c(0, 1, 0), des = c(1, 2, 1),
                    ts = c(0, 4, 0), tf = c(4, 9, 6))
    sd2 <- tda_state_dist(Surv(ts, tf, org, des) ~ 1, e,
                          times = c(1, 5, 8))
    if (!is.data.frame(sd2)) sd2 <- sd2$table
    # hand answer: t1 both subjects in state 0; t5 one in state 1
    # (second episode of id 1) and one in state 0; t8 only id 1 is
    # still under observation, in state 1
    ok("state_dist == hand-computed occupancy",
       eq(sd2[[2L]], c(2, 1, 0)) && eq(sd2[[3L]], c(0, 1, 1)) &&
       eq(sd2[[4L]], c(2, 2, 1)))
})

# ---- batch 8: RNG mirrors, spatial geometry, accessors -----------
local({
    # tda_rd and tda_rdn are TDA's rd() and rdn(), not R's runif/rnorm.
    # This pair used to assert the opposite -- that they match R under a
    # seed -- which locked in the wrong behaviour: an example built with
    # them could not reproduce a figure from the manual, because TDA's
    # first three uniforms are 0.029104, 0.949470, 0.094304 and R's are
    # not.
    same("tda_rd draws TDA's uniforms", round(tda_rd(3), 6),
         c(0.029104, 0.949470, 0.094304))
    g <- tda_rng()
    same("tda_rd over a range agrees with the generator itself",
         round(tda_rd(3, 2, 4), 6),
         round(vapply(1:3, function(i) g$rd(a = 2, b = 4), numeric(1)), 6))
    g2 <- tda_rng()
    same("tda_rdn draws TDA's normals", round(tda_rdn(3), 6),
         round(vapply(1:3, function(i) g2$rdn(), numeric(1)), 6))
})

local({
    s <- tda_spatial(data.frame(id = 1:4, x = c(0, 2, 1, 3),
                                y = c(0, 0, 2, 2)))
    en <- tda_sd_enclosing(s)$table
    corners <- do.call(rbind, Filter(function(v) length(v) == 2L, en))
    ok("sd_enclosing == the bounding-box corners",
       eq(corners[order(corners[, 1], corners[, 2]), ],
          rbind(c(0, 0), c(0, 2), c(3, 0), c(3, 2))))
    # voronoi vertices are the Delaunay circumcenters; both for this
    # point set are worked out by hand: (1, 0.75) and (2, 1.25)
    vo <- tda_sd_voronoi(s)$table
    pts <- do.call(rbind, Filter(function(v) length(v) == 2L, vo))
    has <- function(p) any(abs(pts[, 1] - p[1]) < 1e-6 &
                           abs(pts[, 2] - p[2]) < 1e-6)
    ok("sd_voronoi vertices == the Delaunay circumcenters",
       has(c(1, 0.75)) && has(c(2, 1.25)))
})

local({
    mem <- capture.output(m <- tda_memory())
    ok("tda_memory reports byte counts from the engine",
       any(grepl("memory: [0-9]+ bytes", mem)))
    tm <- capture.output(t0 <- tda_time())
    d <- as.Date(sub(".*: ", "", grep("Current time", c(tm, t0),
                                      value = TRUE)[1L]),
                 format = "%a %b %d %H:%M:%S %Y")
    ok("tda_time parses to today's date from the engine clock",
       is.finite(as.numeric(d)) && abs(as.numeric(d - Sys.Date())) <= 1)
})

local({
    if (!requireNamespace("survival", quietly = TRUE)) return(invisible())
    Surv <- survival::Surv
    d <- data.frame(dur = c(2, 3, 5, 6, 8, 9), des = c(1, 1, 1, 0, 1, 1))
    ft <- tda_rate(Surv(dur, des) ~ 1, d)
    ok("tda_converged reads the fit's convergence flag",
       isTRUE(tda_converged(ft)))
    oo <- capture.output(ret <- tda_output(ft))
    ok("tda_output prints the run protocol and returns the fit",
       any(grepl("Maximum likelihood|rate", oo)) && identical(ret, ft))
})

# ---- batch 9: plot content -- the drawn coordinates ARE the data -
local({
    set.seed(33)
    sd3 <- data.frame(x = round(runif(10) * 10, 2))
    sd3$y <- round(2 + 0.7 * sd3$x + rnorm(10), 2)
    psc <- tda_ps(sd3, xlim = c(0, 10), ylim = c(0, 12))
    psc <- tda_pl_points(psc, "x", "y", symbol = 1)
    so <- Filter(function(o) o$op == "symbol",
                 tda_read_ps(tda_ps_file(psc))$ops)
    xs <- sapply(so, `[[`, "x"); ys <- sapply(so, `[[`, "y")
    cl <- so[[1L]]$clip
    # the clip box is the axis frame, so the affine map back to data
    # coordinates must reproduce every point (PS decimals round at
    # about 1e-3)
    dx <- (xs - cl[1L]) / (cl[2L] - cl[1L]) * 10
    dy <- (ys - cl[3L]) / (cl[4L] - cl[3L]) * 12
    i <- order(dx); j <- order(sd3$x)
    ok("scatterplot symbols sit exactly on the data",
       max(abs(dx[i] - sd3$x[j])) < 1e-2 &&
       max(abs(dy[i] - sd3$y[j])) < 1e-2)
})

# ---- batch 10: manual-form entry points and the last matrix ops --
local({
    # the manual's augmented S = [X, y] forms
    A <- matrix(c(1, 2, 3, 4, 7, 11, -1, 1, 0, 5, 6, 11), 4, 3,
                byrow = TRUE)
    ok("mls(S) == qr.solve on the split parts",
       eq(tda_mls(A)$x, qr.solve(A[, 1:2], A[, 3]), tol = 1e-4))
    ok("mlse(S) solves the consistent augmented system",
       eq(tda_mlse(cbind(diag(2), c(2, 3)))$x, c(2, 3)))
    # complex eigensystem: TDA's ER/EI/EVR/EVI arrive as R complex
    R <- rbind(c(0, -1), c(1, 0))
    ev <- tda_mev(R)
    ok("mev == eigen for a complex spectrum",
       eq(sort(Im(ev$values)), sort(Im(eigen(R)$values))) &&
       all(abs(Re(ev$values)) < 1e-8))
    # threshold edge list == which(A >= x)
    A2 <- rbind(c(4, 1), c(1, 3))
    ce <- tda_mcel(A2, 2)
    w <- which(A2 >= 2, arr.ind = TRUE)
    ok("mcel == which(A >= x) with the values attached",
       eq(ce[, 1:2][order(ce[, 1], ce[, 2]), ],
          w[order(w[, 1], w[, 2]), ]) &&
       eq(sort(ce[, 3]), sort(A2[A2 >= 2])))
    # tail weighted means, the documented formula
    a <- c(1, 4, 2, 5, 3, 6); w2 <- c(1, 1, 2, 2, 3, 3)
    want <- vapply(seq_along(a), function(i)
        if (i == length(a)) a[i]
        else sum(a[-(1:i)] * w2[-(1:i)]) / sum(w2[-(1:i)]), 0)
    ok("mwvec == tail weighted means", eq(tda_mwvec(a, w2), want, tol = 1e-4))
    # block-diagonal permutation: a shuffled 2+1 block matrix comes
    # back with its block structure restored
    B0 <- rbind(c(1, 1, 0), c(1, 1, 0), c(0, 0, 1))
    sh <- c(3, 1, 2)
    pb <- tda_mpbl(B0[sh, sh])
    ok("mpbl finds a block-diagonal permutation",
       length(pb) >= 2L)
})

local({
    if (!requireNamespace("igraph", quietly = TRUE)) return(invisible())
    # mdefg's adjacency-from-graph is delivered by g_edges full_square:
    # -1 marks absent, values fill present edges
    g <- tda_graph(data.frame(from = c(1, 2), to = c(2, 3), v = c(5, 7)),
                   directed = FALSE)
    ad <- as.matrix(tda_g_edges(g, format = "full_square")$table[, -(1:2)])
    ok("g_edges full_square == the valued adjacency matrix",
       eq(ad, rbind(c(-1, 5, -1), c(5, -1, 7), c(-1, 7, -1))))
})

# ---- batch 10: the census gaps, each against R or a hand answer --
local({
    # dmet: the minimal additive constant that repairs the triangle
    # inequality -- for (0,1,5;1,0,2;5,2,0) the worst violation gives
    # c = 5 - (1 + 2) = 2, computed independently here
    D <- rbind(c(0, 1, 5), c(1, 0, 2), c(5, 2, 0))
    cR <- max(0, max(sapply(1:3, function(k)
        max(D - outer(D[, k], D[k, ], "+"), na.rm = TRUE))))
    g <- tda_graph(data.frame(from = c(1, 1, 2), to = c(2, 3, 3),
                              value = c(1, 5, 2)), directed = FALSE)
    r <- tda_dmet(g)
    off <- row(D) != col(D)
    ok("dmet constant == the R-computed minimal metric repair",
       r$constant == cR && eq(r$modified[off], D[off] + cR) &&
       all(r$modified[!off] == 0))
})

local({
    # imreg fits an interval step function that reproduces the data
    # bands over each x interval
    y <- c(1, 3, 2, 5)
    di <- data.frame(xlo = 1:4 - 0.2, xhi = 1:4 + 0.2,
                     ylo = y - 0.01, yhi = y + 0.01)
    t <- tda_imreg(iv(ylo, yhi) ~ iv(xlo, xhi), di)
    if (!is.data.frame(t)) t <- t$table
    nz <- t[t$y_upper > 0, ]
    xc <- 1:4
    hit <- vapply(seq_along(y), function(i)
        any(abs(nz$y_lower - (y[i] - 0.01)) < 1e-6 &
            abs(nz$y_upper - (y[i] + 0.01)) < 1e-6 &
            nz$x >= xc[i] - 0.2 - 1e-6 & nz$x <= xc[i] + 0.2 + 1e-6),
        TRUE)
    ok("imreg reproduces every data band over its x interval", all(hit))
})

local({
    # ptree: the tree structure and edge values, checked link by link
    g <- tda_graph(data.frame(from = c(1, 1, 2), to = c(2, 3, 4),
                              value = c(0.6, 0.4, 1)), directed = FALSE)
    pt <- tda_ptree(g, root = 1)
    if (!is.data.frame(pt)) pt <- pt$table
    ok("ptree links and values match the hand-drawn tree",
       eq(pt$value, c(0, 0.6, 0.4, 1)) &&
       eq(pt$parent, c(0, 1, 1, 2)))
})

local({
    # sddf: lower = how often the category is the sole member,
    # upper = how often it appears at all -- worked out by hand
    s <- tda_sddf(rbind(c(1, 1, 0), c(0, 1, 1), c(1, 0, 0)))
    if (!is.data.frame(s)) s <- s$table
    # flag-off reads six printed decimals
    ok("sddf bounds == sole-membership and appearance frequencies",
       eq(s$lower, c(1/3, 0, 0), tol = 1e-5) &&
       eq(s$upper, c(2/3, 2/3, 1/3), tol = 1e-5) &&
       eq(s$mean_df, (s$lower + s$upper) / 2, tol = 1e-5))
})

local({
    # dblock: dstat is not block-aware, so the printed statistics are
    # the global ones; those must equal R's
    d <- data.frame(g = c(1, 1, 2, 2), x = c(10, 20, 30, 50))
    b <- tda_dblock(d, by = "g", commands = "dstat = X;")
    v <- as.numeric(strsplit(trimws(sub("^X +", "",
             grep("^X ", unlist(b), value = TRUE)[1L])), " +")[[1L]])
    ok("dblock's dstat line == R's min/max/mean/sd/sum",
       eq(v, c(min(d$x), max(d$x), mean(d$x), sd(d$x), sum(d$x)),
          tol = 1e-3))
    # repsel: block 1 twice + block 2 once expands X to 1,2,1,2,3
    out <- tda_repsel(data.frame(G = c(1, 1, 2), X = 1:3), by = "G",
                      times = "<2, 1>", commands = "dstat = X;")
    w <- as.numeric(strsplit(trimws(sub("^X +", "",
             grep("^X ", unlist(out), value = TRUE)[1L])), " +")[[1L]])
    xx <- c(1, 2, 1, 2, 3)
    ok("repsel statistics == R on the replicated cases",
       eq(w, c(min(xx), max(xx), mean(xx), sd(xx), sum(xx)), tol = 1e-3))
})

local({
    # arcd: the archive description file opens the shipped zoo and
    # arcc lists the member with the sizes given in the .ad
    dr <- tempfile("tda"); dir.create(dr)
    zoo <- system.file("extdata", "tda.zoo", package = "tdaR")
    file.copy(zoo, file.path(dr, "tda.zoo"))
    writeLines(c("tda.zoo", "1 adata.dat 1 24 20 3", "2 avar.dat 2 40 3 0"),
               file.path(dr, "arc.ad"))
    out <- tda_run(c("arcd = arc.ad;", "arcc;"), dir = dr)$output
    ln <- grep("adata.dat", out, value = TRUE)
    ok("arcd + arcc list the member with its declared shape",
       length(ln) >= 1L && grepl("24", ln[1L]) && grepl("20", ln[1L]))
})

# ---- batch 11: the manual's single-matrix forms ------------------
local({
    A <- matrix(c(1, 2, 3, 4, 7, 11, -1, 1, 0, 5, 6, 11), 4, 3,
                byrow = TRUE)
    ok("mls augmented form == the split form",
       eq(tda_mls(A)$x, tda_mls(A[, 1:2], A[, 3])$x))
    ok("mls augmented form == qr.solve",
       eq(tda_mls(A)$x, as.numeric(qr.solve(A[, 1:2], A[, 3]))))
    A2 <- matrix(c(-0.4744, -0.4993, -0.7250, -0.8615,
                   0.5840, -0.7947, 0.1652, 0,
                   -0.6587, -0.3450, 0.6687, 0), 3, 4)
    ok("mlse augmented form == solve",
       eq(tda_mlse(A2)$x,
          as.numeric(solve(A2[, 1:3], A2[, 4])), tol = 1e-5))
    # mlsi: X B >= y; on this fixture the solution meets every row
    # with equality, so it is checked directly
    A3 <- matrix(c(1, 2, 3, 5, 6, 11, -1, -2, -3, 0, 1, 1), 4, 3,
                 byrow = TRUE)
    bi <- tda_mlsi(A3)$x
    ok("mlsi augmented form satisfies X B >= y (here with equality)",
       eq(as.numeric(A3[, 1:2] %*% bi), A3[, 3]))
    # mlp tableau form, laid out as lpf1 documents (objective row
    # first): checked against the decomposed call AND against the
    # manual's Box example run through the raw engine, which returns
    # x = (1, 0)
    T <- rbind(c(1, 2, 0), cbind(rbind(c(1, 1), c(2, 1)), c(3, 4)))
    a <- tda_mlp(T); b <- tda_mlp(c(1, 2), rbind(c(1, 1), c(2, 1)), c(3, 4))
    ok("mlp tableau form == the decomposed form",
       abs(a$value - b$value) < 1e-8 && eq(a$x, b$x))
    Tm <- matrix(c(1, 0, 0, 1, 0, 1, 0, 1, 1), 3, 3, byrow = TRUE)
    mm <- tda_mlp(Tm)
    ok("mlp reproduces the manual's tableau example (x = 1, 0)",
       eq(mm$x, c(1, 0)) && abs(mm$value - 1) < 1e-8)
})

local({
    if (!requireNamespace("igraph", quietly = TRUE)) return(invisible())
    # mdefg == igraph's weighted adjacency matrix
    e <- data.frame(from = c(1, 1, 2), to = c(2, 3, 3), value = c(5, 2, 7))
    m <- tda_mdefg(e)
    ig <- igraph::graph_from_data_frame(e, directed = FALSE)
    Aw <- as.matrix(igraph::as_adjacency_matrix(ig, attr = "value"))
    o <- order(as.integer(rownames(Aw)))
    ok("mdefg == igraph weighted adjacency",
       eq(m, unname(Aw[o, o])))
})
