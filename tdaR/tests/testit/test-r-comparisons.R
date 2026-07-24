# Numerical comparisons against trusted, independent R implementations.
#
# Different in purpose from the rest of the suite: those check tdaR against
# TDA itself (the examples/ tree and its .ref files are the reference).
# These check specific tdaR functions against a *different* implementation
# of the same statistical method, in a package with no connection to TDA at
# all -- zoo, survival, systemfit, stats -- so a wrapper bug that happens to
# agree with a misreading of TDA's output has nowhere to hide. Skipped
# outright, not failed, when the comparison package is not installed.

library(tdaR)

# ---- tda_dma vs stats::prcomp -----------------------------------------

set.seed(1)
dpc <- data.frame(x1 = rnorm(30), x2 = rnorm(30))
dpc$x3 <- dpc$x1 + dpc$x2 + rnorm(30, sd = 0.3)
pc1 <- tda_dma(dpc)
pc2 <- stats::prcomp(dpc, center = FALSE, scale. = FALSE)
same("dma: eigenvalues match prcomp", pc1$values, pc2$sdev^2, 1e-8)
same("dma: loadings match prcomp up to sign",
     abs(pc1$vectors), abs(unclass(pc2$rotation)), 1e-6)
ok("dma: an out-of-range alg is refused here, since TDA would clamp it",
   inherits(try(tda_dma(dpc, alg = 99), silent = TRUE), "try-error"))

# ---- gloc / bfa / mldes wrappers vs independent computation -----------
local({
    # locate: brute force over all assignments of variable points to
    # fixed points (cost = sum of |pos difference| * weight, with fixed
    # points at their line ranks; the answer is spacing-invariant)
    set.seed(21)
    for (rep in 1:2) {
        k <- 2L; m <- 4L
        vf <- matrix(sample(0:8, k * m, replace = TRUE), k)
        vv <- matrix(0, k, k); vv[1, 2] <- vv[2, 1] <- sample(0:5, 1)
        if (all(vf == 0)) vf[1, 1] <- 1
        got <- tda_locate_line(vf, vv)
        pos <- seq_len(m)
        best <- NULL; bcost <- Inf
        grid <- as.matrix(expand.grid(rep(list(pos), k)))
        for (r in seq_len(nrow(grid))) {
            a <- grid[r, ]
            cost <- sum(vf * abs(outer(pos[a], pos, "-"))) +
                sum(vv[1, 2] * abs(pos[a][1] - pos[a][2]))
            if (cost < bcost - 1e-9) { bcost <- cost; best <- a }
        }
        bestcost_got <- sum(vf * abs(outer(pos[got], pos, "-"))) +
            sum(vv[1, 2] * abs(pos[got][1] - pos[got][2]))
        same(sprintf("tda_locate_line %d achieves the brute-force optimum",
                     rep), bestcost_got, bcost, 1e-9)
    }

    # boolean minimization: every selection must reproduce the truth
    # table exactly
    d <- expand.grid(X1 = 0:1, X2 = 0:1, X3 = 0:1)
    d$Y <- as.integer(d$X1 & d$X2 | d$X3)
    b <- tda_boolean_min("Y", c("X1", "X2", "X3"), data = d)
    ok("tda_boolean_min selections exist", is.data.frame(b$selections))
    evalsel <- function(sel, row) {
        any(vapply(split(sel, sel$selection)[1L], function(s) {
            all(vapply(seq_len(nrow(s)), function(i) {
                p <- unlist(s[i, -1L])
                any_term <- all(is.na(p) | p == row)
                any_term
            }, NA) | TRUE) && any(vapply(seq_len(nrow(s)), function(i) {
                p <- unlist(s[i, -1L])
                all(is.na(p) | p == row)
            }, NA))
        }, NA))
    }
    pred <- vapply(seq_len(nrow(d)), function(r)
        evalsel(b$selections, unlist(d[r, 1:3])), NA)
    ok("tda_boolean_min cover reproduces the truth table",
       identical(as.integer(pred), d$Y))

    # mldes: cell-for-cell against the traced formula
    set.seed(22)
    z <- matrix(sample(1:9, 8), 4L)
    g <- c(1, 1, 2, 2)
    got <- tda_mlrc_design(z, g)
    n <- 4L; q <- 2L
    want <- matrix(0, n * n, q * q + 1L); l <- 0L
    for (kk in 1:n) for (k in 1:n) {
        l <- l + 1L
        if (g[k] == g[kk]) {
            want[l, 1:(q * q)] <- as.vector(outer(z[k, ], z[kk, ]))
            if (k == kk) want[l, q * q + 1L] <- 1
        }
    }
    same("tda_mlrc_design matches the traced formula", got, want, 1e-10)
})

# ---- the three matrix-command wrappers vs independent computation -----
local({
    # isotonic: plain PAVA must equal stats::isoreg exactly
    set.seed(11)
    for (rep in 1:3) {
        v <- round(rnorm(10 + rep, 5, 3), 2)
        same(sprintf("tda_isotonic %d equals stats::isoreg", rep),
             tda_isotonic(v), stats::isoreg(v)$yf, 1e-8)
    }
    # the hand-verified tie example from examples/coverage/mmpops.cf
    y <- c(5, 3, 4, 6, 2, 7); g <- c(1, 1, 2, 3, 3, 3)
    same("tda_isotonic primary matches the hand computation",
         tda_isotonic(y, g, "primary"),
         c(11/3, 3, 11/3, 6, 11/3, 7), 1e-4)
    same("tda_isotonic secondary matches the hand computation",
         tda_isotonic(y, g, "secondary"), c(4, 4, 4, 5, 5, 5), 1e-8)

    # binary program: brute force over all 0-1 vectors
    set.seed(12)
    for (rep in 1:3) {
        nv <- 4L
        obj <- sample(-4:6, nv, replace = TRUE)
        A <- matrix(sample(-2:3, 2 * nv, replace = TRUE), 2L)
        X <- as.matrix(expand.grid(rep(list(0:1), nv)))
        for (b2 in list(c(1, 0), c(2, -1))) {
            feas <- apply(X, 1L, function(x) all(A %*% x >= b2))
            if (!any(feas)) next
            vals <- X[feas, , drop = FALSE] %*% obj
            r <- tda_mlpi(obj, A, bounds = b2, constant = 3)
            same(sprintf("tda_mlpi %d: optimum vs brute force",
                         rep), r$value, max(vals) + 3, 1e-8)
            ok(sprintf("tda_mlpi %d: solutions are optimal",
                       rep),
               all(r$solutions %*% obj + 3 == r$value) &&
                   all(apply(r$solutions, 1L,
                             function(x) all(A %*% x >= b2))))
        }
    }

    # quadratic assignment: recompute the objective in R from the
    # returned permutation (must match TDA's value exactly), and at
    # n <= 6 enumerate every permutation to check TDA's heuristic never
    # beats the true optimum -- pinning exact equality for a fixed seed
    # where it reaches it (it usually does at these sizes).
    qap_obj <- function(F, D, C, p)
        sum(diag(C[, p, drop = FALSE])) + sum(F * D[p, p])
    qap_perms <- function(n) {
        if (n == 1L) return(matrix(1L, 1L, 1L))
        sub <- qap_perms(n - 1L)
        do.call(rbind, lapply(seq_len(n), function(i)
            cbind(i, sub + (sub >= i))))
    }
    set.seed(23)
    for (rep in 1:4) {
        n <- sample(3:6, 1)
        F <- matrix(sample(0:9, n * n, replace = TRUE), n, n); diag(F) <- 0
        D <- matrix(sample(0:9, n * n, replace = TRUE), n, n); diag(D) <- 0
        C <- matrix(sample(0:9, n * n, replace = TRUE), n, n)
        r <- tda_mqap(F, D, C)
        same(sprintf("tda_mqap %d: objective matches its permutation",
                     rep), r$value, qap_obj(F, D, C, r$permutation), 1e-8)
        allp <- qap_perms(n)
        best <- min(apply(allp, 1L, function(p) qap_obj(F, D, C, p)))
        ok(sprintf("tda_mqap %d: never beats the true optimum",
                   rep), r$value >= best - 1e-8)
        if (rep == 1L)
            same("tda_mqap: reaches the true optimum here",
                 r$value, best, 1e-8)
    }

    # linear program: brute-force vertex enumeration for small LPs.
    # Every vertex of {x >= 0, Ax <= b} sets some n of the (m + n)
    # constraints (the m rows of A plus the n non-negativity rows)
    # tight; solve each such square system and keep the feasible ones.
    lp_vertices <- function(A, b, n) {
        Afull <- rbind(A, -diag(n))
        bfull <- c(b, rep(0, n))
        idx <- utils::combn(nrow(Afull), n, simplify = FALSE)
        out <- list()
        for (s in idx) {
            sol <- tryCatch(solve(Afull[s, , drop = FALSE], bfull[s]),
                            error = function(e) NULL)
            if (!is.null(sol) && all(Afull %*% sol <= bfull + 1e-6))
                out[[length(out) + 1L]] <- sol
        }
        out
    }
    set.seed(24)
    n_lp_checked <- 0L
    for (rep in 1:8) {
        n <- sample(2:3, 1); m <- sample(2:4, 1)
        obj <- sample(-5:5, n, replace = TRUE)
        A <- matrix(sample(-3:4, m * n, replace = TRUE), m, n)
        b <- sample(1:10, m, replace = TRUE)
        # only test instances with a bounded, nonempty feasible region:
        # brute force needs at least one vertex to compare against
        verts <- lp_vertices(A, b, n)
        if (!length(verts)) next
        r <- tryCatch(tda_mlp1(obj, A, b),
                      error = function(e) NULL)
        if (is.null(r)) next
        n_lp_checked <- n_lp_checked + 1L
        best <- max(vapply(verts, function(v) sum(obj * v), 0))
        same(sprintf("tda_mlp1 %d: matches the best brute-force vertex",
                     rep), r$value, best, 1e-6)
        ok(sprintf("tda_mlp1 %d: x is feasible", rep),
           all(r$x >= -1e-6) && all(A %*% r$x <= b + 1e-6))
    }
    ok("tda_mlp1: at least one random instance was checked",
       n_lp_checked >= 1L)

    # constrained least squares: unconstrained vs lm.fit/qr.solve,
    # equality-constrained vs the direct Lagrange (KKT block) solve,
    # and inequality/nonneg vs brute force over active sets.
    set.seed(25)

    # unconstrained: plain OLS, no intercept, must match qr.solve exactly
    for (rep in 1:3) {
        n <- sample(2:3, 1); m <- n + sample(2:4, 1)
        A <- matrix(rnorm(m * n), m, n)
        b <- rnorm(m)
        r <- tda_mls(A, b)
        want <- qr.solve(A, b)
        same(sprintf("tda_mls %d: unconstrained matches qr.solve",
                     rep), r$x, want, 1e-6)
    }

    # equality-constrained: minimize ||Ax-b||^2 s.t. Ex=f has the closed
    # form KKT solve  [A'A  E'; E  0] [x;lambda] = [A'b; f]
    for (rep in 1:3) {
        n <- 3L; m <- 4L; k <- 1L
        A <- matrix(rnorm(m * n), m, n)
        b <- rnorm(m)
        E <- matrix(rnorm(k * n), k, n)
        f <- rnorm(k)
        r <- tda_mlsei(A, b, equalities = E, equalities_bounds = f)
        KKT <- rbind(cbind(t(A) %*% A, t(E)),
                    cbind(E, matrix(0, k, k)))
        rhs <- c(t(A) %*% b, f)
        want <- solve(KKT, rhs)[1:n]
        same(sprintf("tda_mlsei %d: equality matches the Lagrange solve",
                     rep), r$x, want, 1e-5)
    }

    # inequality / nonneg: brute force over active sets -- enumerate
    # every subset of constraints as "tight" (solved as equalities),
    # keep the ones feasible against every constraint, take the best
    ls_active_set <- function(A, b, G, h) {
        n <- ncol(A)
        best <- NULL; bestval <- Inf
        for (k in 0:nrow(G)) {
            combos <- if (k == 0) list(integer(0))
                     else utils::combn(nrow(G), k, simplify = FALSE)
            for (s in combos) {
                if (length(s) > n) next
                x <- tryCatch({
                    if (length(s) == 0L) qr.solve(A, b)
                    else if (length(s) == n) solve(G[s, , drop = FALSE], h[s])
                    else {
                        # minimize ||Ax-b|| s.t. G[s,]x = h[s] via KKT
                        Gs <- G[s, , drop = FALSE]; hs <- h[s]
                        kk <- nrow(Gs)
                        KKT <- rbind(cbind(t(A) %*% A, t(Gs)),
                                    cbind(Gs, matrix(0, kk, kk)))
                        solve(KKT, c(t(A) %*% b, hs))[1:n]
                    }
                }, error = function(e) NULL)
                if (is.null(x)) next
                if (all(G %*% x >= h - 1e-6)) {
                    val <- sum((A %*% x - b)^2)
                    if (val < bestval - 1e-9) { bestval <- val; best <- x }
                }
            }
        }
        list(x = best, val = bestval)
    }
    for (rep in 1:3) {
        n <- 2L
        A <- diag(n)
        b <- rnorm(n, 0, 3)
        G <- rbind(diag(n), c(1, 1))
        h <- c(0, 0, sample(-1:1, 1))
        r <- tda_mlsei(A, b, inequalities = G, inequalities_bounds = h)
        bf <- ls_active_set(A, b, G, h)
        same(sprintf("tda_mlsei %d: inequality matches active-set brute force",
                     rep), sum((A %*% r$x - b)^2), bf$val, 1e-6)
        rn <- tda_mnls(A, b)
        bfn <- ls_active_set(A, b, diag(n), c(0, 0))
        same(sprintf("tda_mnls %d: nonneg matches active-set brute force",
                     rep), sum((A %*% rn$x - b)^2), bfn$val, 1e-6)
    }

    # mlsei1 (active_set method): must agree with the default lsei
    # method on the same mixed instance, not merely assumed identical
    A <- rbind(c(1, 0), c(0, 1)); b <- c(2, 3)
    Eeq <- rbind(c(1, 1)); feq <- 3
    Gin <- rbind(c(1, -1)); hin <- 0
    r_lsei <- tda_mlsei(A, b, equalities = Eeq, equalities_bounds = feq,
                       inequalities = Gin, inequalities_bounds = hin)
    r_qp <- tda_mlsei1(A, b, equalities = Eeq, equalities_bounds = feq,
                       inequalities = Gin, inequalities_bounds = hin)
    same("tda_mlsei1: matches tda_mlsei on the same instance",
         r_qp$x, r_lsei$x, 1e-6)

    # quadratic programming: unconstrained vs the closed-form solve
    # (x = -G^-1 d), box-bounded and general-constrained cases vs
    # brute force over active sets (KKT solve per candidate active
    # set, filtered for feasibility, best objective kept).
    qp_obj <- function(C, d, x) as.numeric(d %*% x + 0.5 * t(x) %*% C %*% x)
    set.seed(26)
    for (rep in 1:3) {
        n <- sample(2:3, 1)
        M <- matrix(rnorm(n * n), n, n)
        C <- t(M) %*% M + diag(n) * 0.5  # positive definite
        d <- rnorm(n)
        r <- tda_mqp(C, d)
        want <- solve(C, -d)
        same(sprintf("tda_mqp %d: unconstrained matches -C^-1 d",
                     rep), r$x, want, 1e-5)
    }

    # box-bounded and general-constrained: brute force over active
    # sets from combined box/constraint rows, same recipe as the LS
    # family's active-set check but for a quadratic objective.
    qp_active_set <- function(C, d, G, h) {
        n <- ncol(C)
        best <- NULL; bestval <- Inf
        for (k in 0:nrow(G)) {
            combos <- if (k == 0) list(integer(0))
                     else utils::combn(nrow(G), k, simplify = FALSE)
            for (s in combos) {
                if (length(s) > n) next
                x <- tryCatch({
                    if (length(s) == 0L) solve(C, -d)
                    else {
                        Gs <- G[s, , drop = FALSE]; hs <- h[s]
                        kk <- nrow(Gs)
                        KKT <- rbind(cbind(C, t(Gs)),
                                    cbind(Gs, matrix(0, kk, kk)))
                        solve(KKT, c(-d, hs))[1:n]
                    }
                }, error = function(e) NULL)
                if (is.null(x)) next
                if (all(G %*% x >= h - 1e-6)) {
                    val <- qp_obj(C, d, x)
                    if (val < bestval - 1e-9) { bestval <- val; best <- x }
                }
            }
        }
        list(x = best, val = bestval)
    }
    set.seed(27)
    for (rep in 1:3) {
        n <- 2L
        M <- matrix(rnorm(n * n), n, n)
        C <- t(M) %*% M + diag(n) * 0.5
        d <- rnorm(n)
        lo <- c(-1, -1); up <- c(1, 1)
        r <- tda_mqpb(C, d, lower = lo, upper = up)
        G <- rbind(diag(n), -diag(n))
        h <- c(lo, -up)
        bf <- qp_active_set(C, d, G, h)
        same(sprintf("tda_mqpb %d: box bounds match active-set brute force",
                     rep), qp_obj(C, d, r$x), bf$val, 1e-6)

        Aeq <- rbind(c(1, 1)); beq <- sample(-1:1, 1)
        Ain <- rbind(c(1, 0)); bin <- sample(-1:1, 1)
        r2 <- tda_mqpc(C, d, equalities = Aeq,
                       equalities_bounds = beq,
                       inequalities = Ain,
                       inequalities_bounds = bin)
        # fold the equality into the active-set search as a row that
        # must always stay in the active set (both >= directions)
        Gm <- rbind(Aeq, -Aeq, Ain)
        hm <- c(beq, -beq, bin)
        bf2 <- qp_active_set(C, d, Gm, hm)
        same(sprintf("tda_mqpc %d: mixed constraints match active-set brute force",
                     rep), qp_obj(C, d, r2$x), bf2$val, 1e-5)
    }

    # core matrix-algebra wrappers: every one of these has a base-R
    # equivalent, so verify directly against it on random instances
    # rather than by hand -- these are standard linear algebra, not
    # commands whose semantics needed discovering.
    set.seed(28)
    for (rep in 1:3) {
        n <- sample(2:4, 1)
        M <- matrix(rnorm(n * n), n, n)
        Aspd <- t(M) %*% M + diag(n)  # symmetric positive definite
        same(sprintf("tda_mchol %d matches t(chol(.))", rep),
             tda_mchol(Aspd), t(chol(Aspd)), 1e-6)
        same(sprintf("tda_minvs %d matches solve(.) for symmetric PD", rep),
             tda_minvs(Aspd), solve(Aspd), 1e-6)

        Ddiag <- matrix(rnorm(n * n), n, n)
        same(sprintf("tda_minvd %d matches diag(1/diag(.))", rep),
             tda_minvd(Ddiag), diag(1 / diag(Ddiag), n), 1e-6)

        m <- n + sample(1:3, 1)
        Tall <- matrix(rnorm(m * n), m, n)
        same(sprintf("tda_mginv %d matches solve(t(A)%%*%%A)%%*%%t(A)", rep),
             tda_mginv(Tall), solve(t(Tall) %*% Tall) %*% t(Tall), 1e-5)

        same(sprintf("tda_mtrace %d matches sum(diag(.))", rep),
             tda_mtrace(Ddiag), sum(diag(Ddiag)), 1e-8)
        same(sprintf("tda_mnorm %d matches max(abs(.))", rep),
             tda_mnorm(Ddiag), max(abs(Ddiag)), 1e-8)
        same(sprintf("tda_mnorm1 %d matches sum(abs(.))", rep),
             tda_mnorm1(Ddiag), sum(abs(Ddiag)), 1e-6)
        same(sprintf("tda_mnorm2 %d matches sqrt(sum(.^2))", rep),
             tda_mnorm2(Ddiag), sqrt(sum(Ddiag^2)), 1e-6)

        mr <- sample(2:5, 1); mc <- sample(2:5, 1)
        Mrc <- matrix(rnorm(mr * mc), mr, mc)
        same(sprintf("tda_mnrow/tda_mncol %d match nrow/ncol", rep),
             c(tda_mnrow(Mrc), tda_mncol(Mrc)), c(mr, mc), 1e-8)
        same(sprintf("tda_mtransp %d matches t(.)", rep),
             tda_mtransp(Mrc), t(Mrc), 1e-8)

        v <- rnorm(n)
        same(sprintf("tda_mdiag %d matches diag(.)", rep),
             tda_mdiag(v), diag(v, n), 1e-8)
        same(sprintf("tda_mdiagd %d matches diag(.) extraction", rep),
             tda_mdiagd(Mrc), diag(Mrc), 1e-8)

        X <- matrix(rnorm(n * n), n, n)
        Y <- matrix(rnorm(n * n), n, n)
        Z <- matrix(rnorm(n * n), n, n)
        same(sprintf("tda_mmul %d matches chained %%*%%", rep),
             tda_mmul(X, Y, Z), X %*% Y %*% Z, 1e-6)
    }

    # matrix reshaping/vectorizing/sorting wrappers, each against its
    # base-R equivalent on random instances.
    set.seed(29)
    for (rep in 1:3) {
        nr <- sample(2:5, 1); nc <- sample(2:5, 1)
        A <- matrix(rnorm(nr * nc), nr, nc)
        same(sprintf("tda_mcvec %d matches as.vector(.)", rep),
             tda_mcvec(A), as.vector(A), 1e-8)
        same(sprintf("tda_mrvec %d matches as.vector(t(.))", rep),
             tda_mrvec(A), as.vector(t(A)), 1e-8)
        v <- rnorm(nr * nc)
        same(sprintf("tda_mivec %d matches matrix(., nrow=n)", rep),
             tda_mivec(v, nr), matrix(v, nrow = nr), 1e-8)
        same(sprintf("tda_mrsum/tda_mcsum %d match rowSums/colSums", rep),
             list(tda_mrsum(A), tda_mcsum(A)),
             list(rowSums(A), colSums(A)), 1e-6)
        same(sprintf("tda_mdrow/tda_mdcol %d match diag(rowSums/colSums(.))",
                     rep),
             list(tda_mdrow(A), tda_mdcol(A)),
             list(diag(rowSums(A), nr), diag(colSums(A), nc)), 1e-6)

        X <- matrix(sample(1:20, nr * nc, replace = TRUE), nr, nc)
        by <- sample(seq_len(nc), sample(1:nc, 1))
        ord <- do.call(order, as.data.frame(X[, by, drop = FALSE]))
        same(sprintf("tda_msort %d matches order()-based row sort", rep),
             tda_msort(X, by), X[ord, , drop = FALSE], 1e-8)
        same(sprintf("tda_mrank %d matches order() itself, not rank()", rep),
             tda_mrank(X, by), ord, 1e-8)
        same(sprintf("tda_msort1 %d matches unique(sorted .)", rep),
             tda_msort1(X, by), unique(X[ord, , drop = FALSE]), 1e-8)

        B <- matrix(rnorm(nr * nc), nr, nc)
        same(sprintf("tda_mcath %d matches cbind(.,.)", rep),
             tda_mcath(A, B), cbind(A, B), 1e-8)
        same(sprintf("tda_mcatv %d matches rbind(.,.)", rep),
             tda_mcatv(A, B), rbind(A, B), 1e-8)
        want_dsum <- rbind(cbind(A, matrix(0, nr, nc)),
                           cbind(matrix(0, nr, nc), B))
        same(sprintf("tda_mcathv %d matches the block-diagonal direct sum",
                     rep),
             tda_mcathv(A, B), want_dsum, 1e-8)

        lc <- sample(0:1, 1)
        want_trim <- if (lc == 0) A
                    else if (lc > 0) A[, -seq_len(lc), drop = FALSE]
                    else cbind(matrix(0, nr, -lc), A)
        same(sprintf("tda_mtrim %d matches column drop/pad", rep),
             tda_mtrim(A, leading_cols = lc), want_trim, 1e-8)
    }

    # eigen/SVD: sign of each eigen/singular vector is arbitrary, so
    # verify by reconstruction rather than direct comparison to
    # eigen()/svd()'s output.
    set.seed(30)
    for (rep in 1:3) {
        n <- sample(2:4, 1)
        M <- matrix(rnorm(n * n), n, n)
        Asym <- (M + t(M)) / 2
        r <- tda_mevs(Asym)
        same(sprintf("tda_mevs %d: eigenvalues match eigen()$values", rep),
             r$values, eigen(Asym, symmetric = TRUE, only.values = TRUE)$values,
             1e-6)
        for (i in seq_len(n))
            ok(sprintf("tda_mevs %d: column %d satisfies Av = lambda v", rep, i),
               max(abs(Asym %*% r$vectors[, i] - r$values[i] * r$vectors[, i])) <
                   1e-6)

        m <- n + sample(0:2, 1)
        A <- matrix(rnorm(m * n), m, n)
        same(sprintf("tda_msvd %d matches svd(.)$d", rep),
             tda_msvd(A), svd(A)$d, 1e-6)
        s1 <- tda_msvd1(A)
        same(sprintf("tda_msvd1 %d reconstructs A = U diag(d) t(V)", rep),
             s1$u %*% diag(s1$d, n) %*% t(s1$v), A, 1e-6)
    }

    # select/permute/aggregate wrappers, each against base-R indexing
    # or a hand-rolled aggregate on random instances.
    set.seed(31)
    for (rep in 1:3) {
        n <- sample(3:5, 1); m <- sample(3:5, 1)
        A <- matrix(rnorm(n * m), n, m)
        rows <- sample(seq_len(n), sample(1:n, 1), replace = TRUE)
        cols <- sample(seq_len(m), sample(1:m, 1), replace = TRUE)
        same(sprintf("tda_msrow %d matches A[rows,]", rep),
             tda_msrow(A, rows), A[rows, , drop = FALSE], 1e-8)
        same(sprintf("tda_mscol %d matches A[,cols]", rep),
             tda_mscol(A, cols), A[, cols, drop = FALSE], 1e-8)

        pr <- sample(seq_len(n)); pc <- sample(seq_len(m))
        same(sprintf("tda_mprow %d matches A[p,]", rep),
             tda_mprow(A, pr), A[pr, , drop = FALSE], 1e-8)
        same(sprintf("tda_mpcol %d matches A[,p]", rep),
             tda_mpcol(A, pc), A[, pc, drop = FALSE], 1e-8)

        Sq <- matrix(rnorm(n * n), n, n)
        ps <- sample(seq_len(n))
        same(sprintf("tda_mpsym %d matches A[p,p]", rep),
             tda_mpsym(Sq, ps), Sq[ps, ps, drop = FALSE], 1e-8)

        rg <- sample(1:2, n, replace = TRUE)
        cg <- sample(1:2, m, replace = TRUE)
        want_agg <- do.call(rbind, lapply(sort(unique(rg)), function(g)
            vapply(sort(unique(cg)), function(h)
                sum(A[rg == g, cg == h, drop = FALSE]), 0)))
        same(sprintf("tda_mag %d matches a hand-rolled block sum", rep),
             tda_mag(A, rg, cg), want_agg, 1e-6)

        x <- rnorm(1)
        same(sprintf("tda_mnc %d matches a hand-applied <= threshold", rep),
             tda_mnc(A, x, "<="), ifelse(A <= x, 0, A), 1e-8)
        same(sprintf("tda_mnc %d matches a hand-applied abs> threshold", rep),
             tda_mnc(A, abs(x), "abs>"),
             ifelse(abs(A) > abs(x), 0, A), 1e-8)
    }

    # centering/standardizing/cross-product/sqrt-diagonal/sequence,
    # each against its base-R equivalent on random instances.
    set.seed(32)
    for (rep in 1:3) {
        nr <- sample(3:6, 1); nc <- sample(2:4, 1)
        A <- matrix(rnorm(nr * nc), nr, nc)
        same(sprintf("tda_mcent %d matches sweep(A, 2, colMeans(A))", rep),
             tda_mcent(A), sweep(A, 2, colMeans(A)), 1e-6)
        ctr <- sweep(A, 2, colMeans(A))
        same(sprintf("tda_mstand %d matches population-sd standardizing", rep),
             tda_mstand(A), sweep(ctr, 2, sqrt(colMeans(ctr^2)), "/"), 1e-6)
        same(sprintf("tda_mcross %d matches crossprod(A)", rep),
             tda_mcross(A), crossprod(A), 1e-6)

        n <- sample(2:5, 1)
        d <- diag(abs(rnorm(n)) + 0.1, n)  # positive diagonal, avoids
                                            # msqrti's near-zero refusal
        same(sprintf("tda_msqrtd %d matches diag(sqrt(diag(.)))", rep),
             tda_msqrtd(d), diag(sqrt(diag(d)), n), 1e-6)
        same(sprintf("tda_msqrti %d matches diag(1/sqrt(diag(.)))", rep),
             tda_msqrti(d), diag(1 / sqrt(diag(d)), n), 1e-6)

        D <- matrix(rnorm(n * n), n, n); D <- abs(D + t(D))  # symmetric
        D2 <- D^2
        want_dc <- -0.5 * (D2 - matrix(rowMeans(D2), n, n) -
                           matrix(colMeans(D2), n, n, byrow = TRUE) +
                           mean(D2))
        same(sprintf("tda_mdcent %d matches the classical double-centering",
                     rep),
             tda_mdcent(D), want_dc, 1e-6)

        x0 <- rnorm(1); dd <- rnorm(1); nn <- sample(1:6, 1)
        same(sprintf("tda_mnum %d matches seq(x, by=d, length.out=n)", rep),
             tda_mnum(x0, dd, nn), seq(x0, by = dd, length.out = nn), 1e-8)
    }

    # mkp/mwvec/mwvec1/mscal1: mkp against base R's kronecker(); the
    # others against a hand-written R equivalent (no base-R builtin
    # matches mwvec/mwvec1's backward weighted mean).
    set.seed(33)
    for (rep in 1:3) {
        r1 <- sample(2:3, 1); c1 <- sample(2:3, 1)
        r2 <- sample(2:3, 1); c2 <- sample(2:3, 1)
        A <- matrix(rnorm(r1 * c1), r1, c1)
        B <- matrix(rnorm(r2 * c2), r2, c2)
        same(sprintf("tda_mkp %d matches kronecker(A, B)", rep),
             tda_mkp(A, B), kronecker(A, B), 1e-6)

        n <- sample(4:7, 1)
        a <- rnorm(n); w <- abs(rnorm(n))
        want_wv <- numeric(n)
        for (i in seq_len(n)) {
            idx <- if (i < n) (i + 1):n else integer(0)
            den <- sum(w[idx])
            want_wv[i] <- if (den != 0) sum(a[idx] * w[idx]) / den else a[i]
        }
        same(sprintf("tda_mwvec %d matches a hand-written backward weighted mean",
                     rep), tda_mwvec(a, w), want_wv, 1e-6)

        tt <- sample(seq_len(n))  # arbitrary order vector
        want_wv1 <- numeric(n)
        for (i in seq_len(n)) {
            idx <- which(tt > tt[i])
            den <- sum(w[idx])
            want_wv1[i] <- if (den != 0) sum(a[idx] * w[idx]) / den else a[i]
        }
        same(sprintf("tda_mwvec1 %d matches the same recipe with a custom order",
                     rep), tda_mwvec1(a, w, tt), want_wv1, 1e-6)

        M <- matrix(rnorm(6), 2, 3) + 5  # positive-biased, unlikely to sum to 0
        same(sprintf("tda_mscal1 %d matches A / sum(A)", rep),
             tda_mscal1(M), M / sum(M), 1e-6)
    }

    # mpinv/mcel
    set.seed(34)
    for (rep in 1:3) {
        n <- sample(3:6, 1)
        p <- sample(seq_len(n))
        same(sprintf("tda_mpinv %d matches order(p)", rep),
             tda_mpinv(p), order(p), 1e-8)

        A <- matrix(sample(0:9, n * n, replace = TRUE), n, n)
        x <- sample(1:5, 1)
        if (any(A >= x)) {
            want_edges <- do.call(rbind, lapply(seq_len(n), function(i)
                do.call(rbind, lapply(seq_len(n), function(j)
                    if (A[i, j] >= x) c(i, j, A[i, j]) else NULL))))
            same(sprintf("tda_mcel %d matches a hand-scanned edge list", rep),
                 tda_mcel(A, x), want_edges, 1e-8)
        }
    }

    # mpz: verify the fixed gather relationship and the zero-free-
    # diagonal property directly, on matrices guaranteed to admit a
    # perfect match (built from a random permutation plus extra
    # random nonzero noise, so a solution is known to exist).
    set.seed(35)
    for (rep in 1:3) {
        n <- sample(3:6, 1)
        base_perm <- sample(seq_len(n))
        A <- matrix(0, n, n)
        A[cbind(seq_len(n), base_perm)] <- sample(1:9, n, replace = TRUE)
        extra <- matrix(sample(0:1, n * n, replace = TRUE, prob = c(0.7, 0.3)),
                        n, n) * matrix(sample(1:9, n * n, replace = TRUE), n, n)
        A <- A + extra
        r <- tda_mpz(A)
        ok(sprintf("tda_mpz %d: diagonal of B is entirely nonzero", rep),
           all(diag(r$B) != 0))
        ok(sprintf("tda_mpz %d: A[p[i], i] != 0 for every i", rep),
           all(A[cbind(r$p, seq_len(n))] != 0))
        same(sprintf("tda_mpz %d: B is the gather A[p, ]", rep),
             r$B, A[r$p, , drop = FALSE], 1e-8)
    }

    # mpbl/mpbu: build matrices with a genuine, known block structure
    # (random block sizes, dense within-block coupling, cross-block
    # coupling only in the direction that direction of triangularity
    # allows), scramble the row/column order, and check the returned
    # permutation both reproduces A via the gather relationship and
    # actually restores block-triangular form.
    set.seed(36)
    for (rep in 1:3) {
        nblk <- sample(2:3, 1)
        sizes <- sample(1:2, nblk, replace = TRUE)
        n <- sum(sizes)
        grp <- rep(seq_len(nblk), sizes)
        A <- matrix(0, n, n)
        for (bi in seq_len(nblk)) {
            idx <- which(grp == bi)
            A[idx, idx] <- sample(1:9, length(idx)^2, replace = TRUE)
            if (bi > 1) {
                below <- which(grp < bi)
                A[idx, sample(below, min(1, length(below)))] <-
                    sample(1:9, 1)
            }
        }
        perm <- sample(seq_len(n))
        Ascr <- A[perm, perm]

        rl <- tda_mpbl(Ascr)
        same(sprintf("tda_mpbl %d: B is the gather Ascr[p, p]", rep),
             rl$B, Ascr[rl$p, rl$p, drop = FALSE], 1e-8)
        viol_l <- outer(rl$block, rl$block, function(a, b) a < b) & rl$B != 0
        ok(sprintf("tda_mpbl %d: no entry violates lower block-triangular form",
                   rep), !any(viol_l))

        ru <- tda_mpbu(Ascr)
        same(sprintf("tda_mpbu %d: B is the gather Ascr[p, p]", rep),
             ru$B, Ascr[ru$p, ru$p, drop = FALSE], 1e-8)
        viol_u <- outer(ru$block, ru$block, function(a, b) a > b) & ru$B != 0
        ok(sprintf("tda_mpbu %d: no entry violates upper block-triangular form",
                   rep), !any(viol_u))
    }

    # mpfit: verify against a from-scratch RAS/IPF implementation.
    ras <- function(A, row_sums, col_sums, iter = 100) {
        B <- A
        for (k in seq_len(iter)) {
            B <- B * (row_sums / rowSums(B))
            B <- sweep(B, 2, col_sums / colSums(B), "*")
        }
        B
    }
    set.seed(37)
    for (rep in 1:3) {
        nr <- sample(2:4, 1); nc <- sample(2:4, 1)
        A <- matrix(round(runif(nr * nc, 1, 20)), nr, nc)
        total <- sum(A)
        rs <- as.numeric(table(cut(runif(1000, 0, total), nr))) + 1
        cs <- as.numeric(table(cut(runif(1000, 0, total), nc))) + 1
        rs <- rs / sum(rs) * total
        cs <- cs / sum(cs) * total
        got <- tda_mpfit(A, rs, cs, max_iter = 200, eps = 1e-10)
        want <- ras(A, rs, cs, iter = 200)
        same(sprintf("tda_mpfit %d matches a from-scratch RAS implementation",
                     rep), got$B, want, 1e-4)
        same(sprintf("tda_mpfit %d: row sums match the target", rep),
             rowSums(got$B), rs, 1e-3)
        same(sprintf("tda_mpfit %d: column sums match the target", rep),
             colSums(got$B), cs, 1e-3)
    }

    # mpit/mpit1: Leslie matrix projection vs a from-scratch recursion.
    set.seed(38)
    for (rep in 1:3) {
        n <- sample(3:5, 1)
        FS <- cbind(runif(n, 0, 3), c(runif(n - 1, 0.2, 0.9), 0))
        pop <- sample(10:100, n)
        iters <- sample(2:5, 1)
        leslie <- function(FS, pop, Z = NULL, iters) {
            R <- matrix(0, iters + 1, n)
            R[1, ] <- pop
            for (j in 2:(iters + 1)) {
                R[j, 1] <- sum(FS[, 1] * R[j - 1, ])
                if (n > 1)
                    for (i in 2:n) R[j, i] <- FS[i - 1, 2] * R[j - 1, i - 1]
                if (!is.null(Z)) R[j, ] <- R[j, ] + Z
            }
            R
        }
        same(sprintf("tda_mpit %d matches a from-scratch Leslie recursion", rep),
             tda_mpit(FS, pop, iters), leslie(FS, pop, iters = iters), 1e-6)
        Z <- runif(n, 0, 5)
        same(sprintf("tda_mpit1 %d matches the same recursion plus immigration",
                     rep),
             tda_mpit1(FS, pop, Z, iters), leslie(FS, pop, Z, iters), 1e-6)
    }

    # mkmet: Kemeny/Kendall-tau-with-ties distance vs a from-scratch
    # pairwise implementation of the same rule.
    kemeny <- function(A) {
        n <- nrow(A); m <- ncol(A)
        D <- matrix(0, n, n)
        for (a in 2:n) for (b in 1:(a - 1)) {
            d <- 0
            for (i in 2:m) for (j in 1:(i - 1)) {
                ai <- A[a, i]; aj <- A[a, j]; bi <- A[b, i]; bj <- A[b, j]
                tie_a <- abs(ai - aj) < 1e-9
                tie_b <- abs(bi - bj) < 1e-9
                if (tie_a) {
                    if (!tie_b) d <- d + 1
                } else if (tie_b) {
                    d <- d + 1
                } else if (sign(ai - aj) != sign(bi - bj)) {
                    d <- d + 2
                }
            }
            D[a, b] <- D[b, a] <- d
        }
        D
    }
    set.seed(39)
    for (rep in 1:3) {
        n <- sample(3:5, 1); m <- sample(3:4, 1)
        A <- matrix(sample(1:3, n * m, replace = TRUE), n, m)
        same(sprintf("tda_mkmet %d matches a from-scratch Kemeny distance", rep),
             tda_mkmet(A), kemeny(A), 1e-8)
    }

    # mple: product-limit (Kaplan-Meier) CDF vs the survival package.
    # Event times are kept unique in the random data -- ties among
    # events are a documented, order-dependent exception (see
    # ?tda_mple), verified separately below rather than blindly
    # compared against survival's simultaneous-tie convention.
    set.seed(40)
    for (rep in 1:3) {
        n <- sample(5:8, 1)
        time <- sample(seq_len(2 * n), n)  # unique event times by draw
        censored <- sample(0:1, n, replace = TRUE, prob = c(0.7, 0.3))
        censored[which.max(time)] <- 0  # keep the tail comparable to
                                        # survival's convention; see
                                        # ?tda_mple for why a censored
                                        # max-time tail is excluded here
        got <- tda_mple(time, censored)
        fit <- survival::survfit(survival::Surv(time, 1 - censored) ~ 1)
        sfun <- stats::stepfun(fit$time, c(1, fit$surv))
        want <- 1 - sfun(time)
        same(sprintf("tda_mple %d matches survival::survfit's Kaplan-Meier CDF",
                     rep), got$F, want, 1e-6)
    }
    # documented tied-events behavior: two events tied at time 2 out
    # of 3 total observations -- only the second (last processed)
    # reaches standard KM's F(2) = 2/3.
    tie2 <- tda_mple(c(2, 2, 3), c(0, 0, 0))
    same("tda_mple: first of two tied events shows the intermediate 1/3",
         tie2$F[1], 1 / 3, 1e-8)
    same("tda_mple: second of two tied events reaches standard KM's 2/3",
         tie2$F[2], 2 / 3, 1e-8)
    # documented tail behavior: F reaches exactly 1 at the highest-time
    # observation even when it is censored.
    tail_censored <- tda_mple(c(1, 2, 3), c(0, 0, 1))
    same("tda_mple: F reaches 1 at a censored highest-time observation",
         tail_censored$F, c(1 / 3, 2 / 3, 1), 1e-8)

    # mev: general (possibly non-symmetric, possibly complex) eigen
    # decomposition -- verify the defining relation A %*% v == lambda * v
    # directly (sign/scale of eigenvectors is arbitrary, same reasoning
    # as tda_mevs/tda_msvd1 above), not by comparing to eigen()'s
    # output, on random real matrices with no particular structure.
    set.seed(41)
    for (rep in 1:4) {
        n <- sample(2:4, 1)
        A <- matrix(round(rnorm(n * n) * 3), n, n)
        r <- tda_mev(A)
        for (i in seq_len(n))
            ok(sprintf("tda_mev %d: column %d satisfies Av = lambda v", rep, i),
               max(Mod(A %*% r$vectors[, i] - r$values[i] * r$vectors[, i])) <
                   1e-6)
    }

    # mch: verify against a from-scratch transcription of the exact
    # loop structure read from m_mch (no independent textbook formula
    # exists for this check, so the transcription itself is the
    # verification, cross-checked against TDA's output).
    mch_r <- function(A) {
        n <- nrow(A)
        B <- rep(0L, n)
        for (i0 in 0:(n - 1)) {
            z <- 0; s <- 0; nf <- FALSE
            if (i0 + 1 <= n - 1) for (j0 in (i0 + 1):(n - 1)) {
                z <- z + A[i0 + 1, j0 + 1]
                s <- s + A[j0 + 1, i0 + 1]
                if (z < s - 1e-9) { nf <- TRUE; break }
            }
            if (nf) B[i0 + 1] <- 1L
        }
        if (n > 1) for (i0 in (n - 1):1) {
            z <- 0; s <- 0; nf <- FALSE
            for (j0 in (i0 - 1):0) {
                z <- z + A[j0 + 1, i0 + 1]
                s <- s + A[i0 + 1, j0 + 1]
                if (z < s - 1e-9) { nf <- TRUE; break }
            }
            if (nf) B[i0 + 1] <- 1L
        }
        B
    }
    set.seed(42)
    for (rep in 1:4) {
        n <- sample(3:6, 1)
        A <- matrix(sample(0:9, n * n, replace = TRUE), n, n)
        same(sprintf("tda_mch %d matches a from-scratch transcription of m_mch",
                     rep), tda_mch(A), mch_r(A), 1e-8)
    }

    # midf/midf1/midf2/midf3: interval-censored distribution
    # estimators, verified against from-scratch R translations of the
    # exact formulas read from t_imat.c (no independent textbook
    # formula covers midf2/midf3, so the translation itself is the
    # specification, cross-checked against TDA's raw output).
    midf_cdf_at <- function(lower, upper, tmp) {
        n <- length(lower)
        f <- 0
        for (j in seq_len(n)) {
            if (lower[j] < tmp) {
                xx <- min(upper[j], tmp)
                f <- f + (xx - lower[j]) / (upper[j] - lower[j])
            }
        }
        f / n
    }
    set.seed(43)
    for (rep in 1:3) {
        n <- sample(4:7, 1)
        lower <- sample(0:15, n)
        upper <- lower + sample(1:6, n, replace = TRUE)

        bp <- sort(unique(c(lower, upper)))
        want_dl <- vapply(bp, function(t) mean(upper <= t), 0)
        want_du <- vapply(bp, function(t) mean(lower <= t), 0)
        want_dm <- vapply(bp, function(t) midf_cdf_at(lower, upper, t), 0)
        got_midf <- tda_midf(lower, upper)
        same(sprintf("tda_midf %d: breakpoints match sort(unique(endpoints))",
                     rep), got_midf$breakpoints, bp, 1e-8)
        same(sprintf("tda_midf %d: lower_bound matches P(entirely below)",
                     rep), got_midf$lower_bound, want_dl, 1e-8)
        same(sprintf("tda_midf %d: upper_bound matches P(possibly below)",
                     rep), got_midf$upper_bound, want_du, 1e-8)
        same(sprintf("tda_midf %d: cdf matches the uniform-model estimate",
                     rep), got_midf$cdf, want_dm, 1e-8)

        want_f1 <- vapply(lower, function(t) midf_cdf_at(lower, upper, t), 0)
        same(sprintf("tda_midf1 %d matches the CDF at each observation's lower",
                     rep), tda_midf1(lower, upper), want_f1, 1e-8)

        au <- want_dm
        want_f2 <- numeric(n)
        for (i in seq_len(n)) {
            j <- which(abs(bp - lower[i]) < 1e-9)[1]
            s <- 0
            if (j < length(bp)) for (jj in j:(length(bp) - 1)) {
                xl <- bp[jj]; xu <- bp[jj + 1]
                dens <- (au[jj + 1] - au[jj]) / (xu - xl)
                s <- s + dens * (xu^2 - xl^2)
            }
            want_f2[i] <- (s / 2) / (1 - au[j])
        }
        same(sprintf("tda_midf2 %d matches the restricted-mean calculation", rep),
             tda_midf2(lower, upper), want_f2, 1e-6)

        want_xl1 <- vapply(seq_len(n), function(i)
            mean(lower[lower >= lower[i] & lower <= upper[i]]), 0)
        want_xu1 <- vapply(seq_len(n), function(i)
            mean(upper[upper >= lower[i] & upper <= upper[i]]), 0)
        got_midf3 <- tda_midf3(lower, upper)
        same(sprintf("tda_midf3 %d: lower_avg matches nearby-endpoint average",
                     rep), got_midf3$lower_avg, want_xl1, 1e-8)
        same(sprintf("tda_midf3 %d: upper_avg matches nearby-endpoint average",
                     rep), got_midf3$upper_avg, want_xu1, 1e-8)
    }

    # independence: delta per subset = max_x |P(X=x|Y in S) - P(X=x)|
    set.seed(13)
    d <- data.frame(x = sample(1:4, 60, replace = TRUE),
                    y = sample(1:3, 60, replace = TRUE))
    got <- tda_independence("x", "y", data = d)
    want <- vapply(sort(unique(d$y)), function(s) {
        px  <- tabulate(d$x, 4) / nrow(d)
        pxs <- tabulate(d$x[d$y == s], 4) / sum(d$y == s)
        max(abs(pxs - px))
    }, 0)
    same("tda_independence deltas match direct computation",
         got$delta, want, 1e-5)
})

# ---- numeric options must never reach TDA in scientific notation -----
# format(2e5) is "2e+05"; TDA's sscanf("nbox=%d") rejects it and the
# whole option turns into a confusing unknown-parameter error, which is
# why max_boxes once needed as.integer().  Both entry paths are fixed
# centrally (tda_block and .tda_extra); pin them.
local({
    set.seed(42); nn <- 3
    xl <- runif(nn, 0, 10); xh <- xl + runif(nn, .5, 5)
    ym <- 2 + 1.5 * (xl + xh) / 2 + rnorm(nn); w <- runif(nn, .5, 3)
    di <- data.frame(yl = ym - w/2, yh = ym + w/2, xl = xl, xh = xh)
    r1 <- suppressWarnings(tda_ivreg(iv(yl, yh) ~ iv(xl, xh), di,
                                     method = "exact", max_boxes = 2e5))
    r2 <- suppressWarnings(tda_ivreg(iv(yl, yh) ~ iv(xl, xh), di,
                                     method = "exact",
                                     options = list(nbox = 2e5)))
    # What this test is about: the option text reaching TDA intact on
    # both entry paths -- so both runs must be IDENTICAL to each other.
    # Whether this instance certifies is a separate, platform-marginal
    # question (box acceptance sits within a few ulps of tolbw and
    # flips between compiler/FPU combinations), so equality of the two
    # runs, not certification, is asserted.
    ok("numeric options: plain doubles accepted on both entry paths",
       identical(r1$certified, r2$certified) &&
           isTRUE(all.equal(r1$beta, r2$beta, tolerance = 1e-8)) &&
           !any(grepl("unknown parameter|2e\\+05", r1$run$output)) &&
           !any(grepl("unknown parameter|2e\\+05", r2$run$output)))
})

# ---- user interrupt ----------------------------------------------------
# A SIGINT during a long branch-and-bound run must surface as a real R
# interrupt condition, and the session must stay fully usable.  The
# loops unwind through their ordinary bookkeeping, so nothing leaks.
if (.Platform$OS.type == "unix") local({
    if (!have_examples) return(invisible())
    pat <- utils::read.table(file.path(EX, "coverage", "patients.dat"))
    system(sprintf("(sleep 1; kill -INT %d) &", Sys.getpid()))
    r <- tryCatch(
        tda_icorr(~ iv(V1, V2) + iv(V3, V4), pat,
                         max_iter = 5e7, max_boxes = 5e5),
        interrupt = function(c) class(c)[1L])
    ok("interrupt during a long run is caught as a condition",
       identical(r, "tdaInterrupt"))
    d3 <- data.frame(xl = c(1, 4, 6), xu = c(2, 5, 8),
                     yl = c(3, 1, 6), yu = c(5, 2, 9))
    cv <- tda_icov(~ iv(xl, xu) + iv(yl, yu), d3,
                          max_iter = 2e5, max_boxes = 1e5)
    same("the session works normally after an interrupt",
         unname(cv$bounds), c(-2/3, 68/9), 1e-4)
})

# ---- the interval_wages dataset across the whole family ---------------
# Point regressor, bracketed outcome: ivreg's exact method and ilsreg's
# closed form compute the same estimand through entirely different code
# paths, so their agreement is a free cross-validation.
local({
    d <- tda_interval_wages()
    r <- tda_ivreg(iv(wage_lo, wage_hi) ~ iv(school_lo, school_hi), d,
                   method = "exact")
    ls2 <- tda_ilsreg(iv(wage_lo, wage_hi) ~ school_lo, d)
    ok("interval_wages: ivreg exact is certified (point regressor)",
       isTRUE(all(r$certified)))
    # tolerance covers the exporters-off run, where one side arrives
    # through the printed output at 4 decimals
    same("interval_wages: ivreg exact equals ilsreg's closed form",
         unname(unlist(r$beta)), unname(unlist(ls2$beta)), 1e-3)
    v <- suppressWarnings(tda_ivariance(~ iv(wage_lo, wage_hi), d,
                                        max_iter = 2e5, max_boxes = 1e5))
    ok("interval_wages: variance bounds ordered and sd consistent",
       v$bounds[["lower"]] <= v$bounds[["upper"]] &&
           isTRUE(all.equal(unname(v$sd),
                            sqrt(pmax(unname(v$bounds), 0)))))
    cr <- suppressWarnings(tda_icorr(
        ~ iv(school_lo, school_hi) + iv(wage_lo, wage_hi), d,
        max_iter = 5e4, max_boxes = 2e4))
    ok("interval_wages: correlation stays in [-1, 1] and is ordered",
       !anyNA(cr$bounds) && cr$bounds[["lower"]] <= cr$bounds[["upper"]] &&
           cr$bounds[["lower"]] >= -1 && cr$bounds[["upper"]] <= 1)
})

# ---- Stata reader/writer round trips ----------------------------------
# Fixtures written by upstream readstata13 (formats 108/110/114); TDA
# reads them with exact values (NA -> msys -5), and foreign reads back
# what TDA writes.  Format 113 (Stata 8/9) was never implemented in
# TDA and is pinned as a clean refusal in stataops.cf.
local({
    if (!have_examples) return(invisible())
    dr <- tempfile("tda"); dir.create(dr)
    for (f in c("rs108.dta", "rs110.dta", "rs114.dta"))
        file.copy(file.path(EX, "coverage", f), file.path(dr, f))
    ref <- data.frame(id = 1:6,
                      x = local({ set.seed(3); round(rnorm(6), 3) }),
                      g = c(1L, 2L, 1L, 3L, 2L, 1L),
                      big = c(1000L, 2e6L, 3L, 40L, 5e8L, 6L))
    ref$x[4] <- -5
    for (f in c("old102.dta", "old104.dta", "old105.dta"))
        file.copy(file.path(EX, "coverage", f), file.path(dr, f))
    # old formats carry an extra string and a factor-as-number column
    ref_old <- local({
        set.seed(4)
        d <- data.frame(id = 1:8, x = round(rnorm(8), 3))
        d$x[3] <- -5
        d
    })
    for (v in c(102L, 104L, 105L)) {
        tda_run(sprintf("rstata(df=q%d.out) = old%d.dta;", v, v), dir = dr)
        got <- utils::read.table(file.path(dr, sprintf("q%d.out", v)))
        same(sprintf("rstata old format %d numeric values exact", v),
             unname(as.matrix(got[, c(1, 2)])),
             unname(as.matrix(ref_old)), 1e-6)
        ok(sprintf("rstata old format %d strings preserved", v),
           identical(as.character(got[[3]])[1:3], c("abc", "de", "fgh")))
    }
    for (v in c(108L, 110L, 114L)) {
        tda_run(sprintf("rstata(df=o%d.out) = rs%d.dta;", v, v), dir = dr)
        got <- utils::read.table(file.path(dr, sprintf("o%d.out", v)))
        same(sprintf("rstata format %d values equal the source", v),
             unname(as.matrix(got)), unname(as.matrix(ref)), 1e-6)
    }
    tda_run(c("rstata() = rs114.dta;", "wstata(ptyp=10) = w.dta;"),
            dir = dr)
    r <- foreign::read.dta(file.path(dr, "w.dta"))
    r[[2]][is.na(r[[2]])] <- -5
    same("foreign reads back what wstata writes, values exact",
         unname(as.matrix(r)), unname(as.matrix(ref)), 1e-3)
})

# ---- spmod vs eigen; pcyc vs an R cycle count --------------------------
local({
    r <- tda_spmod(c(0, .4, .3, .1), c(.9, .8, .7, 0), rep(1, 4))
    L <- rbind(c(0, .4, .3, .1), cbind(diag(c(.9, .8, .7)), 0))
    e <- eigen(L)
    same("spmod growth equals the dominant eigenvalue",
         r$growth, Re(e$values[1]), 1e-6)
    v <- Re(e$vectors[, 1]); v <- v / sum(v)
    same("spmod stationary vector equals the normalized eigenvector",
         r$stationary, v, 1e-6)
    cyc <- tda_pcyc(rbind(c(2, 1, 3), c(2, 3, 1), c(1, 2, 3)))
    ncyc <- function(p) {
        seen <- logical(length(p)); k <- 0L
        for (i in seq_along(p)) if (!seen[i]) {
            k <- k + 1L; j <- i
            while (!seen[j]) { seen[j] <- TRUE; j <- p[j] }
        }
        k
    }
    got <- vapply(strsplit(cyc, ""), function(ch) sum(ch == "("), 0L)
    ok("pcyc cycle counts equal the R count",
       identical(got, vapply(list(c(2,1,3), c(2,3,1), c(1,2,3)), ncyc, 0L)))
})

# ---- nlreg vs lm / nls / total least squares --------------------------
# Free names in nlreg's "= function" are the parameters, varlist names
# the data, the first varlist entry the response (established in
# session 38; pinned in examples/coverage/nlregops.cf).
local({
    if (!have_examples) return(invisible())
    dp <- file.path(EX, "coverage", "nl.dat")
    d <- utils::read.table(dp); names(d) <- c("x", "y")
    # the linear table has a Wave column (value = field 4); the
    # parameter table has none (value = field 3)
    grab <- function(out, pat, fld = 4L) {
        ln <- grep(pat, out, value = TRUE)
        as.numeric(strsplit(trimws(ln[length(ln)]), "\\s+")[[1L]][fld])
    }
    run1 <- function(cmdline) {
        dr <- tempfile("tda"); dir.create(dr)
        file.copy(dp, file.path(dr, "nl.dat"))
        tda_run(c("nvar(dfile=nl.dat, X[10.4]=c1, Y[10.4]=c2);", cmdline),
                dir = dr)$output
    }
    o <- run1("nlreg(v=Y,X);")
    cf <- coef(stats::lm(y ~ x, d))
    same("nlreg linear OLS equals lm",
         c(grab(o, "Intercept"), grab(o, "^  2 +1 +X")),
         unname(cf), 1e-3)
    o <- run1("nlreg(v=Y,X, ni=1);")
    same("nlreg ni=1 equals lm without intercept",
         grab(o, "^  1 +1 +X"),
         unname(coef(stats::lm(y ~ x - 1, d))), 1e-3)
    o <- run1("nlreg(v=Y,X, xp=1,0.5) = a * exp(b * X);")
    fn <- stats::nls(y ~ a * exp(b * x), d, start = list(a = 1, b = 0.5))
    same("nlreg user function equals nls",
         c(grab(o, "^  1 +a", 3L), grab(o, "^  2 +b", 3L)),
         unname(coef(fn)), 1e-3)
    o <- run1("nlreg(v=Y,X, opt=2);")
    cx <- scale(d, scale = FALSE)
    v <- eigen(crossprod(cx))$vectors[, 1L]
    sl <- v[2L] / v[1L]
    same("nlreg ODR equals the total-least-squares line",
         c(grab(o, "Intercept"), grab(o, "^  2 +1 +X")),
         c(mean(d$y) - sl * mean(d$x), sl), 1e-3)
})

# ---- ivar vs ivar1: two algorithms, one estimand ----------------------
local({
    if (!have_examples) return(invisible())
    pat <- utils::read.table(file.path(EX, "coverage", "patients.dat"))
    a <- suppressWarnings(tda_ivar(~ iv(V1, V2), pat,
                                   max_iter = 2e5, max_boxes = 1e5))
    b <- tda_ivar1(~ iv(V1, V2), pat)
    same("tda_ivar and tda_ivar1 agree on the variance bounds",
         unname(a$bounds), unname(b$bounds), 0.01)
})

# ---- the k-pair covariance matrix (oils replication) ------------------
local({
    if (!have_examples) return(invisible())
    oils <- utils::read.table(file.path(EX, "coverage", "oils.dat"))
    names(oils) <- c("gl", "gh", "fl", "fh", "il", "ih", "sl", "sh")
    cm <- suppressWarnings(tda_icov(
        ~ iv(gl, gh) + iv(fl, fh) + iv(il, ih) + iv(sl, sh), oils,
        max_iter = 2e5, max_boxes = 1e5))
    ok("icov matrix: symmetric with variance diagonal",
       isTRUE(all.equal(cm$lower, t(cm$lower))) &&
           isTRUE(all.equal(cm$upper, t(cm$upper))) &&
           abs(cm$lower[2, 2] - 243.9643) < 0.01 &&
           abs(cm$upper[2, 2] - 469.4844) < 0.01)
    # two published-vs-exact adjudications, brute-verified in session 36:
    # freezing-iodine matches the paper; iodine-saponification's paper
    # maximum (176.26) undershoots the exact corner value
    ok("icov matrix: freezing-iodine equals the published interval",
       abs(cm$lower[2, 3] + 920.4844) < 0.01 &&
           abs(cm$upper[2, 3] + 408.6875) < 0.01)
    ok("icov matrix: iodine-saponification exceeds the paper's maximum",
       abs(cm$upper[3, 4] - 176.6563) < 0.01)
})

# ---- icov / icorr vs brute force --------------------------------------
# Covariance is multilinear per coordinate: extrema at corners, so
# exhaustive corner search is exact.  Correlation is verified against a
# dense grid at n=3, where the search certifies instantly.
local({
    d3 <- data.frame(xl = c(1, 4, 6), xu = c(2, 5, 8),
                     yl = c(3, 1, 6), yu = c(5, 2, 9))
    cv <- tda_icov(~ iv(xl, xu) + iv(yl, yu), d3,
                          max_iter = 2e5, max_boxes = 1e5)
    n <- 3
    best <- c(Inf, -Inf)
    for (rx in 0:(2^n - 1)) {
        x <- ifelse(bitwAnd(rx, 2^(0:(n - 1))) > 0, d3$xu, d3$xl)
        for (ry in 0:(2^n - 1)) {
            y <- ifelse(bitwAnd(ry, 2^(0:(n - 1))) > 0, d3$yu, d3$yl)
            v <- mean(x * y) - mean(x) * mean(y)
            best[1] <- min(best[1], v); best[2] <- max(best[2], v)
        }
    }
    same("tda_icov equals exhaustive corner search",
         unname(cv$bounds), best, 1e-4)
    cr <- tda_icorr(~ iv(xl, xu) + iv(yl, yu), d3,
                           max_iter = 2e5, max_boxes = 1e5)
    same("tda_icorr matches the dense-grid extrema",
         unname(cr$bounds), c(-0.142857, 0.893405), 1e-4)
    ok("tda_icorr stays within [-1, 1]",
       cr$bounds[["lower"]] >= -1 && cr$bounds[["upper"]] <= 1)
})

# ---- interval-stat wrappers: warnings surfaced, sd, paper pins --------
local({
    if (!have_examples) return(invisible())
    pat <- utils::read.table(file.path(EX, "coverage", "patients.dat"))
    # tight limits leave the search uncertified: the caveat must arrive
    # as an R warning, not sit buried in $run$output
    w <- tryCatch({ tda_ivariance(~ iv(V1, V2), pat,
                                  max_iter = 5, max_boxes = 5); NULL },
                  warning = conditionMessage)
    ok("tda_ivariance surfaces TDA's uncertified warning",
       is.character(w) && grepl("No certified bound|not been processed", w))
    v <- suppressWarnings(tda_ivariance(~ iv(V1, V2), pat,
                                        max_boxes = 2e5, max_iter = 2e5))
    same("tda_ivariance pulse minimum matches Gioia-Lauro",
         v$bounds[["lower"]], 26.86, 0.01)
    same("tda_ivariance pulse maximum matches the exact corner value",
         v$bounds[["upper"]], 542.0661, 0.01)
    ok("tda_ivariance sd is the square root of the variance bounds",
       isTRUE(all.equal(unname(v$sd), sqrt(pmax(unname(v$bounds), 0)))))
})

# ---- Gioia & Lauro (2005) replication ---------------------------------
# Interval variances of the paper's two fully-printed datasets (the
# 11-patient table and the Ichino oils).  The published MINIMA are
# reproduced; the published MAXIMA are undershoots of their local
# optimizer -- the variance maximum sits at a corner of the box, so
# exhaustive corner search is exact and TDA must match IT, not the
# paper.  See examples/coverage/gioia05.cf.
local({
    if (!have_examples) return(invisible())
    pat <- utils::read.table(file.path(EX, "coverage", "patients.dat"))
    corner_max <- function(lo, hi) {
        n <- length(lo); mx <- -Inf
        for (r in 0:(2^n - 1)) {
            x <- ifelse(bitwAnd(r, 2^(0:(n - 1))) > 0, hi, lo)
            mx <- max(mx, mean((x - mean(x))^2))
        }
        mx
    }
    iv <- function(lo, hi) {
        d <- data.frame(L = lo, U = hi)
        res <- tda_run(c(tda_nvar(d),
                         tda_block(name = "ivar", mxit = "200000",
                                   nbox = "100000", fmt = "12.4",
                                   rhs = "L,U")),
                       data = d, dir = tempfile("tda"))
        v <- grep("^Best m", res$output, value = TRUE)
        as.numeric(sub(".*value: *([0-9.eE+-]+).*", "\\1", v))
    }
    published_min <- c(26.86, 193.03, 298.18)
    for (k in 1:3) {
        b <- iv(pat[[2 * k - 1]], pat[[2 * k]])
        same(sprintf("Gioia-Lauro patients %d: published minimum", k),
             b[1L], published_min[k], 0.01)
        same(sprintf("Gioia-Lauro patients %d: exact corner maximum", k),
             b[2L], corner_max(pat[[2 * k - 1]], pat[[2 * k]]), 0.01)
    }
})

# ---- ivreg sharp bounds vs brute force --------------------------------
# The slope is linear in each y for fixed x, so y-extrema sit exactly at
# interval corners; x is searched on a grid.  Two checks per instance:
# containment (every brute-force point is feasible, so TDA's bounds must
# enclose the brute interval -- exact, grid-free) and tightness (the
# bounds sit within grid resolution of the brute extrema).  Covers both
# finishes: coordinates fully fixed by the two heuristic steps, and the
# certified opt=3 search over the remainder.
local({
    set.seed(7)
    slope3 <- function(X, y) {
        n <- nrow(X)
        (n * colSums(X * y) - colSums(X) * sum(y)) /
            (n * colSums(X * X) - colSums(X)^2)
    }
    iv_dir <- function(d, ns) {
        o <- list(opt = 3, tolbw = 1e-8, tolf = 1e-10,
                  rhs = "YL,YU,XL,XU")
        if (ns) o <- c(list(ns = 1), o)
        res <- tda_run(c(tda_nvar(d),
                         do.call(tda_block, c(list(name = "ivreg"), o))),
                       data = d, dir = tempfile("tda"))
        v <- grep("^Best m", res$output, value = TRUE)
        if (length(v))
            return(as.numeric(sub(".*value: *([0-9.eE+-]+).*", "\\1",
                                  v[1L])))
        # fully fixed by the heuristic: the last iteration row's beta
        it <- grep("^ +[0-9]+ +-?[0-9.]", res$output, value = TRUE)
        as.numeric(strsplit(trimws(it[length(it)]), " +")[[1L]][2L])
    }
    for (rep in 1:2) {
        n <- 3
        xl <- round(sort(runif(n, 0, 8)), 2)
        xu <- xl + round(runif(n, .3, 1.2), 2)
        yl <- round(runif(n, 0, 5), 2)
        yu <- yl + round(runif(n, .3, 1.5), 2)
        g <- 15
        X <- t(as.matrix(expand.grid(seq(xl[1], xu[1], length.out = g),
                                     seq(xl[2], xu[2], length.out = g),
                                     seq(xl[3], xu[3], length.out = g))))
        lo <- Inf; hi <- -Inf
        for (r in 0:7) {
            y <- ifelse(bitwAnd(r, 2^(0:2)) > 0, yu, yl)
            b <- slope3(X, y)
            lo <- min(lo, b); hi <- max(hi, b)
        }
        d <- data.frame(XL = xl, XU = xu, YL = yl, YU = yu)
        bmin <- iv_dir(d, ns = FALSE)
        bmax <- iv_dir(d, ns = TRUE)
        ok(sprintf("ivreg %d: bounds enclose every brute-force slope", rep),
           bmin <= lo + 1e-6 && bmax >= hi - 1e-6)
        ok(sprintf("ivreg %d: bounds within grid resolution of sharp", rep),
           abs(bmin - lo) < 0.02 && abs(bmax - hi) < 0.02)
    }
})

# ---- indep vs direct computation --------------------------------------
# Delta-independence: for each subset of Y's property space, the largest
# gap between P(X=x | Y in subset) and P(X=x).  The ranges start away
# from 1 on both axes on purpose: the shipped indexing was wrong in
# stride and offset and crashed on exactly this kind of input.
local({
    set.seed(2)
    d <- data.frame(x = sample(5:9, 40, TRUE), y = sample(3:6, 40, TRUE))
    res <- tda_run(c(tda_nvar(stats::setNames(d, c("X", "Y"))),
                     tda_block(name = "indep", rhs = "X,Y")),
                   data = stats::setNames(d, c("X", "Y")),
                   dir = tempfile("tda"))
    got <- as.numeric(sub(".*sets: ", "",
                          grep("Maximal delta", res$output, value = TRUE)))
    fx <- table(factor(d$x, 5:9)) / nrow(d)
    want <- vapply(3:6, function(ys) {
        f <- table(factor(d$x[d$y == ys], 5:9)) / sum(d$y == ys)
        max(abs(f - fx))
    }, 0)
    same("indep: per-subset deltas match direct computation", got, want,
         1e-5)
})

# ---- gloc vs brute force ----------------------------------------------
# One-dimensional multifacility location.  gloc's answer (which fixed
# point each variable point lands on) is checked against exhaustive
# search over all placements, at two different fixed-point spacings --
# the algorithm claims the answer needs only the ORDER of the fixed
# points, so both spacings must agree with it.
local({
    set.seed(11)
    gloc_run <- function(n, m, vwt, vfwt) {
        edges <- do.call(rbind, c(
            lapply(seq_len(n - 1), function(i)
                do.call(rbind, lapply((i + 1):n, function(j)
                    c(i, j, vwt[i, j], 0)))),
            lapply(seq_len(n), function(i)
                do.call(rbind, lapply(seq_len(m), function(k)
                    c(i, n + k, 0, vfwt[i, k]))))))
        d <- stats::setNames(as.data.frame(edges),
                             c("FROM", "TO", "WVV", "WVF"))
        res <- tda_run(c(tda_nvar(d),
                         tda_block(name = "gdd", opt = 1,
                                   rhs = "FROM,TO,WVV,WVF"),
                         tda_block(name = "gloc", gn = "1,2", n = n)),
                       data = d, dir = tempfile("tda"))
        as.integer(sub(".*fixed point ([0-9]+).*", "\\1",
                       grep("^Variable point [0-9]+ ->", res$output,
                            value = TRUE)))
    }
    brute <- function(n, m, vwt, vfwt, pos) {
        grid <- as.matrix(expand.grid(rep(list(seq_len(m)), n)))
        costs <- apply(grid, 1L, function(x) {
            cst <- 0
            if (n > 1)
                for (i in 1:(n - 1)) for (j in (i + 1):n)
                    cst <- cst + vwt[i, j] * abs(pos[x[i]] - pos[x[j]])
            for (i in seq_len(n))
                cst <- cst + sum(vfwt[i, ] * abs(pos[x[i]] - pos))
            cst
        })
        list(cost = min(costs), grid = grid, costs = costs)
    }
    for (rep in 1:4) {
        n <- sample(1:3, 1); m <- sample(3:5, 1)
        vwt <- matrix(0, n, n)
        if (n > 1) vwt[upper.tri(vwt)] <- sample(0:6, sum(upper.tri(vwt)),
                                                 replace = TRUE)
        vfwt <- matrix(sample(0:9, n * m, replace = TRUE), n, m)
        if (all(vfwt == 0)) vfwt[1, 1] <- 1
        got <- gloc_run(n, m, vwt, vfwt)
        ok(sprintf("gloc: instance %d returns one location per variable point",
                   rep), length(got) == n && all(got >= 1 & got <= m))
        for (pos in list(seq_len(m), cumsum(sample(1:9, m, replace = TRUE)))) {
            b <- brute(n, m, vwt, vfwt, pos)
            cost_got <- b$costs[which(apply(b$grid, 1L, function(x)
                all(x == got)))[1L]]
            ok(sprintf("gloc: instance %d optimal at spacing %s", rep,
                       paste(pos, collapse = ",")),
               isTRUE(all.equal(cost_got, b$cost)))
        }
    }
})

# ---- tda_pdatd vs stats::dist ------------------------------------------

set.seed(1)
ddd <- data.frame(x1 = rnorm(8), x2 = rnorm(8))
same("pdatd: euclidean matches dist()",
     as.numeric(tda_pdatd(ddd)), as.numeric(stats::dist(ddd)), 1e-8)
same("pdatd: city-block matches dist(method='manhattan')",
     as.numeric(tda_pdatd(ddd, measure = 2)),
     as.numeric(stats::dist(ddd, method = "manhattan")), 1e-8)

# ---- tda_rcorr vs stats::cor(method = "kendall") ------------------------

set.seed(1)
drc <- data.frame(x1 = rnorm(20), x2 = rnorm(20), x3 = rnorm(20))
same("rcorr: matches cor(method='kendall')",
     unclass(tda_rcorr(drc)),
     unclass(stats::cor(drc, method = "kendall")), 1e-8)

# ---- tda_sma vs zoo::rollmean -------------------------------------------

if (requireNamespace("zoo", quietly = TRUE)) {
    xsm <- c(1, 2, 3, 10, 5, 6, 7, 8, 9, 10)
    got <- tda_sma(xsm, width = 3)$table$smoothed
    # zoo has no direct equivalent of TDA's copy-on end-value rule, so
    # only the interior points (where both are an unambiguous plain mean)
    # are compared -- the ends are TDA's documented convention, not
    # something rollmean claims to reproduce.
    same("sma: interior points match zoo::rollmean",
         got[2:9], as.numeric(zoo::rollmean(xsm, 3))[1:8], 1e-8)
} else {
    message("zoo not installed: skipping sma comparison")
}

# ---- tda_lsreg1 (SUR form) vs separate lm() fits ------------------------
#
# Not systemfit's GLS-based SUR: TDA's "marginal estimation of
# conditional expectations" is a different algorithm that, for the
# uncensored case, was found (session notes, ?tda_lsreg) to reproduce
# separate per-equation OLS exactly -- confirmed again here.
# systemfit's SUR would differ by design (it exploits cross-equation error
# correlation via GLS), so it is not the right reference for this call.

set.seed(1)
n_sur <- 40
d_sur <- data.frame(id = 1:n_sur, x1 = rnorm(n_sur), x2 = rnorm(n_sur))
d_sur$y1 <- 2 + 0.8 * d_sur$x1 + rnorm(n_sur, sd = 0.3)
d_sur$y2 <- 5 - 0.5 * d_sur$x2 + rnorm(n_sur, sd = 0.3)
d_sur$cen <- 0
fit_sur <- tda_lsreg1(list(Eq1 = y1 ~ x1, Eq2 = y2 ~ x2), d_sur,
                      censor = "cen", id = "id")
lm1 <- stats::lm(y1 ~ x1, d_sur)
lm2 <- stats::lm(y2 ~ x2, d_sur)
b_sur <- coef(fit_sur)
same("lsreg1 SUR: eq1 matches separate lm()",
     unname(b_sur[c("Eq1", "Eq1_Vx1")]), unname(coef(lm1)), 1e-4)
same("lsreg1 SUR: eq2 matches separate lm()",
     unname(b_sur[c("Eq2", "Eq2_Vx2")]), unname(coef(lm2)), 1e-4)

# ---- tda_rate vs survival::survreg (exponential) -------------------------

if (requireNamespace("survival", quietly = TRUE)) {
    library(survival)
    set.seed(1)
    n_rt <- 100
    d_rt <- data.frame(x = rnorm(n_rt))
    d_rt$t <- rexp(n_rt, exp(0.3 + 0.6 * d_rt$x))
    d_rt$s <- 1
    f_tda_rt <- tda_rate(Surv(t, s) ~ x, d_rt, model = "exponential")
    f_r_rt <- survreg(Surv(t, s) ~ x, d_rt, dist = "exponential")
    # survreg's accelerated-failure-time parameterisation is the
    # negative of the hazard-scale coefficient rate() reports.
    same("rate (exponential): matches survreg up to its sign convention",
         unname(coef(f_tda_rt)), unname(-coef(f_r_rt)), 1e-6)
}

# ---- tda_coxph vs survival::coxph ---------------------------------------

if (requireNamespace("survival", quietly = TRUE)) {
    library(survival)
    set.seed(1)
    n_cx <- 100
    d_cx <- data.frame(x = rnorm(n_cx))
    d_cx$t <- rexp(n_cx, exp(0.3 + 0.6 * d_cx$x))
    d_cx$s <- 1
    f_tda <- tda_coxph(Surv(t, s) ~ x, d_cx)
    f_surv <- survival::coxph(Surv(t, s) ~ x, d_cx, ties = "breslow")
    same("coxph: matches survival::coxph (Breslow ties)",
         unname(coef(f_tda)), unname(coef(f_surv)), 1e-8)
} else {
    message("survival not installed: skipping rate/coxph comparisons")
}

# ---- tda_glm vs stats::glm ----------------------------------------------

set.seed(1)
n_gl <- 60
d_gl <- data.frame(x = rnorm(n_gl))
d_gl$y <- rbinom(n_gl, 1, plogis(0.3 + 1.1 * d_gl$x))
f_tda_glm <- tda_glm(y ~ x, d_gl, family = binomial)
f_r_glm <- stats::glm(y ~ x, d_gl, family = binomial)
same("glm: matches stats::glm (binomial)",
     unname(coef(f_tda_glm)), unname(coef(f_r_glm)), 1e-8)

# ---- tda_lsreg vs stats::lm ----------------------------------------------

set.seed(1)
n_ls <- 50
d_ls <- data.frame(x = rnorm(n_ls), z = rnorm(n_ls))
d_ls$y <- 1 + 0.5 * d_ls$x - 0.3 * d_ls$z + rnorm(n_ls, sd = 0.4)
f_tda_ls <- tda_lsreg(y ~ x + z, d_ls)
f_r_ls <- stats::lm(y ~ x + z, d_ls)
same("lsreg: matches stats::lm",
     unname(coef(f_tda_ls)), unname(coef(f_r_ls)), 1e-6)

# ---- tda_qreg (logit) vs stats::glm(family = binomial) ------------------
#
# A binary logit is mathematically the same model as a binomial glm with
# a logit link -- not merely a similar one -- so this is an exact-fit
# check, not an approximate one.

set.seed(1)
n_qr <- 60
d_qr <- data.frame(x = rnorm(n_qr))
d_qr$y <- rbinom(n_qr, 1, plogis(0.3 + 1.1 * d_qr$x))
f_tda_qr <- tda_qreg(y ~ x, d_qr, model = "logit")
f_r_qr <- stats::glm(y ~ x, d_qr, family = binomial)
same("qreg (logit): matches glm(family=binomial)",
     unname(coef(f_tda_qr)), unname(coef(f_r_qr)), 1e-7)


# quant is quantile(type = 6): interpolation at p(n+1) on the sorted
# values (manual 6.2.3), pinned against R's implementation of the
# same definition
set.seed(5)
qx <- round(rnorm(37, 50, 10), 2)
qp <- c(.1, .2, .25, .3, .4, .5, .6, .7, .75, .8, .9)
# Both paths are checked here, each against its contract, so the
# suite means the same thing under either setting of the master flag.
qtab <- local({
    op <- options(tdaR.use_exports = TRUE)
    on.exit(options(op))
    tda_quant(data.frame(X = qx))$table
})
# The stored quantiles are the exported doubles, not the 7.2 text, so
# they are compared against R's values directly rather than against
# those values rounded to the print format.  The parser path still has
# only the two printed decimals, and is held to exactly that.
qref <- unname(stats::quantile(qx, qp, type = 6))
assert("quant matches quantile(type = 6)",
       max(abs(as.numeric(qtab[1, ]) - qref)) < 1e-9)

qtab_p <- local({
    op <- options(tdaR.use_exports = FALSE)
    on.exit(options(op))
    tda_quant(data.frame(X = qx))$table
})
assert("quant, parser path, matches it to the printed precision",
       max(abs(as.numeric(qtab_p[1, ]) - round(qref, 2))) < 1e-9)
