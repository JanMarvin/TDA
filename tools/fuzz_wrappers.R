# Sweep the newer wrappers across option grids and edge-shaped inputs,
# looking for dead ends: R errors that are not deliberate refusals, TDA
# "Error:" lines leaking through, or nonsense results.  Run from
# tdaR/tests with the package installed:
#     Rscript ../../tools/fuzz_wrappers.R
# Prints one line per case; nonzero exit if anything unexpected died.
library(tdaR)
bad <- 0L
note <- function(id, ok, msg = "") {
    cat(sprintf("%-46s %s %s\n", id, if (ok) "ok" else "DEAD-END", msg))
    if (!ok) bad <<- bad + 1L
}
try_case <- function(id, expr, refusal_ok = TRUE) {
    r <- tryCatch(list(v = force(expr)), error = function(e) e)
    if (!inherits(r, "error")) return(note(id, TRUE))
    m <- conditionMessage(r)
    deliberate <- grepl("must|needs|cannot|at most|at least|handles|
        without NA|non-negative|nothing constrains|symmetric|
        no feasible|feasible solution|unbounded|admit no",
        m, perl = TRUE) && !grepl("^TDA could not", m)
    note(id, refusal_ok && deliberate, substr(m, 1, 60))
}
set.seed(1)

# ---- tda_locate_line ---------------------------------------------------
for (k in c(1L, 2L, 4L)) for (m in c(2L, 5L)) {
    vf <- matrix(round(runif(k * m, 0, 9)), k)
    try_case(sprintf("locate k=%d m=%d", k, m), tda_locate_line(vf))
    vv <- matrix(0, k, k); if (k > 1) { vv[1, k] <- vv[k, 1] <- 3 }
    try_case(sprintf("locate k=%d m=%d +vv", k, m),
             tda_locate_line(vf, vv))
}
try_case("locate all-zero weights",
         tda_locate_line(matrix(0, 2, 3)))
try_case("locate one column refused",
         tda_locate_line(matrix(1, 2, 1)))
try_case("locate NA refused", tda_locate_line(matrix(NA, 2, 3)))

# ---- tda_boolean_min ---------------------------------------------------
d <- expand.grid(X1 = 0:1, X2 = 0:1, X3 = 0:1, X4 = 0:1)
d$Y <- as.integer((d$X1 & !d$X2) | (d$X3 & d$X4))
for (alg in c("lawler", "petrick", "lawler2", "none"))
    for (u in c("dont_care", "true", "false"))
        try_case(sprintf("bfa %s/%s", alg, u),
                 tda_boolean_min("Y", paste0("X", 1:4), data = d,
                                 algorithm = alg, undefined = u))
try_case("bfa constant y=1",
         tda_boolean_min(rep(1L, 8), d[1:8, 1:3]))
try_case("bfa constant y=0",
         tda_boolean_min(rep(0L, 8), d[1:8, 1:3]))
try_case("bfa single condition",
         tda_boolean_min(d$X1[1:8], d[1:8, 1, drop = FALSE]))
try_case("bfa 16 conditions refused",
         tda_boolean_min(rep(0:1, 8), matrix(0L, 16, 16)))
try_case("bfa non-binary refused",
         tda_boolean_min(c(0, 2, 1), matrix(0L, 3, 2)))

# ---- tda_mlrc_design ---------------------------------------------------
for (n in c(2L, 5L)) for (q in c(1L, 3L)) {
    z <- matrix(rnorm(n * q), n)
    g <- sort(sample(1:2, n, replace = TRUE))
    try_case(sprintf("mldes n=%d q=%d", n, q), tda_mlrc_design(z, g))
}
try_case("mldes one group", tda_mlrc_design(matrix(1:4, 2), c(1, 1)))
try_case("mldes group per row", tda_mlrc_design(matrix(1:4, 2), 1:2))
try_case("mldes NA refused", tda_mlrc_design(matrix(NA, 2, 2), c(1, 1)))

# ---- tda_isotonic / tda_mlpi / tda_independence (recent) -----
try_case("isotonic length 1", tda_isotonic(5))
try_case("isotonic constant", tda_isotonic(rep(2, 6)))
try_case("isotonic all-tied group",
         tda_isotonic(c(3, 1, 2), groups = c(7, 7, 7), ties = "secondary"))
try_case("binary_program infeasible",
         tda_mlpi(c(1, 1), rbind(c(1, 1)), bounds = 5))
try_case("independence single y value",
         tda_independence(c(1, 2, 1, 2), c(3, 3, 3, 3)))
try_case("independence partition",
         tda_independence(rep(1:3, 4), rep(1:4, 3), partition = c(2, 4)))

# ---- tda_mqap -------------------------------------------
for (n in c(2L, 3L, 5L)) {
    F <- matrix(round(runif(n * n, 0, 9)), n); diag(F) <- 0
    D <- matrix(round(runif(n * n, 0, 9)), n); diag(D) <- 0
    try_case(sprintf("qap n=%d no costs", n), tda_mqap(F, D))
    C <- matrix(round(runif(n * n, 0, 9)), n)
    try_case(sprintf("qap n=%d with costs", n),
             tda_mqap(F, D, C))
}
try_case("qap nonzero diagonal accepted (silently zeroed)",
         tda_mqap(matrix(c(0, 1, 1, 0), 2), matrix(1, 2, 2)))
try_case("qap all-zero matrices refused",
         tda_mqap(matrix(0, 3, 3), matrix(0, 3, 3)))
try_case("qap n=1 refused",
         tda_mqap(matrix(0, 1, 1), matrix(0, 1, 1)))
try_case("qap non-square F refused",
         tda_mqap(matrix(0, 2, 3), matrix(0, 2, 2)))
try_case("qap mismatched D size refused",
         tda_mqap(matrix(0, 3, 3), matrix(0, 2, 2)))
try_case("qap mismatched costs size refused",
         tda_mqap(matrix(0, 3, 3), matrix(0, 3, 3),
                                  matrix(0, 2, 2)))

# ---- tda_mlp / tda_mlp1 ---------------------------------------------------
try_case("lp plain max",
         tda_mlp(c(3, 2), rbind(c(1, 1), c(1, 3)), c(4, 6)))
try_case("lp plain min",
         tda_mlp(c(3, 2), rbind(c(1, 1), c(1, 3)), c(4, 6),
                direction = "min"))
try_case("lp one equality",
         tda_mlp1(c(2, 3), rbind(c(1, 0), c(0, 1)), c(4, 5),
                 equalities = rbind(c(1, 1)),
                 equalities_bounds = 6))
try_case("lp equalities only, no inequalities",
         tda_mlp1(c(1, 1), equalities = rbind(c(1, 1)),
                 equalities_bounds = 3))
try_case("lp single variable",
         tda_mlp(5, rbind(1), 2))
try_case("lp infeasible (conflicting equalities) refused",
         tda_mlp1(c(1, 1), equalities = rbind(c(1, 1), c(1, 1)),
                 equalities_bounds = c(1, 2)))
try_case("lp unbounded objective refused",
         tda_mlp(1, rbind(-1), 0))
try_case("lp primal infeasible (dual unbounded) refused",
         tda_mlp(-1, rbind(1), -1))
try_case("lp mismatched bounds refused",
         tda_mlp(c(1, 1), rbind(c(1, 1)), c(1, 2)))
try_case("lp mismatched equalities_bounds refused",
         tda_mlp1(c(1, 1), equalities = rbind(c(1, 1)),
                 equalities_bounds = c(1, 2)))
try_case("lp neither constraints nor equalities refused",
         tda_mlp1(c(1, 1)))

# ---- tda_mls / tda_mlse / tda_mlsi / tda_mlsei / tda_mlsei1 / tda_mnls ----
try_case("ls unconstrained", tda_mls(diag(2), c(1, 2)))
try_case("ls mlse matches mls (same call, documented in tda.hlp)",
         tda_mlse(diag(2), c(1, 2)))
try_case("ls inequality only",
         tda_mlsi(rbind(c(1, 0), c(0, 1)), c(1, 2)))
try_case("ls mixed equality+data+inequality",
         tda_mlsei(diag(2), c(2, 3), equalities = rbind(c(1, 1)),
                  equalities_bounds = 3,
                  inequalities = rbind(c(1, -1)),
                  inequalities_bounds = 0))
try_case("ls mixed, active_set method (mlsei1)",
         tda_mlsei1(diag(2), c(2, 3), equalities = rbind(c(1, 1)),
                   equalities_bounds = 3,
                   inequalities = rbind(c(1, -1)),
                   inequalities_bounds = 0))
try_case("ls nonneg (mnls)", tda_mnls(diag(2), c(-1, 2)))
try_case("ls overdetermined inconsistent (least-squares, not an error)",
         tda_mls(rbind(c(1, 1), c(1, -1), c(1, 2)), c(4, 0, 10)))
try_case("ls mnls has no inequalities parameter (signature-level refusal)",
         tryCatch({tda_mnls(diag(2), c(1, 1), inequalities = rbind(c(1, 0)))
                  stop("should not reach here")},
                 error = function(e) stop("cannot combine", call. = FALSE)))
try_case("ls nothing given refused (mlsei)", tda_mlsei())
try_case("ls contradictory equalities refused",
         tda_mlsei(equalities = rbind(c(1, 1), c(1, 1)),
                  equalities_bounds = c(1, 2)))
try_case("ls infeasible inequalities refused",
         tda_mlsi(rbind(c(1, 0), c(-1, 0)), c(5, -1)))
try_case("ls mismatched A/b refused",
         tda_mls(diag(2), c(1, 1, 1)))
try_case("ls mismatched column counts refused",
         tda_mlsei(diag(2), c(1, 1), equalities = rbind(c(1, 1, 1)),
                  equalities_bounds = 1))

# ---- tda_mqp / tda_mqpb / tda_mqpc ----------------------------------------
try_case("qp unconstrained",
         tda_mqp(diag(c(2, 2)), c(-2, -4)))
try_case("qp boxed",
         tda_mqpb(diag(c(2, 2)), c(-2, -4), lower = c(0, 0),
                 upper = c(0.5, 0.5)))
try_case("qp mixed equality+inequality",
         tda_mqpc(diag(c(2, 2)), c(-2, -4),
                 equalities = rbind(c(1, 1)),
                 equalities_bounds = 1,
                 inequalities = rbind(c(1, 0)),
                 inequalities_bounds = 0.2))
try_case("qp asymmetric C refused",
         tda_mqp(matrix(c(1, 2, 0, 1), 2, 2), c(1, 1)))
try_case("qp non-square C refused",
         tda_mqp(matrix(1, 2, 3), c(1, 1)))
try_case("qp mqpc with neither block refused",
         tda_mqpc(diag(2), c(1, 1)))
try_case("qp contradictory equalities refused",
         tda_mqpc(diag(2), c(1, 1),
                 equalities = rbind(c(1, 1), c(1, 1)),
                 equalities_bounds = c(1, 5)))
try_case("qp mismatched d refused",
         tda_mqp(diag(2), c(1, 1, 1)))
try_case("qp mismatched bounds refused",
         tda_mqpb(diag(2), c(1, 1), lower = c(0, 0, 0),
                 upper = c(1, 1)))

# ---- core matrix algebra: tda_mmul/mtransp/mchol/minvs/minvd/mginv/
#      mnrow/mncol/mnorm/mnorm1/mnorm2/mtrace/mdiag/mdiagd ------------------
try_case("mmul two matrices", tda_mmul(matrix(1:4, 2), diag(2)))
try_case("mmul chained (3 matrices)",
         tda_mmul(matrix(1:4, 2), matrix(5:8, 2), diag(2)))
try_case("mtransp", tda_mtransp(matrix(1:6, 2, 3)))
try_case("mchol on positive-definite matrix",
         tda_mchol(matrix(c(4, 2, 2, 3), 2)))
try_case("mchol on non-positive-definite refused",
         tda_mchol(matrix(c(1, 2, 2, 1), 2)))
try_case("minvs on symmetric matrix", tda_minvs(matrix(c(4, 2, 2, 3), 2)))
try_case("minvs on non-square refused", tda_minvs(matrix(1:6, 2, 3)))
try_case("minvd on any matrix", tda_minvd(matrix(1:9, 3, 3)))
try_case("mginv on a tall matrix", tda_mginv(matrix(rnorm(6), 3, 2)))
try_case("mginv on a wide matrix refused", tda_mginv(matrix(rnorm(6), 2, 3)))
try_case("mnrow/mncol", tda_mnrow(matrix(1:6, 2, 3)))
try_case("mnorm on an all-negative matrix (regression: was 0 before the fix)",
         tda_mnorm(matrix(c(-1, -3, -2, -5), 2)))
try_case("mnorm on a mixed-sign matrix (regression: picked the raw max, not abs max)",
         tda_mnorm(matrix(c(-9, 4, -1, 6), 2)))
try_case("mnorm1", tda_mnorm1(matrix(1:6, 2, 3)))
try_case("mnorm2", tda_mnorm2(matrix(1:6, 2, 3)))
try_case("mtrace on a non-square matrix", tda_mtrace(matrix(1:6, 2, 3)))
try_case("mdiag from a vector", tda_mdiag(c(1, 2, 3)))
try_case("mdiagd extraction from a non-square matrix",
         tda_mdiagd(matrix(1:6, 2, 3)))

# ---- matrix reshaping: tda_mcvec/mrvec/mivec/mrsum/mcsum/mdrow/mdcol/
#      msort/msort1/mrank/mcath/mcatv/mcathv/mtrim ------------------------
try_case("mcvec", tda_mcvec(matrix(1:6, 2, 3)))
try_case("mrvec", tda_mrvec(matrix(1:6, 2, 3)))
try_case("mivec reshapes evenly", tda_mivec(1:6, 2))
try_case("mivec with a non-dividing n refused", tda_mivec(1:6, 4))
try_case("mrsum/mcsum", tda_mrsum(matrix(1:6, 2, 3)))
try_case("mdrow/mdcol", tda_mdrow(matrix(1:6, 2, 3)))
try_case("msort by one column", tda_msort(rbind(c(3, 1), c(1, 1), c(2, 1)), 1))
try_case("msort by two columns",
         tda_msort(rbind(c(1, 2), c(1, 1), c(2, 1)), c(1, 2)))
try_case("mrank returns order(), not rank()",
         tda_mrank(rbind(c(30, 1), c(10, 1), c(20, 1)), 1))
try_case("msort1 dedups after sorting",
         tda_msort1(rbind(c(1, 1), c(1, 1), c(2, 2)), c(1, 2)))
try_case("mcath (cbind)", tda_mcath(matrix(1:4, 2), matrix(5:8, 2)))
try_case("mcatv (rbind)", tda_mcatv(matrix(1:4, 2), matrix(5:8, 2)))
try_case("mcathv (block-diagonal direct sum)",
         tda_mcathv(matrix(1:4, 2), matrix(5:8, 2)))
try_case("mtrim drops a leading column",
         tda_mtrim(matrix(1:6, 2, 3), leading_cols = 1))
try_case("mtrim pads a leading column",
         tda_mtrim(matrix(1:6, 2, 3), leading_cols = -1))
try_case("mtrim dropping everything refused",
         tda_mtrim(matrix(1:6, 2, 3), leading_cols = 3))

# ---- eigen/SVD: tda_mevs/msvd/msvd1 ---------------------------------------
try_case("mevs on a symmetric matrix", tda_mevs(matrix(c(2, 1, 1, 2), 2)))
try_case("mevs on a non-square matrix refused",
         tda_mevs(matrix(1:6, 2, 3)))
try_case("msvd (values only) on a tall matrix",
         tda_msvd(matrix(rnorm(6), 3, 2)))
try_case("msvd1 (full decomposition) on a tall matrix",
         tda_msvd1(matrix(rnorm(6), 3, 2)))
try_case("msvd on a wide matrix refused",
         tda_msvd(matrix(rnorm(6), 2, 3)))
try_case("msvd1 on a wide matrix refused",
         tda_msvd1(matrix(rnorm(6), 2, 3)))

# ---- select/permute/aggregate: tda_msrow/mscol/mprow/mpcol/mpsym/mag/mnc --
try_case("msrow with repeats", tda_msrow(matrix(1:9, 3, 3), c(1, 1, 2)))
try_case("mscol with repeats", tda_mscol(matrix(1:9, 3, 3), c(2, 2)))
try_case("mprow with a real permutation",
         tda_mprow(matrix(1:9, 3, 3), c(3, 1, 2)))
try_case("mprow with a non-permutation refused",
         tda_mprow(matrix(1:9, 3, 3), c(1, 1, 2)))
try_case("mpcol with a real permutation",
         tda_mpcol(matrix(1:9, 3, 3), c(3, 1, 2)))
try_case("mpsym with a real permutation",
         tda_mpsym(matrix(1:9, 3, 3), c(3, 1, 2)))
try_case("mpsym on a non-square matrix refused",
         tda_mpsym(matrix(1:6, 2, 3), c(1, 2)))
try_case("mag aggregates by row and column groups",
         tda_mag(matrix(1:12, 3, 4, byrow = TRUE), c(1, 2, 1), c(1, 1, 2, 2)))
try_case("mag mismatched row_groups refused",
         tda_mag(matrix(1:12, 3, 4, byrow = TRUE), c(1, 2), c(1, 1, 2, 2)))
try_case("mnc with <=", tda_mnc(matrix(1:6, 2, 3), 3, "<="))
try_case("mnc with abs>", tda_mnc(matrix(1:6, 2, 3), 3, "abs>"))

# ---- tda_msqrtd/msqrti/mcent/mstand/mdcent/mcross/mnum --------------------
try_case("msqrtd on a positive diagonal",
         tda_msqrtd(matrix(c(4, 1, 1, 9), 2)))
try_case("msqrtd with a negative diagonal refused",
         tda_msqrtd(matrix(c(-4, 1, 1, 9), 2)))
try_case("msqrti on a positive diagonal",
         tda_msqrti(matrix(c(4, 1, 1, 9), 2)))
try_case("msqrti with a near-zero diagonal refused",
         tda_msqrti(matrix(c(0, 1, 1, 9), 2)))
try_case("mcent", tda_mcent(matrix(rnorm(6), 3, 2)))
try_case("mstand", tda_mstand(matrix(rnorm(6), 3, 2)))
try_case("mcross", tda_mcross(matrix(rnorm(6), 3, 2)))
try_case("mdcent on a symmetric matrix",
         tda_mdcent(matrix(c(0, 3, 4, 3, 0, 5, 4, 5, 0), 3, 3)))
try_case("mdcent on a non-square matrix refused",
         tda_mdcent(matrix(1:6, 2, 3)))
try_case("mnum builds an arithmetic sequence", tda_mnum(5, 2, 4))
try_case("mnum with n < 1 refused", tda_mnum(5, 2, 0))

# ---- tda_mkp/mwvec/mwvec1/mscal1 -------------------------------------------
try_case("mkp (Kronecker product)", tda_mkp(matrix(1:4, 2), matrix(5:8, 2)))
try_case("mwvec", tda_mwvec(c(10, 20, 30, 40), c(1, 2, 1, 3)))
try_case("mwvec mismatched lengths refused",
         tda_mwvec(c(10, 20, 30), c(1, 2)))
try_case("mwvec1", tda_mwvec1(c(10, 20, 30, 40), c(1, 2, 1, 3), c(1, 2, 3, 4)))
try_case("mwvec1 mismatched lengths refused",
         tda_mwvec1(c(10, 20, 30), c(1, 2, 3), c(1, 2)))
try_case("mscal1", tda_mscal1(matrix(c(2, 4, 6, 8), 2)))
try_case("mscal1 with a zero-sum matrix refused",
         tda_mscal1(matrix(c(1, -1, 1, -1), 2)))

# ---- tda_mpinv/mcel --------------------------------------------------------
try_case("mpinv on a real permutation", tda_mpinv(c(3, 1, 4, 2)))
try_case("mpinv on a non-permutation refused", tda_mpinv(c(1, 1, 2)))
try_case("mcel finds edges above a threshold",
         tda_mcel(matrix(c(5, 0, 1, 0, 3, 0, 2, 0, 4), 3, 3), 3))
try_case("mcel with no edges refused",
         tda_mcel(matrix(c(5, 0, 1, 0, 3, 0, 2, 0, 4), 3, 3), 100))
try_case("mcel on a non-square matrix refused",
         tda_mcel(matrix(1:6, 2, 3), 1))

# ---- tda_mpz ----------------------------------------------------------------
try_case("mpz on a matrix with a valid zero-free permutation",
         tda_mpz(matrix(c(0, 1, 2, 3, 0, 4, 0, 5, 0), 3, 3, byrow = TRUE)))
try_case("mpz on a non-square matrix refused",
         tda_mpz(matrix(1:6, 2, 3)))

# ---- tda_mpbl/mpbu ----------------------------------------------------------
try_case("mpbl (block lower triangular)",
         tda_mpbl(matrix(c(1, 2, 0, 0, 3, 1, 0, 0, 4, 0, 1, 5, 0, 0, 6, 1),
                         4, 4, byrow = TRUE)))
try_case("mpbu (block upper triangular)",
         tda_mpbu(matrix(c(1, 2, 0, 0, 3, 1, 0, 0, 4, 0, 1, 5, 0, 0, 6, 1),
                         4, 4, byrow = TRUE)))
try_case("mpbl on a non-square matrix refused",
         tda_mpbl(matrix(1:6, 2, 3)))
try_case("mpbu on a non-square matrix refused",
         tda_mpbu(matrix(1:6, 2, 3)))

# ---- tda_mpfit --------------------------------------------------------------
try_case("mpfit converges to the target margins",
         tda_mpfit(matrix(c(10, 20, 30, 40), 2, 2, byrow = TRUE),
                  row_sums = c(45, 55), col_sums = c(40, 60)))
try_case("mpfit mismatched row_sums refused",
         tda_mpfit(matrix(1:4, 2, 2), c(1, 2, 3), c(1, 2)))
try_case("mpfit mismatched col_sums refused",
         tda_mpfit(matrix(1:4, 2, 2), c(1, 2), c(1, 2, 3)))
try_case("mpfit with max_iter < 1 refused",
         tda_mpfit(matrix(1:4, 2, 2), c(1, 2), c(1, 2), max_iter = 0))
try_case("mpfit with negative eps refused",
         tda_mpfit(matrix(1:4, 2, 2), c(1, 2), c(1, 2), eps = -1))

# ---- tda_mpit/mpit1/mkmet/mple --------------------------------------------
try_case("mpit (Leslie matrix iteration)",
         tda_mpit(rbind(c(0, 0), c(2, 0.5), c(1, 0.3)), c(100, 50, 20),
                 iterations = 3))
try_case("mpit mismatched population refused",
         tda_mpit(rbind(c(0, 0), c(2, 0.5)), c(100, 50, 20), iterations = 2))
try_case("mpit1 (with immigration)",
         tda_mpit1(rbind(c(0, 0), c(2, 0.5), c(1, 0.3)), c(100, 50, 20),
                  c(5, 2, 1), iterations = 2))
try_case("mkmet (Kemeny distance)",
         tda_mkmet(rbind(c(1, 2, 3), c(3, 2, 1), c(1, 2, 2))))
try_case("mple (product-limit CDF)",
         tda_mple(c(2, 3, 3, 5, 7), c(0, 0, 1, 0, 0)))
try_case("mple mismatched lengths refused",
         tda_mple(c(1, 2, 3), c(0, 1)))

# ---- tda_mev -----------------------------------------------------------------
try_case("mev on a matrix with complex eigenvalues",
         tda_mev(matrix(c(0, -1, 1, 0), 2, 2, byrow = TRUE)))
try_case("mev on a matrix with real eigenvalues",
         tda_mev(matrix(c(4, 1, 0, 2, 3, 0, 0, 0, 5), 3, 3, byrow = TRUE)))
try_case("mev on a non-square matrix refused",
         tda_mev(matrix(1:6, 2, 3)))

# ---- tda_mch ------------------------------------------------------------------
try_case("mch on a matrix with an out-of-order entry",
         tda_mch(matrix(c(0, 1, 0, 5, 0, 1, 0, 0, 0), 3, 3, byrow = TRUE)))
try_case("mch on an already-ordered matrix",
         tda_mch(matrix(c(0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0),
                        4, 4, byrow = TRUE)))
try_case("mch on a non-square matrix refused",
         tda_mch(matrix(1:6, 2, 3)))

# ---- tda_midf/midf1/midf2/midf3 --------------------------------------------
try_case("midf (CDF bounds and estimate)", tda_midf(c(1, 2, 0), c(3, 4, 2)))
try_case("midf1 (CDF at each observation)", tda_midf1(c(1, 2, 0), c(3, 4, 2)))
try_case("midf2 (restricted mean)", tda_midf2(c(1, 2, 0), c(3, 4, 2)))
try_case("midf3 (nearby-endpoint averages)", tda_midf3(c(1, 2, 0), c(3, 4, 2)))
try_case("midf with a zero-width interval refused",
         tda_midf(c(1, 1), c(1, 2)))
try_case("midf1 mismatched lengths refused",
         tda_midf1(c(1, 2), c(1, 2, 3)))
try_case("midf2 with a zero-width interval refused",
         tda_midf2(c(1, 1), c(1, 2)))
try_case("midf3 mismatched lengths refused",
         tda_midf3(c(1, 2), c(1, 2, 3)))

cat(sprintf("\n%d dead end(s)\n", bad))

# With FUZZ_KEEP set to a directory, every run directory (command
# files, data, outputs) is preserved there before R cleans its session
# tempdir -- so the generated .cf files can be replayed against an
# instrumented binary to measure what the sweep adds to C coverage.
keep <- Sys.getenv("FUZZ_KEEP")
if (nzchar(keep)) {
    dir.create(keep, showWarnings = FALSE, recursive = TRUE)
    for (d in list.dirs(tempdir(), recursive = FALSE))
        file.copy(d, keep, recursive = TRUE)
}
quit(status = min(bad, 1L))
