# Optimization commands.

#' Quadratic assignment problem
#'
#' Approximate solution of the quadratic assignment problem -- TDA's
#' \code{mqap}, adapted from CACM algorithm 608 (D.H. West). Tries to
#' find a permutation \code{p} of \code{1:n} minimizing
#' \deqn{\sum_i \left( C_{i,p(i)} + \sum_j F_{ij} D_{p(i),p(j)} \right)}{
#' sum_i ( C[i, p(i)] + sum_j F[i,j] * D[p(i), p(j)] )}
#' This is a heuristic, not an exact solver: the value returned is an
#' upper bound on the true minimum, and need not be optimal for larger
#' \code{n} (see \code{examples/coverage/mqapops.cf} for the objective
#' and small hand-verified instances). \code{flows} and \code{distances}
#' must have a zero main diagonal; a nonzero diagonal is silently reset
#' to zero with a warning from TDA itself.
#'
#' @param flows the flow matrix \code{F} (n x n).
#' @param distances the distance matrix \code{D} (n x n).
#' @param costs the fixed placement cost matrix \code{C} (n x n),
#'   default all zero.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A list: \code{value}, the objective at the returned
#'   permutation; \code{permutation}, integer vector \code{p} with
#'   \code{p[i]} the location assigned to facility \code{i}.
#' @family optimization
#' @examples
#' F <- matrix(c(0,5,2,4, 5,0,3,0, 2,3,0,0, 4,0,0,0), 4, 4, byrow = TRUE)
#' D <- matrix(c(0,8,15,13, 8,0,9,10, 15,9,0,17, 13,10,17,0), 4, 4,
#'             byrow = TRUE)
#' tda_mqap(F, D)
#' @export
tda_mqap <- function(flows, distances, costs = NULL,
                     dir = tempfile("tda"), ...) {
    flows <- .tda_check_matrix(flows, "flows", square = TRUE)
    distances <- .tda_check_matrix(distances, "distances", square = TRUE,
                                   same_dim_as = flows,
                                   other_arg = "flows")
    flows <- as.matrix(flows)
    distances <- as.matrix(distances)
    n <- nrow(flows)
    if (n < 2L || ncol(flows) != n)
        stop("`flows` must be a square matrix with at least 2 rows",
             call. = FALSE)
    if (!identical(dim(distances), dim(flows)))
        stop("`distances` must be the same size as `flows`", call. = FALSE)
    if (is.null(costs))
        costs <- matrix(0, n, n)
    else {
        costs <- as.matrix(costs)
        if (!identical(dim(costs), dim(flows)))
            stop("`costs` must be the same size as `flows`", call. = FALSE)
    }
    if (all(costs == 0) && (all(flows == 0) || all(distances == 0)))
        stop("`flows`, `distances` and `costs` cannot all be degenerate ",
             "like this: mqap has nothing to optimize", call. = FALSE)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    utils::write.table(flows, file.path(dir, "f.mat"),
                       row.names = FALSE, col.names = FALSE)
    utils::write.table(distances, file.path(dir, "d.mat"),
                       row.names = FALSE, col.names = FALSE)
    utils::write.table(costs, file.path(dir, "c.mat"),
                       row.names = FALSE, col.names = FALSE)
    res <- tda_run(c("mfmt = 24.16;",
                     sprintf("mdef(F,%d,%d) = f.mat;", n, n),
                     sprintf("mdef(D,%d,%d) = d.mat;", n, n),
                     sprintf("mdef(C,%d,%d) = c.mat;", n, n),
                     "mqap(F,D,C,P);", "mpr(P) = p.out;"),
                   dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run mqap: ", err[1L], call. = FALSE)
    val <- NULL
    perm <- NULL
    if (.use_exports()) {
        vv <- res$exports[["mqap.value"]]
        pm <- res$exports[["mpr.matrix"]]
        if (is.matrix(vv) && length(vv) == 1L)
            val <- as.numeric(vv[1L])
        if (is.matrix(pm))
            perm <- as.integer(pm)
    }
    if (is.null(val)) {
        vline <- grep("^Best value:", res$output)
        if (length(vline) != 1L)
            stop("mqap produced no readable result", call. = FALSE)
        val <- as.numeric(sub("^Best value: *", "", res$output[vline]))
    }
    if (is.null(perm)) {
        m <- tda_file(res, "p.out")
        if (is.null(m))
            stop("mqap produced no permutation", call. = FALSE)
        perm <- as.integer(as.matrix(m)[, 1L])
    }
    list(value = val, permutation = perm)
}

#' Linear programming
#'
#' Maximize (or minimize) a linear objective subject to linear
#' inequality and/or equality constraints -- TDA's \code{mlp} and
#' \code{mlp1} (Salazar & Sen's MINIT, CACM algorithm 333).
#' \code{tda_mlp} takes only inequality constraints, matching
#' \code{mlp(T,X,Y)}; \code{tda_mlp1} adds equality constraints,
#' matching \code{mlp1(T,p,X,Y)} -- the two are the same underlying
#' call with \code{p} fixed at 0 for \code{tda_mlp}. The problem solved,
#' after negating the objective internally for \code{direction =
#' "min"}, is
#' \deqn{\max\; c'x \quad \mathrm{s.t.}\; x \ge 0,\; Ax \le b,\; A_e x = b_e}{
#' max c'x  s.t.  x >= 0, A x <= b, Ae x = be}
#' with the dual \eqn{\min\; b'y \;\mathrm{s.t.}\; y \ge 0,\; A'y \ge c}{
#' min b'y s.t. y >= 0, A'y >= c} (this is all stated directly in
#' \code{lpf1}'s header in \code{t_lp.c}; nothing here was
#' guessed). Equality-constraint rows report a dual value of 0 --
#' TDA's convention, not something this wrapper computes. Verified
#' by hand for a plain and an equality instance, and against
#' brute-force vertex enumeration in \code{test-r-comparisons.R}.
#'
#' @param objective coefficients \code{c}, one per variable (at least
#'   one).
#' @param constraints the inequality constraint matrix \code{A}
#'   (\code{Ax <= bounds}), one row per constraint, one column per
#'   variable. For \code{tda_mlp1}, may be omitted if \code{equalities}
#'   is given instead.
#' @param bounds the right-hand sides \code{b} for \code{constraints}.
#' @param equalities an equality constraint matrix \code{Ae}
#'   (\code{Ae x = equalities_bounds}), same number of columns as
#'   \code{constraints}. \code{tda_mlp1} only.
#' @param equalities_bounds the right-hand sides \code{be} for
#'   \code{equalities}. \code{tda_mlp1} only.
#' @param direction \code{"max"} (default) or \code{"min"}.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A list: \code{value}, the objective at the optimum;
#'   \code{x}, the primal solution; \code{y}, the dual solution (one
#'   entry per constraint row, inequalities then equalities, equality
#'   rows always 0).
#' @details
#' \code{tda_mlp} also accepts TDA's single-tableau form,
#' \code{mlp(T, X, Y)}: the first row is the objective
#' \code{(c, 0)} and the remaining rows are the constraints
#' \code{[A, b]}, exactly the layout the solver documents.
#' \code{tda_mlp(T)} decomposes it accordingly.
#'
#' The manual's Box example writes \code{mlp(T, X, Y, 1)} with a
#' fourth operand.  That form does not exist in this source tree:
#' the shipped 6.4q parser (checked against a pristine copy of the
#' upstream source) accepts exactly \code{mlp(T, X, Y)} and
#' \code{mlp1(T, p, X, Y)}, so the box was produced by a different
#' release.  On the box's problem the three-operand form gives
#' the same value (1) and the same dual vector (1, 0); only the
#' primal solution differs -- (1, 0) here against the manual's
#' (1, 1) -- and both are vertices of the same optimal face, since
#' maximising x1 under x1 <= 1, x2 <= 1 leaves x2 free at the
#' optimum.  Nothing disagrees beyond the tie-break.
#' @family optimization
#' @examples
#' # maximize 3x1 + 2x2 s.t. x1 + x2 <= 4, x1 + 3x2 <= 6
#' tda_mlp(c(3, 2), rbind(c(1, 1), c(1, 3)), c(4, 6))
#' # same, plus x1 + x2 = 3 exactly
#' tda_mlp1(c(3, 2), rbind(c(1, 3)), c(6), equalities = rbind(c(1, 1)),
#'          equalities_bounds = 3)
#' @export
tda_mlp <- function(objective, constraints, bounds,
                    direction = c("max", "min"),
                    dir = tempfile("tda"), ...) {
    if (missing(constraints) && is.matrix(objective) &&
        nrow(objective) >= 2L && ncol(objective) >= 2L) {
        # the manual's form: a single simplex tableau T laid out
        # as lpf1 documents it -- first row (c, 0) is the objective,
        # the remaining rows are the constraints [A, b]
        T <- objective; n <- ncol(T) - 1L
        return(tda_mlp1(T[1L, seq_len(n)],
                        T[-1L, seq_len(n), drop = FALSE],
                        T[-1L, n + 1L],
                        direction = direction, dir = dir, ...))
    }
    tda_mlp1(objective, constraints, bounds, direction = direction,
            dir = dir, ...)
}

#' @rdname tda_mlp
#' @export
tda_mlp1 <- function(objective, constraints = NULL, bounds = NULL,
                     equalities = NULL, equalities_bounds = NULL,
                     direction = c("max", "min"),
                     dir = tempfile("tda"), ...) {
    n <- length(objective)
    if (n < 1L)
        stop("`objective` needs at least one variable", call. = FALSE)
    if (!is.null(constraints)) {
        constraints <- rbind(constraints)
        if (ncol(constraints) != n)
            stop("`constraints` must have one column per objective ",
                 "coefficient", call. = FALSE)
        if (nrow(constraints) != length(bounds))
            stop("`bounds` must have one entry per row of `constraints`",
                 call. = FALSE)
    }
    if (!is.null(equalities)) {
        equalities <- rbind(equalities)
        if (ncol(equalities) != n)
            stop("`equalities` must have one column per objective ",
                 "coefficient", call. = FALSE)
        if (nrow(equalities) != length(equalities_bounds))
            stop("`equalities_bounds` must have one entry per row of ",
                 "`equalities`", call. = FALSE)
    }
    A <- rbind(constraints, equalities)
    b <- c(bounds, equalities_bounds)
    m <- nrow(A)
    if (is.null(m) || m < 1L)
        stop("at least one of `constraints` or `equalities` is required",
             call. = FALSE)
    p <- if (is.null(equalities)) 0L else nrow(equalities)
    direction <- match.arg(direction)
    cobj <- if (direction == "min") -objective else objective
    tab <- rbind(c(cobj, 0), cbind(A, b))
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    utils::write.table(tab, file.path(dir, "t.mat"),
                       row.names = FALSE, col.names = FALSE)
    res <- tda_run(c("silent = -1;", "mfmt = 24.16;",
                     sprintf("mdef(T,%d,%d) = t.mat;", m + 1L, n + 1L),
                     sprintf("mlp1(T,%d,X,Y);", p),
                     "mpr(X) = x.out;", "mpr(Y) = y.out;"),
                   dir = dir, ...)
    if (any(grepl("Dual objective function is unbounded", res$output)))
        stop("no solution: the primal has no feasible point ",
             "(the dual is unbounded)", call. = FALSE)
    if (any(grepl("Primal objective function is unbounded", res$output)))
        stop("the objective is unbounded over the feasible region",
             call. = FALSE)
    if (any(grepl("No solution\\.", res$output))) {
        hint <- if (p >= 2L)
            " (possibly linear dependent equality constraints)" else ""
        stop("the constraints admit no feasible solution", hint,
             call. = FALSE)
    }
    err <- grep("^Error|insufficient memory", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run mlp: ", err[1L], call. = FALSE)
    val <- NULL
    xv <- NULL
    yv <- NULL
    if (.use_exports()) {
        vv <- res$exports[["mlp.value"]]
        xm <- res$exports[["mpr.matrix"]]
        ym <- res$exports[["mpr.matrix.2"]]
        if (is.matrix(vv) && length(vv) == 1L)
            val <- as.numeric(vv[1L])
        if (is.matrix(xm))
            xv <- as.numeric(xm)
        if (is.matrix(ym))
            yv <- as.numeric(ym)
    }
    if (is.null(val)) {
        vline <- grep("^Value is:", res$output)
        if (length(vline) != 1L)
            stop("mlp produced no readable result", call. = FALSE)
        val <- as.numeric(sub("^Value is: *", "", res$output[vline]))
    }
    if (is.null(xv)) {
        mx <- tda_file(res, "x.out")
        if (is.null(mx))
            stop("mlp produced no primal solution", call. = FALSE)
        xv <- as.numeric(as.matrix(mx)[, 1L])
    }
    if (is.null(yv)) {
        my <- tda_file(res, "y.out")
        if (is.null(my))
            stop("mlp produced no dual solution", call. = FALSE)
        yv <- as.numeric(as.matrix(my)[, 1L])
    }
    if (direction == "min")
        val <- -val
    list(value = val, x = xv, y = yv)
}

#' Constrained least squares
#'
#' Least squares with optional linear equality constraints, linear
#' inequality constraints, and non-negativity -- TDA's \code{mls},
#' \code{mlse}, \code{mlsi}, \code{mlsei}, \code{mlsei1}, and
#' \code{mnls}, each as its thin wrapper sharing one
#' implementation. The problem solved is
#' \deqn{\min\; \|Ax - b\|_2 \quad \mathrm{s.t.}\;
#' A_e x = b_e,\; A_i x \ge b_i,\; (x \ge 0)}{
#' min ||Ax - b||  s.t.  Ae x = be, Ai x >= bi, (x >= 0)}
#' These are TDA's documented semantics (\code{tda.hlp}, entries
#' \code{##mls}/\code{##mlse}/\code{##mlsi}/\code{##mlsei}/\code{##mnls}):
#' an earlier working note had assumed \code{mlse}/\code{mlsi} took
#' their constraint counts, which the source does not support --
#' \code{tda_mlse} and \code{tda_mls} turn out to make the identical
#' call (confirmed empirically: an overdetermined, inconsistent
#' instance returns the plain least-squares fit rather than an error),
#' and \code{tda_mlsi} treats its entire input as inequality rows.
#' \code{tda_mlsei} and \code{tda_mlsei1} solve the same general
#' equality+least-squares+inequality problem via two different TDA
#' algorithms (Lawson-Hanson vs. an undocumented active-set QP,
#' confirmed to agree to solver tolerance on two independent
#' instances, not merely assumed from the similar name). \code{tda_mnls}
#' is the only one of the six with non-negativity, applied to every
#' variable (TDA's \code{k} split point fixed at 0) and cannot be
#' combined with \code{inequalities} (\code{mnls} has no general
#' inequality block of its own).
#'
#' @param A the least-squares data matrix (\code{ma} x n). Required
#'   for \code{tda_mls}/\code{tda_mlse}/\code{tda_mnls}; optional for
#'   \code{tda_mlsei}/\code{tda_mlsei1} (\code{NULL} for a pure
#'   feasibility problem, together with \code{equalities} and/or
#'   \code{inequalities}); not used by \code{tda_mlsi}.
#' @param b the target vector for \code{A x ~ b}, length \code{ma}.
#' @param equalities an equality constraint matrix
#'   \code{Ae x = equalities_bounds}. \code{tda_mlsei}/
#'   \code{tda_mlsei1}/\code{tda_mnls} only.
#' @param equalities_bounds the right-hand sides for \code{equalities}.
#' @param inequalities an inequality constraint matrix
#'   \code{Ai x >= inequalities_bounds}. \code{tda_mlsi} (required,
#'   its only input) and \code{tda_mlsei}/\code{tda_mlsei1} (optional).
#' @param inequalities_bounds the right-hand sides for
#'   \code{inequalities}.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A list: \code{x}, the solution; \code{residual}, the
#'   Euclidean norm of \code{A \%*\% x - b} (\code{NA} when \code{A}
#'   is absent), computed directly in R from \code{A}, \code{x} and
#'   \code{b} rather than trusted from TDA's text, since
#'   \code{tda_mlsei1} prints no residual at all; \code{rank}, the
#'   reported rank of \code{A} (\code{NA} for \code{tda_mnls} and
#'   \code{tda_mlsei1}, neither of which report one).
#' @details
#' \code{tda_mls}, \code{tda_mlse} and \code{tda_mlsi} also accept the
#' manual's single-matrix form,
#' \code{mls(S, B)} with \code{S = [X, y]}: when the right-hand
#' side is not given, the last column of the first argument is
#' split off as \code{y}.
#' @family optimization
#' @examples
#' # unconstrained: fit x1, x2 to two noisy targets
#' tda_mls(diag(2), c(2.1, 2.9))
#' # the manual's augmented form S = [X, y] works too
#' tda_mls(cbind(diag(2), c(2.1, 2.9)))
#' # x1 + x2 = 3, x1 close to 2, x2 close to 3
#' tda_mlsei(diag(2), c(2, 3), equalities = rbind(c(1, 1)),
#'          equalities_bounds = 3)
#' @export
tda_mls <- function(A, b, dir = tempfile("tda"), ...) {
    if (missing(b)) {
        # the manual's form: a single augmented matrix S = [X, y]
        b <- A[, ncol(A)]
        A <- A[, -ncol(A), drop = FALSE]
    }
    .tda_ls_impl(A = A, b = b, dir = dir, ...)
}

#' @rdname tda_mls
#' @export
tda_mlse <- function(A, b, dir = tempfile("tda"), ...) {
    if (missing(b)) {
        b <- A[, ncol(A)]
        A <- A[, -ncol(A), drop = FALSE]
    }
    .tda_ls_impl(A = A, b = b, dir = dir, ...)
}

#' @rdname tda_mls
#' @export
tda_mlsi <- function(inequalities, inequalities_bounds,
                     dir = tempfile("tda"), ...) {
    if (missing(inequalities_bounds)) {
        # the manual's form: a single augmented matrix S = [X, y]
        inequalities_bounds <- inequalities[, ncol(inequalities)]
        inequalities <- inequalities[, -ncol(inequalities), drop = FALSE]
    }
    .tda_ls_impl(inequalities = inequalities,
                inequalities_bounds = inequalities_bounds, dir = dir, ...)
}

#' @rdname tda_mls
#' @export
tda_mlsei <- function(A = NULL, b = NULL, equalities = NULL,
                      equalities_bounds = NULL, inequalities = NULL,
                      inequalities_bounds = NULL, dir = tempfile("tda"), ...)
    .tda_ls_impl(A = A, b = b, equalities = equalities,
                equalities_bounds = equalities_bounds,
                inequalities = inequalities,
                inequalities_bounds = inequalities_bounds,
                method = "lsei", dir = dir, ...)

#' @rdname tda_mls
#' @export
tda_mlsei1 <- function(A = NULL, b = NULL, equalities = NULL,
                       equalities_bounds = NULL, inequalities = NULL,
                       inequalities_bounds = NULL, dir = tempfile("tda"), ...)
    .tda_ls_impl(A = A, b = b, equalities = equalities,
                equalities_bounds = equalities_bounds,
                inequalities = inequalities,
                inequalities_bounds = inequalities_bounds,
                method = "active_set", dir = dir, ...)

#' @rdname tda_mls
#' @export
tda_mnls <- function(A, b, equalities = NULL, equalities_bounds = NULL,
                     dir = tempfile("tda"), ...)
    .tda_ls_impl(A = A, b = b, equalities = equalities,
                equalities_bounds = equalities_bounds,
                nonneg = TRUE, dir = dir, ...)

.tda_ls_impl <- function(A = NULL, b = NULL, equalities = NULL,
                         equalities_bounds = NULL,
                         inequalities = NULL,
                         inequalities_bounds = NULL,
                         nonneg = FALSE,
                         method = c("lsei", "active_set"),
                         dir = tempfile("tda"), ...) {
    method <- match.arg(method)
    if (nonneg && !is.null(inequalities))
        stop("`nonneg` and `inequalities` cannot be combined: TDA's ",
             "mnls has no general inequality block", call. = FALSE)
    if (nonneg && method == "active_set")
        stop("`method = \"active_set\"` cannot combine with `nonneg`; ",
             "use `method = \"lsei\"` for `nonneg = TRUE`", call. = FALSE)
    n <- if (!is.null(A)) ncol(rbind(A))
         else if (!is.null(equalities)) ncol(rbind(equalities))
         else if (!is.null(inequalities)) ncol(rbind(inequalities))
         else stop("at least one of `A`, `equalities` or `inequalities` ",
                   "is required", call. = FALSE)
    check_block <- function(mat, rhs, mat_name, rhs_name) {
        if (is.null(mat)) return(NULL)
        mat <- rbind(mat)
        if (ncol(mat) != n)
            stop(sprintf("`%s` must have %d column(s)", mat_name, n),
                 call. = FALSE)
        if (nrow(mat) != length(rhs))
            stop(sprintf("`%s` must have one entry per row of `%s`",
                         rhs_name, mat_name), call. = FALSE)
        mat
    }
    A <- check_block(A, b, "A", "b")
    equalities <- check_block(equalities, equalities_bounds,
                              "equalities", "equalities_bounds")
    inequalities <- check_block(inequalities, inequalities_bounds,
                                "inequalities", "inequalities_bounds")
    me <- if (is.null(equalities)) 0L else nrow(equalities)
    ma <- if (is.null(A)) 0L else nrow(A)
    mi <- if (is.null(inequalities)) 0L else nrow(inequalities)
    S <- rbind(cbind(equalities, equalities_bounds),
              cbind(A, b),
              cbind(inequalities, inequalities_bounds))
    m <- me + ma + mi
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    utils::write.table(S, file.path(dir, "s.mat"),
                       row.names = FALSE, col.names = FALSE)
    cmd_mat <- sprintf("mdef(S,%d,%d) = s.mat;", m, n + 1L)
    if (nonneg)
        cmds <- c("silent = -1;", "mfmt = 24.16;", cmd_mat,
                  sprintf("mnls(S,%d,0,X);", me), "mpr(X) = x.out;")
    else if (method == "active_set")
        cmds <- c("silent = -1;", "mfmt = 24.16;", cmd_mat,
                  sprintf("mlsei1(S,%d,%d,X);", me, mi), "mpr(X) = x.out;")
    else
        cmds <- c("silent = -1;", "mfmt = 24.16;", cmd_mat,
                  sprintf("mlsei(S,%d,%d,X);", me, mi), "mpr(X) = x.out;")
    res <- tda_run(cmds, dir = dir, ...)
    if (any(grepl("contradictory|Cannot satisfy|Inconsistent constraints",
                 res$output)))
        stop("the constraints admit no feasible solution", call. = FALSE)
    if (any(grepl("Exceeded maximal number of iterations|Insufficient accuracy for convergence",
                 res$output)))
        stop("TDA's active-set solver did not converge", call. = FALSE)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run this: ", err[1L], call. = FALSE)
    xv <- NULL
    rank <- NA_integer_
    if (.use_exports()) {
        xm <- res$exports[["mpr.matrix"]]
        if (is.matrix(xm))
            xv <- as.numeric(xm)
        if (!nonneg && method == "lsei") {
            rm_ <- res$exports[["mls.rank"]]
            if (is.matrix(rm_) && length(rm_) == 1L)
                rank <- as.integer(rm_[1L])
        }
    }
    if (is.null(xv)) {
        mx <- tda_file(res, "x.out")
        if (is.null(mx))
            stop("this produced no solution", call. = FALSE)
        xv <- as.numeric(as.matrix(mx)[, 1L])
        if (!nonneg && method == "lsei") {
            rline <- grep("^Rank of left-hand side:", res$output)
            if (length(rline) == 1L)
                rank <- as.integer(sub(".*: *", "", res$output[rline]))
        }
    }
    residual <- if (is.null(A)) NA_real_
               else sqrt(sum((A %*% xv - b)^2))
    list(x = xv, residual = residual, rank = rank)
}

#' Quadratic programming
#'
#' Minimize a convex quadratic subject to box bounds and/or linear
#' equality/inequality constraints -- TDA's \code{mqp}, \code{mqpb},
#' and \code{mqpc} (undocumented in \code{tda.hlp}; identified from
#' \code{qld}'s header, the Powell/Schittkowski active-set solver
#' \code{ZQPCVX} shared with \code{\link{tda_mlsei1}}). The problem
#' solved is
#' \deqn{\min\; d'x + \tfrac12 x'Cx \quad \mathrm{s.t.}\;
#' A_e x = b_e,\; A_i x \ge b_i,\; l \le x \le u}{
#' min d'x + 0.5 x'Cx  s.t.  Ae x = be, Ai x >= bi, l <= x <= u}
#' \code{tda_mqp} is unconstrained (box bounds fixed at +-huge
#' internally by TDA); \code{tda_mqpb} adds box bounds; \code{tda_mqpc}
#' adds general equality/inequality constraints instead of bounds --
#' no single TDA command accepts both bounds and general constraints.
#' The suite checks all three commands against hand-solved instances: an
#' unconstrained case, a box-clipped case, and a mixed
#' equality+inequality case -- see \code{examples/coverage/mqpops2.cf}.
#'
#' @param C the symmetric objective matrix (n x n).
#' @param d the linear objective term, length n.
#' @param lower,upper box bounds, each length n. \code{tda_mqpb} only.
#' @param equalities an equality constraint matrix
#'   \code{Ae x = equalities_bounds}. \code{tda_mqpc} only.
#' @param equalities_bounds the right-hand sides for \code{equalities}.
#' @param inequalities an inequality constraint matrix
#'   \code{Ai x >= inequalities_bounds}. \code{tda_mqpc} only.
#' @param inequalities_bounds the right-hand sides for
#'   \code{inequalities}.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A list: \code{x}, the solution; \code{value}, the objective
#'   \code{d'x + 0.5 x'Cx} at \code{x}, computed directly in R since
#'   this command family prints no value of its own.
#' @family optimization
#' @examples
#' # minimize x1^2 + x2^2 - 2x1 - 4x2 (unconstrained optimum (1, 2))
#' tda_mqp(diag(c(2, 2)), c(-2, -4))
#' tda_mqpb(diag(c(2, 2)), c(-2, -4), lower = c(0, 0), upper = c(0.5, 0.5))
#' @export
tda_mqp <- function(C, d, dir = tempfile("tda"), ...)
    .tda_qp_impl(C = C, d = d, dir = dir, ...)

#' @rdname tda_mqp
#' @export
tda_mqpb <- function(C, d, lower, upper, dir = tempfile("tda"), ...)
    .tda_qp_impl(C = C, d = d, lower = lower, upper = upper, dir = dir, ...)

#' @rdname tda_mqp
#' @export
tda_mqpc <- function(C, d, equalities = NULL, equalities_bounds = NULL,
                     inequalities = NULL, inequalities_bounds = NULL,
                     dir = tempfile("tda"), ...) {
    if (is.null(equalities) && is.null(inequalities))
        stop("`tda_mqpc` needs `equalities` and/or `inequalities`; use ",
             "`tda_mqp` for an unconstrained problem", call. = FALSE)
    .tda_qp_impl(C = C, d = d, equalities = equalities,
                equalities_bounds = equalities_bounds,
                inequalities = inequalities,
                inequalities_bounds = inequalities_bounds, dir = dir, ...)
}

.tda_qp_impl <- function(C, d, lower = NULL, upper = NULL,
                         equalities = NULL,
                         equalities_bounds = NULL,
                         inequalities = NULL,
                         inequalities_bounds = NULL,
                         dir = tempfile("tda"), ...) {
    C <- as.matrix(C)
    n <- nrow(C)
    if (n < 1L || ncol(C) != n)
        stop("`C` must be a square matrix", call. = FALSE)
    if (max(abs(C - t(C))) > 1e-8 * max(1, max(abs(C))))
        stop("`C` must be a symmetric matrix", call. = FALSE)
    if (length(d) != n)
        stop("`d` must have one entry per row of `C`", call. = FALSE)
    has_bounds <- !is.null(lower) || !is.null(upper)
    has_constraints <- !is.null(equalities) || !is.null(inequalities)
    if (has_bounds && has_constraints)
        stop("`lower`/`upper` cannot be combined with `equalities`/",
             "`inequalities`: no single TDA command accepts both",
             call. = FALSE)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    utils::write.table(C, file.path(dir, "c.mat"),
                       row.names = FALSE, col.names = FALSE)
    utils::write.table(matrix(d), file.path(dir, "d.mat"),
                       row.names = FALSE, col.names = FALSE)
    cmds <- c("silent = -1;", "mfmt = 24.16;",
             sprintf("mdef(C,%d,%d) = c.mat;", n, n),
             sprintf("mdef(D,%d,%d) = d.mat;", n, 1L))
    if (has_bounds) {
        if (length(lower) != n || length(upper) != n)
            stop("`lower` and `upper` must each have one entry per ",
                 "row of `C`", call. = FALSE)
        utils::write.table(matrix(lower), file.path(dir, "xl.mat"),
                           row.names = FALSE, col.names = FALSE)
        utils::write.table(matrix(upper), file.path(dir, "xu.mat"),
                           row.names = FALSE, col.names = FALSE)
        cmds <- c(cmds, sprintf("mdef(XL,%d,%d) = xl.mat;", n, 1L),
                 sprintf("mdef(XU,%d,%d) = xu.mat;", n, 1L),
                 "mqpb(C,D,XL,XU,X);")
    }
    else if (has_constraints) {
        check_block <- function(mat, rhs, mat_name, rhs_name) {
            if (is.null(mat)) return(NULL)
            mat <- rbind(mat)
            if (ncol(mat) != n)
                stop(sprintf("`%s` must have %d column(s)", mat_name, n),
                     call. = FALSE)
            if (nrow(mat) != length(rhs))
                stop(sprintf("`%s` must have one entry per row of `%s`",
                             rhs_name, mat_name), call. = FALSE)
            mat
        }
        equalities <- check_block(equalities, equalities_bounds,
                                  "equalities", "equalities_bounds")
        inequalities <- check_block(inequalities, inequalities_bounds,
                                    "inequalities", "inequalities_bounds")
        me <- if (is.null(equalities)) 0L else nrow(equalities)
        Amat <- rbind(equalities, inequalities)
        bvec <- c(equalities_bounds, inequalities_bounds)
        m <- me + (if (is.null(inequalities)) 0L else nrow(inequalities))
        utils::write.table(Amat, file.path(dir, "a.mat"),
                           row.names = FALSE, col.names = FALSE)
        utils::write.table(matrix(bvec), file.path(dir, "b.mat"),
                           row.names = FALSE, col.names = FALSE)
        cmds <- c(cmds, sprintf("mdef(A,%d,%d) = a.mat;", m, n),
                 sprintf("mdef(B,%d,%d) = b.mat;", m, 1L),
                 sprintf("mqpc(C,D,A,B,%d,X);", me))
    }
    else
        cmds <- c(cmds, "mqp(C,D,X);")
    cmds <- c(cmds, "mpr(X) = x.out;")
    res <- tda_run(cmds, dir = dir, ...)
    if (any(grepl("Inconsistent constraints|Exceeded maximal number of iterations|Insufficient accuracy for convergence",
                 res$output)))
        stop("TDA's active-set QP solver could not find a feasible ",
             "solution", call. = FALSE)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run this: ", err[1L], call. = FALSE)
    xv <- NULL
    if (.use_exports()) {
        xm <- res$exports[["mpr.matrix"]]
        if (is.matrix(xm))
            xv <- as.numeric(xm)
    }
    if (is.null(xv)) {
        mx <- tda_file(res, "x.out")
        if (is.null(mx))
            stop("this produced no solution", call. = FALSE)
        xv <- as.numeric(as.matrix(mx)[, 1L])
    }
    list(x = xv, value = as.numeric(d %*% xv + 0.5 * t(xv) %*% C %*% xv))
}

#' 0-1 linear programming
#'
#' Maximize a linear objective over 0-1 variables under linear
#' constraints -- TDA's \code{mlpi} (CACM algorithm 449). The problem
#' solved is
#' \deqn{\max\; c'x + k \quad \mathrm{s.t.}\; A x \ge b,\; x \in \{0,1\}^n}{
#' max c'x + k  subject to  A x >= b,  x in {0,1}^n}
#' with integer coefficients (negative ones are fine). The command's
#' source states no objective; this is what it computes, established
#' against exhaustive search -- see \code{examples/coverage/mlpiops.cf}.
#'
#' @param objective integer coefficients \code{c}, one per variable (at
#'   least two variables).
#' @param constraints the constraint matrix \code{A}, one row per
#'   constraint (at least one), \code{length(objective)} columns.
#' @param bounds the right-hand sides \code{b}, one per constraint row.
#' @param constant \code{k}, added to the reported value (default 0).
#' @param max_solutions how many optimal solutions to return at most.
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A list: \code{value}, the maximum; \code{solutions}, a matrix
#'   with one optimal 0-1 vector per row.
#' @family optimization
#' @examples
#' # maximize 2 x1 - 3 x2 + 4 x3 + x4  s.t.  sum(x) >= 2,
#' # -x1 + 2 x2 + x4 >= 1
#' tda_mlpi(c(2, -3, 4, 1),
#'         rbind(c(1, 1, 1, 1), c(-1, 2, 0, 1)),
#'         bounds = c(2, 1))
#' @export
tda_mlpi <- function(objective, constraints, bounds,
                     constant = 0, max_solutions = 100,
                     dir = tempfile("tda"), ...) {
    constraints <- rbind(constraints)
    if (length(objective) < 2L)
        stop("`objective` needs at least two variables: mlpi's ",
             "minimum", call. = FALSE)
    if (ncol(constraints) != length(objective))
        stop("`constraints` must have one column per objective ",
             "coefficient", call. = FALSE)
    if (nrow(constraints) != length(bounds))
        stop("`bounds` must have one entry per constraint row",
             call. = FALSE)
    A <- rbind(objective, constraints)
    B <- c(constant, bounds)
    if (any(A != round(A)) || any(B != round(B)))
        stop("mlpi works on integer coefficients only", call. = FALSE)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    utils::write.table(A, file.path(dir, "a.mat"), row.names = FALSE,
                       col.names = FALSE)
    utils::write.table(matrix(B), file.path(dir, "b.mat"),
                       row.names = FALSE, col.names = FALSE)
    res <- tda_run(c(sprintf("mdef(A,%d,%d) = a.mat;", nrow(A), ncol(A)),
                     sprintf("mdef(B,%d,1) = b.mat;", nrow(A)),
                     sprintf("mlpi(A,B,%d);", as.integer(max_solutions))),
                   dir = dir, ...)
    if (any(grepl("inconsistent constraints", res$output)))
        stop("the constraints cannot be satisfied by any 0-1 vector",
             call. = FALSE)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not solve this: ", err[1L], call. = FALSE)
    if (.use_exports()) {
        vv <- res$exports[["mlpi.value"]]
        sm <- res$exports[["mlpi.solutions"]]
        if (is.matrix(vv) && length(vv) == 1L && is.matrix(sm))
            return(list(value = as.numeric(vv[1L]),
                        solutions = matrix(as.integer(sm), nrow(sm))))
    }
    vline <- grep("^Value:", res$output)
    nline <- grep("^Number of solutions:", res$output)
    if (length(vline) != 1L || length(nline) != 1L)
        stop("mlpi produced no readable result (no feasible 0-1 vector?)",
             call. = FALSE)
    k <- as.integer(sub(".*: *", "", res$output[nline]))
    sol <- do.call(rbind, lapply(res$output[nline + seq_len(k)],
                                 function(l) as.integer(strsplit(
                                     trimws(l), " +")[[1L]])))
    list(value = as.numeric(sub(".*: *", "", res$output[vline])),
         solutions = sol)
}

#' Boolean function minimization
#'
#' Which combinations of binary conditions produce the outcome? TDA's
#' \code{bfa} builds the truth table of \code{y} over the condition
#' variables, finds all prime implicants (Quine-McCluskey), and selects
#' minimal covers. This is the machinery behind crisp-set QCA-style
#' analyses.
#'
#' @param y the outcome: a 0/1 vector, or a column name in \code{data}.
#' @param conditions the condition variables: a 0/1 matrix or data
#'   frame, or column names in \code{data} (at most 15).
#' @param data optional data frame supplying the columns.
#' @param undefined how rows never observed are treated:
#'   \code{"dont_care"} (default), \code{"true"}, or \code{"false"}.
#' @param algorithm cover selection: \code{"lawler"} (default),
#'   \code{"petrick"}, \code{"lawler2"}, or \code{"none"} for prime
#'   implicants only.
#' @param max_implicants storage cap.
#' @param options a named list of further TDA options for \code{bfa}
#'   (\code{ptab}, \code{prot}, ...).
#' @param dir working directory.
#' @param ... passed to \code{\link{tda_run}}.
#' @return A list: \code{$selections}, a data frame with one row per
#'   term of each minimal cover -- \code{selection} number and one
#'   column per condition holding 1, 0, or NA for "does not matter" --
#'   and \code{$expression}, the printed symbolic forms.
#' @family optimization
#' @examples
#' d <- expand.grid(X1 = 0:1, X2 = 0:1, X3 = 0:1)
#' d$Y <- as.integer(d$X1 & d$X2 | d$X3)
#' tda_boolean_min("Y", c("X1", "X2", "X3"), data = d)$selections
#' @export
tda_boolean_min <- function(y, conditions, data = NULL,
                            undefined = c("dont_care", "true", "false"),
                            algorithm = c("lawler", "petrick", "lawler2",
                                          "none"),
                            max_implicants = NULL, options = list(),
                            dir = tempfile("tda"), ...) {
    yv <- if (is.character(y) && length(y) == 1L) data[[y]] else y
    xv <- if (is.character(conditions)) data[conditions]
          else as.data.frame(conditions)
    if (anyNA(yv) || any(!yv %in% 0:1) ||
        any(vapply(xv, function(v) anyNA(v) || any(!v %in% 0:1), NA)))
        stop("`y` and every condition must be 0/1 without NA",
             call. = FALSE)
    if (length(xv) > 15L)
        stop("bfa handles at most 15 condition variables", call. = FALSE)
    xn <- .tda_names(names(xv))
    d <- stats::setNames(cbind(Y = as.numeric(yv),
                               as.data.frame(lapply(xv, as.numeric))),
                         c("Y", xn))
    o <- c(list(opt = match(match.arg(undefined),
                            c("dont_care", "true", "false")),
                alg = match(match.arg(algorithm),
                            c("none", "lawler", "petrick",
                              "lawler2")) - 1L),
           .tda_extra(options))
    if (!is.null(max_implicants))
        o$max <- format(max_implicants, scientific = FALSE)
    res <- tda_run(c(tda_nvar(d),
                     do.call(tda_block,
                             c(list(name = "bfa"), o,
                               list(rhs = paste(c("Y", xn),
                                                collapse = ","))))),
                   data = d, dir = dir, ...)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run bfa: ", err[1L], call. = FALSE)
    sel <- NULL
    em <- if (.use_exports()) res$exports[["bfa.selection"]]
    if (is.matrix(em) && ncol(em) == length(xn) + 1L) {
        sel <- as.data.frame(em)
        names(sel) <- c("selection", names(xv))
        sel[sel == -1] <- NA
    } else {
        i <- grep("^Final selection", res$output)
        if (length(i)) {
            pat <- res$output[seq.int(i[1L] + 1L, length(res$output))]
            pat <- pat[grepl(sprintf("^[01-]{%d}( \\+ [01-]{%d})*$",
                                     length(xn), length(xn)), trimws(pat))]
            if (length(pat)) {
                terms <- strsplit(trimws(pat[1L]), " \\+ ")[[1L]]
                sel <- do.call(rbind, lapply(terms, function(t) {
                    v <- strsplit(t, "")[[1L]]
                    as.integer(ifelse(v == "-", NA, v))
                }))
                sel <- data.frame(selection = 1L, sel)
                names(sel) <- c("selection", names(xv))
            }
        }
    }
    # the symbolic form sits between "Final selection" and "Final
    # check": the lines there that contain letters (the pattern lines
    # are only 0, 1, - and +)
    i1 <- grep("^Final selection", res$output)
    i2 <- grep("^Final check", res$output)
    ex <- character()
    if (length(i1) && length(i2) && i2[1L] > i1[1L]) {
        w <- res$output[seq.int(i1[1L] + 1L, i2[1L] - 1L)]
        ex <- w[grepl("[A-Za-z]", w)]
    }
    structure(list(call = match.call(), n = length(yv),
                   table = sel, selections = sel,
                   expression = if (length(ex)) ex[length(ex)]
                                else NA_character_,
                   n_implicants = suppressWarnings(as.integer(sub(
                       "^Found ([0-9]+) prime.*", "\\1",
                       grep("^Found [0-9]+ prime", res$output,
                            value = TRUE)[1L]))),
                   run = res),
              class = c("tda_boolean_min", "tda_table"))
}
