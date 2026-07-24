# Thin, named wrappers around individual TDA matrix-language commands.
#
# Each of these is a one-line call into tda_mat(), which already does
# the file-writing/tda_run/export-reading work generically and is
# tested on its own. They exist so that a command someone would
# reasonably look for by its TDA name -- mmul, mtransp, minvs, ... --
# has a matching R function, without needing to know tda_mat()'s
# generic op="..." interface at all. Grouped into a handful of Rd
# pages by family (see @rdname) rather than one page per function.

.tda_scalar <- function(m) if (is.matrix(m) && length(m) == 1L) as.numeric(m) else m
.tda_vec <- function(m) as.numeric(m)

#' Matrix multiplication, transposition, Cholesky decomposition
#'
#' \code{tda_mmul} multiplies any number of matrices in order
#' (\code{A1 \%*\% A2 \%*\% ... \%*\% Ak}) -- TDA's \code{mmul}.
#' \code{tda_mtransp} transposes -- TDA's \code{mtransp}.
#' \code{tda_mchol} returns the lower-triangular Cholesky factor
#' \code{L} such that \code{L \%*\% t(L) == A} (i.e. \code{t(chol(A))}
#' in base R terms; base R's \code{chol} returns the upper
#' factor) -- TDA's \code{mchol}, requiring a positive-definite
#' \code{A}.
#'
#' @param ... for \code{tda_mmul}, two or more numeric matrices to
#'   multiply in order.
#' @param A a numeric matrix.
#' @param dir working directory.
#' @return A numeric matrix.
#' @family matrix algebra
#' @examples
#' tda_mmul(matrix(1:4, 2), diag(2))
#' tda_mtransp(matrix(1:6, 2))
#' tda_mchol(matrix(c(4, 2, 2, 3), 2))
#' @export
tda_mmul <- function(..., dir = tempfile("tda"))
    tda_mat("mmul", ..., out = "R", dir = dir)

#' @rdname tda_mmul
#' @export
tda_mtransp <- function(A, dir = tempfile("tda"))
    tda_mat("mtransp", A, out = "R", dir = dir)

#' @rdname tda_mmul
#' @export
tda_mchol <- function(A, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (nrow(A) != ncol(A))
        stop("`A` must be a square matrix", call. = FALSE)
    if (inherits(try(chol(A), silent = TRUE), "try-error"))
        stop("`A` must be positive definite", call. = FALSE)
    tda_mat("mchol", A, out = "R", dir = dir)
}

#' Matrix inversion
#'
#' Three different things, despite the similar names -- read the one
#' you need. \code{tda_minvs} is the actual matrix inverse of a
#' symmetric positive-definite matrix (TDA's \code{minvs}, via a
#' Cholesky-based algorithm; refuses if not positive definite).
#' \code{tda_minvd} inverts only the diagonal elements of \code{A},
#' element by element -- \emph{not} a matrix inverse at all, and
#' off-diagonal entries are ignored entirely (TDA's \code{minvd}; a
#' zero diagonal entry is silently left as 0 rather than inverted).
#' \code{tda_mginv} is the Moore-Penrose pseudoinverse, requiring at
#' least as many rows as columns and full column rank (TDA's
#' \code{mginv}; prints the pseudorank as a side effect, discarded
#' here -- read it from \code{attr(x, "run")$output} if needed).
#'
#' @param A a numeric matrix: square for \code{tda_minvs}; any shape
#'   for \code{tda_minvd}; at least as many rows as columns for
#'   \code{tda_mginv}.
#' @param dir working directory.
#' @return A numeric matrix.
#' @family matrix algebra
#' @examples
#' A <- matrix(c(4, 2, 2, 3), 2)
#' tda_minvs(A)            # == solve(A) for symmetric positive-definite A
#' tda_minvd(matrix(c(2, 1, 4, 5, 3, 8, 9, 7, 1), 3))  # diag(1/diag(.))
#' tda_mginv(matrix(c(1, 1, 1, 1, 2, 3), 3))
#' @export
tda_minvs <- function(A, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (nrow(A) != ncol(A))
        stop("`A` must be a square matrix", call. = FALSE)
    if (inherits(try(chol(A), silent = TRUE), "try-error"))
        stop("`A` must be positive definite", call. = FALSE)
    tda_mat("minvs", A, out = "R", dir = dir)
}

#' @rdname tda_minvs
#' @export
tda_minvd <- function(A, dir = tempfile("tda"))
    tda_mat("minvd", A, out = "R", dir = dir)

#' @rdname tda_minvs
#' @export
tda_mginv <- function(A, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (nrow(A) < ncol(A))
        stop("`A` must have at least as many rows as columns",
             call. = FALSE)
    tda_mat("mginv", A, out = "R", dir = dir)
}

#' Matrix size and norms
#'
#' Scalar summaries of a matrix -- TDA's \code{mnrow}, \code{mncol},
#' \code{mnorm} (max absolute element), \code{mnorm1} (sum of
#' absolute elements), \code{mnorm2} (Frobenius/Euclidean norm,
#' \code{sqrt(sum(A^2))}), and \code{mtrace} (\code{sum(diag(A))},
#' using \code{min(nrow(A), ncol(A))} diagonal elements for a
#' non-square \code{A}, same as base R's \code{sum(diag(A))}).
#' Each returns a plain number, not a 1x1 matrix.
#'
#' @param A a numeric matrix.
#' @param dir working directory.
#' @return A single number.
#' @family matrix algebra
#' @examples
#' A <- matrix(c(2, 1, 4, 5, 3, 8, 9, 7, 1), 3)
#' tda_mnrow(A); tda_mncol(A)
#' tda_mnorm(A); tda_mnorm1(A); tda_mnorm2(A); tda_mtrace(A)
#' @export
tda_mnrow <- function(A, dir = tempfile("tda"))
    .tda_scalar(tda_mat("mnrow", A, out = "R", dir = dir))

#' @rdname tda_mnrow
#' @export
tda_mncol <- function(A, dir = tempfile("tda"))
    .tda_scalar(tda_mat("mncol", A, out = "R", dir = dir))

#' @rdname tda_mnrow
#' @export
tda_mnorm <- function(A, dir = tempfile("tda"))
    .tda_scalar(tda_mat("mnorm", A, out = "R", dir = dir))

#' @rdname tda_mnrow
#' @export
tda_mnorm1 <- function(A, dir = tempfile("tda"))
    .tda_scalar(tda_mat("mnorm1", A, out = "R", dir = dir))

#' @rdname tda_mnrow
#' @export
tda_mnorm2 <- function(A, dir = tempfile("tda"))
    .tda_scalar(tda_mat("mnorm2", A, out = "R", dir = dir))

#' @rdname tda_mnrow
#' @export
tda_mtrace <- function(A, dir = tempfile("tda"))
    .tda_scalar(tda_mat("mtrace", A, out = "R", dir = dir))

#' Diagonal construction and extraction
#'
#' \code{tda_mdiag} builds an n x n diagonal matrix from a length-n
#' vector -- TDA's \code{mdiag}, matching base R's \code{diag(x)}
#' when \code{x} has length > 1. \code{tda_mdiagd} extracts the
#' diagonal of a matrix as a plain vector, using
#' \code{min(nrow(A), ncol(A))} entries for a non-square \code{A} --
#' TDA's \code{mdiagd}, matching base R's \code{diag(A)}.
#'
#' @param x a numeric vector.
#' @param A a numeric matrix.
#' @param dir working directory.
#' @return \code{tda_mdiag} returns a matrix; \code{tda_mdiagd}
#'   returns a plain numeric vector.
#' @family matrix algebra
#' @examples
#' tda_mdiag(c(5, 7, 9))
#' tda_mdiagd(matrix(1:6, 2, 3))
#' @export
tda_mdiag <- function(x, dir = tempfile("tda"))
    tda_mat("mdiag", matrix(x), out = "R", dir = dir)

#' @rdname tda_mdiag
#' @export
tda_mdiagd <- function(A, dir = tempfile("tda"))
    .tda_vec(tda_mat("mdiagd", A, out = "R", dir = dir))

#' Vectorizing a matrix and its inverse
#'
#' \code{tda_mcvec} stacks the columns of \code{A} into a single
#' vector -- TDA's \code{mcvec}, the classic matrix-algebra
#' \code{vec(A)} operator, matching base R's \code{as.vector(A)}
#' (R matrices are column-major already). \code{tda_mrvec} stacks the
#' rows instead -- TDA's \code{mrvec}, matching
#' \code{as.vector(t(A))}. \code{tda_mivec} is the inverse of
#' \code{tda_mcvec}: given a vector and a row count \code{n}, it
#' rebuilds an n-row matrix by filling column by column -- TDA's
#' \code{mivec}, matching base R's \code{matrix(x, nrow = n)}.
#'
#' @param A a numeric matrix.
#' @param x a numeric vector, whose length must divide evenly by
#'   \code{n}.
#' @param n the number of rows for \code{tda_mivec}'s result.
#' @param dir working directory.
#' @return \code{tda_mcvec}/\code{tda_mrvec} return a plain numeric
#'   vector; \code{tda_mivec} returns a matrix.
#' @family matrix reshaping
#' @examples
#' A <- matrix(1:6, 2, 3)
#' tda_mcvec(A)               # == as.vector(A)
#' tda_mrvec(A)               # == as.vector(t(A))
#' tda_mivec(1:6, 2)          # == matrix(1:6, nrow = 2)
#' @export
tda_mcvec <- function(A, dir = tempfile("tda"))
    .tda_vec(tda_mat("mcvec", A, out = "R", dir = dir))

#' @rdname tda_mcvec
#' @export
tda_mrvec <- function(A, dir = tempfile("tda"))
    .tda_vec(tda_mat("mrvec", A, out = "R", dir = dir))

#' @rdname tda_mcvec
#' @export
tda_mivec <- function(x, n, dir = tempfile("tda")) {
    if (length(x) %% n != 0)
        stop("`length(x)` must be a multiple of `n`", call. = FALSE)
    tda_mat("mivec", matrix(x), matrix(n), out = "R", dir = dir)
}

#' Row and column sums, as vectors or diagonal matrices
#'
#' \code{tda_mrsum}/\code{tda_mcsum} are TDA's \code{mrsum}/
#' \code{mcsum}, matching base R's \code{rowSums}/\code{colSums}.
#' \code{tda_mdrow}/\code{tda_mdcol} return the same sums as a
#' diagonal matrix instead of a vector -- TDA's \code{mdrow}/
#' \code{mdcol}, matching \code{diag(rowSums(A))}/
#' \code{diag(colSums(A))}; useful for building a graph's degree
#' matrix from its adjacency matrix.
#'
#' @param A a numeric matrix.
#' @param dir working directory.
#' @return \code{tda_mrsum}/\code{tda_mcsum} return a plain numeric
#'   vector; \code{tda_mdrow}/\code{tda_mdcol} return a diagonal
#'   matrix.
#' @family matrix reshaping
#' @examples
#' A <- matrix(1:6, 2, 3)
#' tda_mrsum(A); tda_mcsum(A)
#' tda_mdrow(A); tda_mdcol(A)
#' @export
tda_mrsum <- function(A, dir = tempfile("tda"))
    .tda_vec(tda_mat("mrsum", A, out = "R", dir = dir))

#' @rdname tda_mrsum
#' @export
tda_mcsum <- function(A, dir = tempfile("tda"))
    .tda_vec(tda_mat("mcsum", A, out = "R", dir = dir))

#' @rdname tda_mrsum
#' @export
tda_mdrow <- function(A, dir = tempfile("tda"))
    tda_mat("mdrow", A, out = "R", dir = dir)

#' @rdname tda_mrsum
#' @export
tda_mdcol <- function(A, dir = tempfile("tda"))
    tda_mat("mdcol", A, out = "R", dir = dir)

#' Sort or rank the rows of a matrix
#'
#' \code{tda_msort} sorts the rows of \code{X} ascending by the
#' columns named in \code{by}, first column first -- TDA's
#' \code{msort}, matching
#' \code{X[do.call(order, as.data.frame(X[, by])), ]}.
#' \code{tda_msort1} does the same and then drops consecutive
#' duplicate rows (comparing every column, not just \code{by}) --
#' TDA's \code{msort1}, matching \code{unique(sorted_X)}.
#' \code{tda_mrank} -- despite the name -- is \emph{not} a
#' linear-algebra matrix rank (no TDA command returns one as a
#' matrix value; \code{\link{tda_mginv}} prints a pseudorank as a
#' side effect, but never returns it). It returns the sort
#' \emph{permutation} itself: an integer vector \code{p} such that
#' \code{X[p, ]} is sorted by \code{by} -- TDA's \code{mrank},
#' which equals base R's
#' \code{order(X[, by[1]], X[, by[2]], ...)}, not \code{rank()}.
#'
#' @param X a numeric matrix.
#' @param by which columns to sort by, in order (1-based, at least
#'   one).
#' @param dir working directory.
#' @return \code{tda_msort}/\code{tda_msort1} return a matrix;
#'   \code{tda_mrank} returns an integer vector, the same length as
#'   \code{nrow(X)}.
#' @family matrix reshaping
#' @examples
#' X <- rbind(c(30, 1), c(10, 1), c(20, 1), c(40, 1))
#' tda_msort(X, 1)
#' tda_mrank(X, 1)   # == order(X[, 1]), i.e. c(2, 3, 1, 4)
#' @export
tda_msort <- function(X, by, dir = tempfile("tda"))
    tda_mat("msort", X, matrix(by), out = "R", dir = dir)

#' @rdname tda_msort
#' @export
tda_msort1 <- function(X, by, dir = tempfile("tda"))
    tda_mat("msort1", X, matrix(by), out = "R", dir = dir)

#' @rdname tda_msort
#' @export
tda_mrank <- function(X, by, dir = tempfile("tda"))
    as.integer(tda_mat("mrank", X, matrix(by), out = "R", dir = dir))

#' Concatenating matrices
#'
#' \code{tda_mcath} concatenates horizontally -- TDA's \code{mcath},
#' matching \code{cbind}. \code{tda_mcatv} concatenates vertically --
#' TDA's \code{mcatv}, matching \code{rbind}. \code{tda_mcathv} takes
#' the matrix direct sum: the inputs placed as block-diagonal
#' blocks, zero elsewhere -- TDA's \code{mcathv} (all three take two
#' or more matrices).
#'
#' @param ... two or more numeric matrices.
#' @param dir working directory.
#' @return A numeric matrix.
#' @family matrix reshaping
#' @examples
#' A <- matrix(1:4, 2); B <- matrix(5:8, 2)
#' tda_mcath(A, B)    # == cbind(A, B)
#' tda_mcatv(A, B)    # == rbind(A, B)
#' tda_mcathv(A, B)   # block-diagonal direct sum
#' @export
tda_mcath <- function(..., dir = tempfile("tda"))
    tda_mat("mcath", ..., out = "R", dir = dir)

#' @rdname tda_mcath
#' @export
tda_mcatv <- function(..., dir = tempfile("tda"))
    tda_mat("mcatv", ..., out = "R", dir = dir)

#' @rdname tda_mcath
#' @export
tda_mcathv <- function(..., dir = tempfile("tda"))
    tda_mat("mcathv", ..., out = "R", dir = dir)

#' Trim or pad a matrix's edges
#'
#' Deletes or adds whole rows/columns at the edges of \code{A} --
#' TDA's \code{mtrim}. A positive count deletes that many rows or
#' columns from that edge; a negative count adds that many zero
#' rows/columns instead. \code{leading_cols}/\code{trailing_cols}
#' act on the first and last columns, \code{leading_rows}/
#' \code{trailing_rows} on the first and last rows.
#'
#' @param A a numeric matrix.
#' @param leading_cols,trailing_cols,leading_rows,trailing_rows
#'   integer counts (default 0); positive deletes, negative pads
#'   with zeros.
#' @param dir working directory.
#' @return A numeric matrix.
#' @family matrix reshaping
#' @examples
#' A <- matrix(1:4, 2)
#' tda_mtrim(A, leading_cols = 1)    # drop the first column
#' tda_mtrim(A, leading_cols = -1)   # add a zero column in front
#' @export
tda_mtrim <- function(A, leading_cols = 0, leading_rows = 0,
                      trailing_cols = 0, trailing_rows = 0,
                      dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (nrow(A) - leading_rows - trailing_rows <= 0 ||
        ncol(A) - leading_cols - trailing_cols <= 0)
        stop("the result must have at least one row and one column",
             call. = FALSE)
    tda_mat("mtrim", A, matrix(leading_cols), matrix(leading_rows),
           matrix(trailing_cols), matrix(trailing_rows),
           out = "R", dir = dir)
}

#' Eigen decomposition of a symmetric matrix
#'
#' Eigenvalues and eigenvectors of a symmetric matrix -- TDA's
#' \code{mevs}. Eigenvalues come back in \emph{descending} order,
#' same as base R's \code{eigen(A, symmetric = TRUE)}. The
#' eigenvectors come back as the \emph{columns} of \code{vectors},
#' same convention as \code{eigen()$vectors} (column \code{i} satisfies
#' \code{A \%*\% vectors[, i] == values[i] * vectors[, i]}); an
#' initial 2x2 test matrix was checked first and satisfied this
#' equally well as rows, which turned out to be a coincidence of
#' that specific matrix's eigenvector matrix happening to be
#' symmetric, not a real property of TDA's convention -- worth
#' recording since it is exactly the kind of false confirmation a
#' too-simple test case can produce. As with any eigenvector, the
#' sign of each one is arbitrary; TDA's and R's choices need not
#' agree.
#'
#' @param A a symmetric numeric matrix (only the lower triangle is
#'   read).
#' @param dir working directory.
#' @return A list: \code{values}, the eigenvalues in descending
#'   order; \code{vectors}, a matrix whose \emph{columns} are the
#'   corresponding unit-norm eigenvectors.
#' @family matrix algebra
#' @examples
#' A <- matrix(c(4, 1, 0,
#'               1, 3, 1,
#'               0, 1, 2), 3)
#' e <- tda_mevs(A)
#' e$values                      # descending, as eigen(A)$values
#' A %*% e$vectors[, 1] - e$values[1] * e$vectors[, 1]   # zero: an eigenpair
#' @export
tda_mevs <- function(A, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (nrow(A) != ncol(A))
        stop("`A` must be a square matrix", call. = FALSE)
    r <- tda_mat("mevs", A, out = c("E", "V"), dir = dir)
    list(values = as.numeric(r$E), vectors = r$V)
}

#' Singular value decomposition
#'
#' \code{tda_msvd} returns just the singular values -- TDA's
#' \code{msvd}. \code{tda_msvd1} returns the full decomposition
#' \code{A = U \%*\% diag(d) \%*\% t(V)} -- TDA's \code{msvd1},
#' matching base R's \code{svd()} convention exactly (singular
#' vectors as columns of \code{u}/\code{v}). \code{A} must have at
#' least as many rows as columns for either. As with any singular
#' vector, each column's sign is arbitrary and need not match
#' \code{svd()}'s choice.
#'
#' @param A a numeric matrix with at least as many rows as columns.
#' @param dir working directory.
#' @return \code{tda_msvd} returns a plain numeric vector of
#'   singular values. \code{tda_msvd1} returns a list: \code{d}, the
#'   singular values; \code{u}, \code{v}, the singular vectors as
#'   columns.
#' @family matrix algebra
#' @examples
#' A <- matrix(c(1, 3, 5, 2, 4, 7), 3, 2)
#' tda_msvd(A)     # == svd(A)$d
#' tda_msvd1(A)    # == svd(A), up to per-column sign
#' @export
tda_msvd <- function(A, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (nrow(A) < ncol(A))
        stop("`A` must have at least as many rows as columns",
             call. = FALSE)
    .tda_vec(tda_mat("msvd", A, out = "Q", dir = dir))
}

#' @rdname tda_msvd
#' @export
tda_msvd1 <- function(A, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (nrow(A) < ncol(A))
        stop("`A` must have at least as many rows as columns",
             call. = FALSE)
    r <- tda_mat("msvd1", A, out = c("Q", "U", "V"), dir = dir)
    list(d = as.numeric(r$Q), u = r$U, v = r$V)
}

#' Selecting or permuting rows and columns
#'
#' \code{tda_msrow}/\code{tda_mscol} select (and may reorder or
#' repeat) rows/columns by index -- TDA's \code{msrow}/\code{mscol},
#' matching \code{X[rows, ]}/\code{X[, cols]}. \code{tda_mprow}/
#' \code{tda_mpcol}/\code{tda_mpsym} permute rows, columns, or both
#' by the same permutation -- TDA's \code{mprow}/\code{mpcol}/
#' \code{mpsym}, matching \code{A[p, ]}/\code{A[, p]}/\code{A[p, p]}
#' (\code{m_mperm}'s header states the convention,
#' \code{b(i,j)=a(p(i),j)} etc.).
#' \code{mprow}/\code{mpcol} require a permutation of
#' \code{1:nrow(A)}/\code{1:ncol(A)}, so unlike \code{msrow}/
#' \code{mscol} they cannot repeat or drop indices; \code{mpsym}
#' additionally requires \code{A} to be square.
#'
#' @param X,A a numeric matrix.
#' @param rows,cols integer indices (1-based; \code{tda_msrow}/
#'   \code{tda_mscol} allow repeats and omissions).
#' @param p an integer permutation of \code{1:nrow(A)}
#'   (\code{tda_mprow}, \code{tda_mpsym}) or \code{1:ncol(A)}
#'   (\code{tda_mpcol}).
#' @param dir working directory.
#' @return A numeric matrix.
#' @family matrix reshaping
#' @examples
#' A <- matrix(1:9, 3, 3)
#' tda_msrow(A, c(2, 1))       # == A[c(2, 1), ]
#' tda_mprow(A, c(2, 3, 1))    # == A[c(2, 3, 1), ]
#' tda_mpsym(A, c(2, 3, 1))    # == A[c(2, 3, 1), c(2, 3, 1)]
#' @export
tda_msrow <- function(X, rows, dir = tempfile("tda"))
    tda_mat("msrow", X, matrix(rows), out = "R", dir = dir)

#' @rdname tda_msrow
#' @export
tda_mscol <- function(X, cols, dir = tempfile("tda"))
    tda_mat("mscol", X, matrix(cols), out = "R", dir = dir)

#' @rdname tda_msrow
#' @export
tda_mprow <- function(A, p, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (!setequal(p, seq_len(nrow(A))))
        stop("`p` must be a permutation of 1:nrow(A)", call. = FALSE)
    tda_mat("mprow", A, matrix(p), out = "R", dir = dir)
}

#' @rdname tda_msrow
#' @export
tda_mpcol <- function(A, p, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (!setequal(p, seq_len(ncol(A))))
        stop("`p` must be a permutation of 1:ncol(A)", call. = FALSE)
    tda_mat("mpcol", A, matrix(p), out = "R", dir = dir)
}

#' @rdname tda_msrow
#' @export
tda_mpsym <- function(A, p, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (nrow(A) != ncol(A))
        stop("`A` must be a square matrix", call. = FALSE)
    if (!setequal(p, seq_len(nrow(A))))
        stop("`p` must be a permutation of 1:nrow(A)", call. = FALSE)
    tda_mat("mpsym", A, matrix(p), out = "R", dir = dir)
}

#' Aggregating a matrix by row and column groups
#'
#' Sums \code{A}'s entries within each (row-group, column-group)
#' block -- TDA's \code{mag}. \code{row_groups}/\code{col_groups}
#' assign each row/column of \code{A} to a group (any positive
#' integers; groups need not be contiguous or start at 1). The
#' result has one row per distinct value in \code{row_groups} and
#' one column per distinct value in \code{col_groups}, in ascending
#' order of group number.
#'
#' @param A a numeric matrix.
#' @param row_groups integer group index for each row, length
#'   \code{nrow(A)}.
#' @param col_groups integer group index for each column, length
#'   \code{ncol(A)}.
#' @param dir working directory.
#' @return A numeric matrix.
#' @family matrix reshaping
#' @examples
#' A <- matrix(1:12, 3, 4, byrow = TRUE)
#' tda_mag(A, row_groups = c(1, 2, 1), col_groups = c(1, 1, 2, 2))
#' @export
tda_mag <- function(A, row_groups, col_groups, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (length(row_groups) != nrow(A))
        stop("`row_groups` must have one entry per row of `A`",
             call. = FALSE)
    if (length(col_groups) != ncol(A))
        stop("`col_groups` must have one entry per column of `A`",
             call. = FALSE)
    tda_mat("mag", A, matrix(row_groups), t(col_groups),
           out = "R", dir = dir)
}

#' Zero out matrix elements by a threshold test
#'
#' Sets \code{A[i,j]} to 0 wherever the chosen comparison against
#' \code{x} holds, leaving every other element unchanged -- TDA's
#' \code{mnc} (documented in \code{tda.hlp}). \code{test} chooses the
#' comparison: \code{"<="}, \code{"<"}, \code{">="}, \code{">"}
#' (against the raw value), or the same four prefixed \code{"abs"}
#' (against \code{abs(A[i,j])}).
#'
#' @param A a numeric matrix.
#' @param x the threshold.
#' @param test which comparison zeroes an element; see above.
#' @param dir working directory.
#' @return A numeric matrix.
#' @family matrix algebra
#' @examples
#' A <- matrix(c(1, 5, 3, -2, 4, 6), 2, 3)
#' tda_mnc(A, 3, "<=")   # zero out every element <= 3
#' @export
tda_mnc <- function(A, x, test = c("<=", "<", ">=", ">",
                                   "abs<=", "abs<", "abs>=", "abs>"),
                    dir = tempfile("tda")) {
    test <- match.arg(test)
    opt <- match(test, c("<=", "<", ">=", ">", "abs<=", "abs<", "abs>=", "abs>"))
    tda_mat("mnc", A, matrix(x), matrix(opt), out = "B", dir = dir)
}

#' Diagonal square root and its reciprocal
#'
#' Given a square matrix \code{A}, returns a diagonal matrix (zero
#' off-diagonal) built from its diagonal -- TDA's \code{msqrtd}
#' (\code{diag(sqrt(diag(A)))}; refuses a negative diagonal element)
#' and \code{msqrti} (\code{diag(1/sqrt(diag(A)))}; refuses an
#' (almost) zero diagonal element). Off-diagonal entries of \code{A}
#' are ignored entirely, same as \code{\link{tda_minvd}}.
#'
#' @param A a square numeric matrix.
#' @param dir working directory.
#' @return A diagonal numeric matrix.
#' @family matrix algebra
#' @examples
#' tda_msqrtd(matrix(c(4, 1, 1, 9), 2))    # diag(2, 3)
#' tda_msqrti(matrix(c(4, 1, 1, 9), 2))    # diag(0.5, 1/3)
#' @export
tda_msqrtd <- function(A, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (nrow(A) != ncol(A))
        stop("`A` must be a square matrix", call. = FALSE)
    if (any(diag(A) < 0))
        stop("`A` must have a non-negative diagonal", call. = FALSE)
    tda_mat("msqrtd", A, out = "R", dir = dir)
}

#' @rdname tda_msqrtd
#' @export
tda_msqrti <- function(A, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (nrow(A) != ncol(A))
        stop("`A` must be a square matrix", call. = FALSE)
    if (any(diag(A) < 0))
        stop("`A` must have a non-negative diagonal", call. = FALSE)
    if (any(sqrt(diag(A)) <= .Machine$double.eps))
        stop("`A` must have no (near) zero diagonal element",
             call. = FALSE)
    tda_mat("msqrti", A, out = "R", dir = dir)
}

#' Centering, standardizing, and double-centering
#'
#' \code{tda_mcent} subtracts each column's mean -- TDA's
#' \code{mcent}, matching \code{sweep(A, 2, colMeans(A))}.
#' \code{tda_mstand} additionally divides each column by its
#' \emph{population} standard deviation (denominator \code{n}, not
#' base R's \code{sd()}'s \code{n-1}) -- TDA's \code{mstand}.
#' \code{tda_mdcent}
#' double-centers a symmetric matrix -- TDA's \code{mdcent}, the
#' classical (Torgerson) transformation multidimensional scaling
#' uses to turn a matrix of \emph{squared} distances into one whose
#' eigendecomposition gives point coordinates.
#'
#' @param A,D a numeric matrix (\code{D}, for \code{tda_mdcent},
#'   should be symmetric -- typically a distance matrix).
#' @param dir working directory.
#' @return A numeric matrix.
#' @family matrix algebra
#' @examples
#' A <- matrix(c(1, 2, 6, 10, 30, 20), 3, 2)
#' tda_mcent(A)
#' tda_mstand(A)
#' D <- matrix(c(0, 3, 4, 3, 0, 5, 4, 5, 0), 3, 3)
#' tda_mdcent(D)
#' @export
tda_mcent <- function(A, dir = tempfile("tda"))
    tda_mat("mcent", A, out = "R", dir = dir)

#' @rdname tda_mcent
#' @export
tda_mstand <- function(A, dir = tempfile("tda"))
    tda_mat("mstand", A, out = "R", dir = dir)

#' @rdname tda_mcent
#' @export
tda_mdcent <- function(D, dir = tempfile("tda")) {
    D <- as.matrix(D)
    if (nrow(D) != ncol(D))
        stop("`D` must be a square matrix", call. = FALSE)
    tda_mat("mdcent", D, out = "R", dir = dir)
}

#' Cross-product and an arithmetic sequence as a column
#'
#' \code{tda_mcross} is TDA's \code{mcross}, matching base R's
#' \code{crossprod(A)} (\code{t(A) \%*\% A}). \code{tda_mnum} builds
#' the length-\code{n} arithmetic sequence starting at \code{x} with
#' step \code{d} -- TDA's \code{mnum}, matching
#' \code{seq(x, by = d, length.out = n)}.
#'
#' @param A a numeric matrix.
#' @param x the starting value, for \code{tda_mnum}.
#' @param d the step, for \code{tda_mnum}.
#' @param n the length, for \code{tda_mnum}.
#' @param dir working directory.
#' @return \code{tda_mcross} returns a matrix; \code{tda_mnum}
#'   returns a plain numeric vector.
#' @family matrix algebra
#' @examples
#' tda_mcross(matrix(1:6, 3, 2))    # == crossprod(matrix(1:6, 3, 2))
#' tda_mnum(5, 2, 4)                # == seq(5, by = 2, length.out = 4)
#' @export
tda_mcross <- function(A, dir = tempfile("tda"))
    tda_mat("mcross", A, out = "R", dir = dir)

#' @rdname tda_mcross
#' @export
tda_mnum <- function(x, d, n, dir = tempfile("tda")) {
    if (n < 1)
        stop("`n` must be at least 1", call. = FALSE)
    .tda_vec(tda_mat("mnum", matrix(x), matrix(d), matrix(n),
                     out = "R", dir = dir))
}

#' Kronecker product
#'
#' TDA's \code{mkp}, matching base R's \code{kronecker(A, B)}
#' exactly.
#'
#' @param A,B numeric matrices.
#' @param dir working directory.
#' @return A numeric matrix, \code{nrow(A)*nrow(B)} by
#'   \code{ncol(A)*ncol(B)}.
#' @family matrix algebra
#' @examples
#' tda_mkp(matrix(1:4, 2), matrix(5:8, 2))
#' @export
tda_mkp <- function(A, B, dir = tempfile("tda"))
    tda_mat("mkp", A, B, out = "R", dir = dir)

#' Weighted mean of the later elements
#'
#' For each position \code{i}, averages \code{a[j]} weighted by
#' \code{weights[j]} over every later position (\code{tda_mwvec}:
#' every \code{j > i}; \code{tda_mwvec1}: every \code{j} with
#' \code{t[j] > t[i]}, letting an arbitrary order vector \code{t}
#' stand in for position -- \code{tda_mwvec} is the special case
#' \code{t = seq_along(a)}) -- TDA's \code{mwvec}/\code{mwvec1}.
#' Falls back to \code{a[i]} unchanged wherever the weights beyond
#' \code{i} sum to zero (including the last position, which always
#' has nothing after it). No base-R builtin matches this directly.
#'
#' @param a,weights numeric vectors of the same length.
#' @param t for \code{tda_mwvec1}, a numeric order vector the same
#'   length as \code{a}; ties are included on neither side (only
#'   strictly later positions count).
#' @param dir working directory.
#' @return A numeric vector, the same length as \code{a}.
#' @family matrix algebra
#' @examples
#' a <- c(10, 20, 30, 40)
#' w <- c(1, 1, 1, 1)
#' # position 1 averages 20, 30, 40; position 3 averages 40 alone;
#' # the last position has nothing after it and keeps its own value
#' tda_mwvec(a, w)
#' # with an order vector, "later" means a larger t: positions 1 and 2
#' # share t = 1, so both average the values at t = 2 and t = 3
#' tda_mwvec1(a, w, t = c(1, 1, 2, 3))
#' @export
tda_mwvec <- function(a, weights, dir = tempfile("tda")) {
    if (length(a) != length(weights))
        stop("`a` and `weights` must have the same length", call. = FALSE)
    .tda_vec(tda_mat("mwvec", matrix(a), matrix(weights), out = "R",
                     dir = dir))
}

#' @rdname tda_mwvec
#' @export
tda_mwvec1 <- function(a, weights, t, dir = tempfile("tda")) {
    if (length(a) != length(weights) || length(a) != length(t))
        stop("`a`, `weights` and `t` must all have the same length",
             call. = FALSE)
    .tda_vec(tda_mat("mwvec1", matrix(a), matrix(weights), matrix(t),
                     out = "R", dir = dir))
}

#' Scale a matrix by its total sum
#'
#' Divides every element of \code{A} by the sum of all its
#' elements -- TDA's \code{mscal1}, matching \code{A / sum(A)}.
#'
#' @param A a numeric matrix whose elements do not sum to zero.
#' @param dir working directory.
#' @return A numeric matrix.
#' @family matrix algebra
#' @examples
#' A <- matrix(c(1, 2, 3, 4), 2)
#' tda_mscal1(A)   # A / sum(A): the cells now sum to one
#' @export
tda_mscal1 <- function(A, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (sum(A) == 0)
        stop("`A`'s elements must not sum to zero", call. = FALSE)
    tda_mat("mscal1", A, out = "R", dir = dir)
}

#' Invert a permutation
#'
#' TDA's \code{mpinv}: given a permutation \code{p} of \code{1:n},
#' returns \code{q} such that \code{q[p[i]] == i} -- matching base
#' R's \code{order(p)}.
#'
#' @param p an integer permutation of \code{1:length(p)}.
#' @param dir working directory.
#' @return An integer vector, the inverse permutation.
#' @family matrix algebra
#' @examples
#' tda_mpinv(c(3, 1, 4, 2))   # == order(c(3, 1, 4, 2))
#' @export
tda_mpinv <- function(p, dir = tempfile("tda")) {
    if (!setequal(p, seq_along(p)))
        stop("`p` must be a permutation of 1:length(p)", call. = FALSE)
    as.integer(tda_mat("mpinv", matrix(p), out = "Q", dir = dir))
}

#' Build an edge list from an adjacency-style matrix
#'
#' Lists every \code{(i, j)} with \code{A[i, j] >= x} as one row
#' \code{(i, j, A[i, j])} -- TDA's \code{mcel}, scanning row by row.
#' Refuses if no element of \code{A} reaches \code{x} (there would be
#' nothing to return).
#'
#' @param A a square numeric matrix.
#' @param x the threshold; an entry is an edge when it is at least
#'   this large.
#' @param dir working directory.
#' @return A 3-column numeric matrix: \code{i}, \code{j}, and the
#'   value \code{A[i, j]}, one row per edge.
#' @family matrix algebra
#' @examples
#' # a valued adjacency matrix; keep the edges of value 2 or more
#' A <- matrix(c(0, 3, 0,
#'               1, 0, 5,
#'               0, 2, 0), 3, byrow = TRUE)
#' tda_mcel(A, 2)   # (1,2) 3, (2,3) 5, (3,2) 2 -- (2,1) 1 is dropped
#' @export
tda_mcel <- function(A, x, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (nrow(A) != ncol(A))
        stop("`A` must be a square matrix", call. = FALSE)
    if (!any(A >= x))
        stop("at least one element of `A` must reach `x`, or there ",
             "are no edges", call. = FALSE)
    tda_mat("mcel", A, matrix(x), out = "L", dir = dir)
}

#' Permute rows for a zero-free diagonal
#'
#' Finds a row permutation that gives a square matrix a zero-free
#' diagonal wherever one exists (Duff's algorithm, ACM 575) -- TDA's
#' \code{mpz}. \code{p} satisfies \code{A[p[i], i] != 0} for every
#' \code{i} (the exact contract \code{m_rperm}'s header states),
#' and \code{B} is built as the gather \code{B[i, ] == A[p[i], ]}, so
#' \code{diag(B)[i] == A[p[i], i]} by construction. An earlier version
#' of TDA's C code built \code{B} as the opposite (a scatter,
#' \code{B[p[i], ] == A[i, ]}) despite \code{p} itself already
#' satisfying the gather contract -- confirmed wrong on a 3x3 matrix
#' with a unique, forced perfect matching (a pure permutation matrix,
#' one nonzero per row and column): the scatter form returned an
#' all-zero diagonal even though TDA's diagnostic confirmed a
#' complete matching was found. Fixed at the source
#' (\code{m_mpz}, \code{t_matc.c}), not worked around here.
#' If \code{A} is structurally singular, some diagonal entries of
#' \code{B} may still be zero -- check \code{all(diag(B) != 0)} if
#' that matters, rather than assuming a full match always exists.
#'
#' @param A a square numeric matrix.
#' @param dir working directory.
#' @return A list: \code{B}, the permuted matrix; \code{p}, the
#'   permutation, with \code{A[p[i], i] != 0} for every \code{i}
#'   where that is achievable.
#' @family matrix algebra
#' @examples
#' A <- matrix(c(0, 3, 0, 1, 0, 5, 2, 4, 0), 3, 3)
#' r <- tda_mpz(A)
#' r$p
#' all(diag(r$B) != 0)
#' @export
tda_mpz <- function(A, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (nrow(A) != ncol(A))
        stop("`A` must be a square matrix", call. = FALSE)
    r <- tda_mat("mpz", A, out = c("B", "P"), dir = dir)
    list(B = r$B, p = as.integer(r$P))
}

#' Permute to block triangular form
#'
#' Finds a symmetric permutation that puts a square matrix into
#' block lower-triangular (\code{tda_mpbl}) or block upper-triangular
#' (\code{tda_mpbu}) form (Duff & Reid, CACM 529) -- TDA's
#' \code{mpbl}/\code{mpbu}. \code{B[i, j] == A[p[i], p[j]]} for
#' every \code{i, j}; \code{block}
#' gives, for each row/column of \code{B}, which block it belongs to
#' (blocks numbered in the order \code{B} is arranged, so entries
#' between an earlier and a later block are always zero on the side
#' \code{tda_mpbl}/\code{tda_mpbu} promises -- confirmed on an
#' instance with cross-block coupling, not one where every
#' block trivially has size 1).
#'
#' @param A a square numeric matrix.
#' @param dir working directory.
#' @return A list: \code{B}, the permuted, block-triangular matrix;
#'   \code{p}, the permutation (\code{B == A[p, p]}); \code{block},
#'   an integer vector the same length as \code{p} giving each
#'   position's block number.
#' @family matrix algebra
#' @examples
#' A <- matrix(c(1, 2, 0, 0, 0, 1, 0, 0, 0, 0, 1, 3, 0, 0, 0, 1), 4, 4)
#' r <- tda_mpbl(A)
#' r$B
#' @export
tda_mpbl <- function(A, dir = tempfile("tda"))
    .tda_mpb_impl(A, "mpbl", dir)

#' @rdname tda_mpbl
#' @export
tda_mpbu <- function(A, dir = tempfile("tda"))
    .tda_mpb_impl(A, "mpbu", dir)

.tda_mpb_impl <- function(A, op, dir) {
    A <- as.matrix(A)
    n <- nrow(A)
    if (n != ncol(A))
        stop("`A` must be a square matrix", call. = FALSE)
    r <- tda_mat(op, A, out = c("B", "P", "N", "U"), dir = dir)
    p <- as.integer(r$P)
    nblk <- as.integer(r$N)
    starts <- as.integer(r$U)
    block <- rep.int(seq_len(nblk), diff(c(starts, n + 1L)))
    list(B = r$B, p = p, block = block)
}

#' Iterative proportional fitting
#'
#' Adjusts a table \code{A} so its row and column sums match
#' prescribed targets while preserving cross-product ratios (the
#' RAS/IPF algorithm) -- TDA's \code{mpfit}. \code{A} must
#' be non-negative with strictly positive row and column sums, as
#' must the target \code{row_sums}/\code{col_sums}, and their totals
#' should agree (TDA warns but does not refuse otherwise).
#'
#' @param A a non-negative numeric matrix (n x m).
#' @param row_sums target row sums, length n.
#' @param col_sums target column sums, length m.
#' @param max_iter maximum number of iterations.
#' @param eps convergence tolerance on the maximum deviation between
#'   the fitted and target sums.
#' @param dir working directory.
#' @return A list: \code{B}, the fitted table; \code{iterations}, how
#'   many were used; \code{accuracy}, the final deviation reached.
#' @family matrix algebra
#' @examples
#' tda_mpfit(matrix(c(10, 20, 30, 40), 2, 2, byrow = TRUE),
#'          row_sums = c(45, 55), col_sums = c(40, 60))
#' @export
tda_mpfit <- function(A, row_sums, col_sums, max_iter = 100, eps = 1e-6,
                      dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (length(row_sums) != nrow(A))
        stop("`row_sums` must have one entry per row of `A`", call. = FALSE)
    if (length(col_sums) != ncol(A))
        stop("`col_sums` must have one entry per column of `A`",
             call. = FALSE)
    if (max_iter < 1)
        stop("`max_iter` must be at least 1", call. = FALSE)
    if (eps < 0)
        stop("`eps` must be non-negative", call. = FALSE)
    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)
    utils::write.table(A, file.path(dir, "a.mat"),
                       row.names = FALSE, col.names = FALSE)
    utils::write.table(t(col_sums), file.path(dir, "u.mat"),
                       row.names = FALSE, col.names = FALSE)
    utils::write.table(matrix(row_sums), file.path(dir, "v.mat"),
                       row.names = FALSE, col.names = FALSE)
    res <- tda_run(c("silent = -1;", "mfmt = 24.16;",
                     sprintf("mdef(A,%d,%d) = a.mat;", nrow(A), ncol(A)),
                     sprintf("mdef(U,1,%d) = u.mat;", ncol(A)),
                     sprintf("mdef(V,%d,1) = v.mat;", nrow(A)),
                     sprintf("mdef(ITER,1,1) = %d;", as.integer(max_iter)),
                     sprintf("mdef(EPS,1,1) = %.15g;", eps),
                     "mpfit(A,U,V,ITER,EPS,B);", "mpr(B) = b.out;"),
                   dir = dir)
    if (any(grepl("error in number of iterations|error in epsilon",
                 res$output)))
        stop("TDA rejected `max_iter` or `eps`", call. = FALSE)
    err <- grep("^Error", res$output, value = TRUE)
    if (length(err))
        stop("TDA could not run mpfit: ", err[1L], call. = FALSE)
    B <- NULL
    if (.use_exports()) {
        bm <- res$exports[["mpr.matrix"]]
        if (is.matrix(bm))
            B <- bm
    }
    if (is.null(B)) {
        mb <- tda_file(res, "b.out")
        if (is.null(mb))
            stop("mpfit produced no result", call. = FALSE)
        B <- as.matrix(mb)
    }
    iterline <- grep("^Number of iterations:", res$output)
    iterations <- NA_integer_
    accuracy <- NA_real_
    if (length(iterline) == 1L) {
        m <- regmatches(res$output[iterline],
                        regexec("iterations: ([0-9]+)\\. Final accuracy: ([^ ]+)",
                                res$output[iterline]))[[1]]
        if (length(m) == 3L) {
            iterations <- as.integer(m[2L])
            accuracy <- as.numeric(m[3L])
        }
    }
    list(B = B, iterations = iterations, accuracy = accuracy)
}

#' Iterate a Leslie matrix
#'
#' Projects an age-structured population forward using a Leslie
#' matrix -- TDA's \code{mpit}/\code{mpit1}. Row \code{i} of
#' \code{fertility_survival} holds \code{(fertility[i], survival[i])}:
#' \code{fertility[i]} is the birth rate contributed by age class
#' \code{i} (used only for class 1 of the next generation);
#' \code{survival[i]} is the fraction of age class \code{i} that
#' survives into age class \code{i + 1} (\code{survival[n]}, the
#' last row, is unused). \code{tda_mpit1} adds a constant vector
#' \code{immigration} to every age class after each projection step.
#'
#' @param fertility_survival an n x 2 numeric matrix, columns
#'   \code{(fertility, survival)} as above.
#' @param population the initial population vector, length n.
#' @param immigration for \code{tda_mpit1}, a constant vector added
#'   to the population after each step, length n.
#' @param iterations how many projection steps to take.
#' @param dir working directory.
#' @return A matrix with \code{iterations + 1} rows (the initial
#'   population, then one row per step) and \code{n} columns.
#' @family matrix algebra
#' @examples
#' F <- rbind(c(0, 0), c(2, 0.5), c(1, 0.3))
#' tda_mpit(F, c(100, 50, 20), iterations = 3)
#' @export
tda_mpit <- function(fertility_survival, population, iterations,
                     dir = tempfile("tda")) {
    FS <- as.matrix(fertility_survival)
    n <- nrow(FS)
    if (ncol(FS) != 2)
        stop("`fertility_survival` must have exactly 2 columns",
             call. = FALSE)
    if (length(population) != n)
        stop("`population` must have one entry per row of ",
             "`fertility_survival`", call. = FALSE)
    tda_mat("mpit", FS, matrix(population), matrix(iterations),
           out = "R", dir = dir)
}

#' @rdname tda_mpit
#' @export
tda_mpit1 <- function(fertility_survival, population, immigration,
                      iterations, dir = tempfile("tda")) {
    FS <- as.matrix(fertility_survival)
    n <- nrow(FS)
    if (ncol(FS) != 2)
        stop("`fertility_survival` must have exactly 2 columns",
             call. = FALSE)
    if (length(population) != n)
        stop("`population` must have one entry per row of ",
             "`fertility_survival`", call. = FALSE)
    if (length(immigration) != n)
        stop("`immigration` must have one entry per row of ",
             "`fertility_survival`", call. = FALSE)
    tda_mat("mpit1", FS, matrix(population), matrix(immigration),
           matrix(iterations), out = "R", dir = dir)
}

#' Kemeny distance between rankings
#'
#' Pairwise Kemeny (Kendall-tau-with-ties) distance between every
#' pair of rows of \code{A}, each row a ranking (or score) over the
#' same set of items -- TDA's \code{mkmet}. For two rankings, each
#' pair of items contributes 0 if both rankings order (or tie) the
#' pair the same way, 1 if one ranking ties the pair and the other
#' doesn't, or 2 if the rankings strictly disagree on the pair's
#' order.
#'
#' @param A a numeric matrix, one ranking per row.
#' @param dir working directory.
#' @return A symmetric numeric matrix with zero diagonal, the Kemeny
#'   distance between every pair of rows.
#' @family matrix algebra
#' @examples
#' # four rankings of the same four items: b swaps the last two of a,
#' # c reverses a, d ties items in pairs
#' R <- rbind(a = c(1, 2, 3, 4),
#'            b = c(1, 2, 4, 3),
#'            c = c(4, 3, 2, 1),
#'            d = c(1, 1, 2, 2))
#' tda_mkmet(R)   # a-b: 2 (one pair reversed); a-c: 12 (all six pairs)
#' @export
tda_mkmet <- function(A, dir = tempfile("tda"))
    tda_mat("mkmet", A, out = "D", dir = dir)

#' Product-limit (Kaplan-Meier) distribution and jumps
#'
#' The product-limit (Kaplan-Meier) estimate of the CDF for a vector
#' of possibly right-censored event times -- TDA's \code{mple}, via
#' Efron's redistribute-to-the-right algorithm, applied one
#' observation at a time in sorted order (events before censoring at
#' the same time). When no two events fall at exactly the same time,
#' this matches the standard product-limit CDF exactly -- confirmed
#' against \code{survival::survfit} across random instances (events
#' may still tie with a censored observation; only ties \emph{among
#' events} are excluded). When two or more events tie,
#' processing them one at a time rather than as a single simultaneous
#' risk-set reduction means only the \emph{last} observation
#' processed within that tied group reaches the value standard
#' Kaplan-Meier software reports for that time -- earlier ones in the
#' same group show smaller, order-dependent intermediate values.
#' Two events tied at \code{time = 2} out of 3
#' total observations give \code{F = (1/3, 2/3, 1)} here, where
#' standard Kaplan-Meier reports \code{F(2) = 2/3} for both tied
#' observations. This is a property of the algorithm as TDA
#' implements it, not a bug and not this wrapper's choice to make.
#' Separately: \code{F} always reaches exactly 1 at the
#' highest-time observation, \emph{even when that observation is
#' censored} -- redistribution only moves probability mass to later
#' observations, never destroys it, and the highest-time observation
#' has nothing later to redistribute its share to, so all of the
#' mass ends up accounted for there by construction. Standard
#' Kaplan-Meier software instead leaves the survival curve at
#' whatever value it last reached, undefined beyond a censored tail.
#' \code{time = c(1, 2, 3)},
#' \code{censored = c(0, 0, 1)} gives \code{F = (1/3, 2/3, 1)} here.
#' This is the low-level building block behind TDA's
#' \code{\link{tda_ple}}, which most users want instead -- this one
#' works directly on plain vectors, with no grouping, formula
#' interface, or standard errors.
#'
#' @param time a numeric vector of event/censoring times.
#' @param censored a numeric vector the same length as \code{time}:
#'   1 where the observation is censored, 0 where it is an observed
#'   event.
#' @param dir working directory.
#' @return A list, both the same length as \code{time}, in the
#'   original (unsorted) order: \code{F}, the estimated CDF at each
#'   observation; \code{jumps}, the increase in \code{F} at each
#'   observation (0 for a censored observation).
#' @family matrix algebra
#' @examples
#' # six durations, two of them censored (still running at 3 and 7)
#' t <- c(2, 3, 3, 5, 7, 8)
#' cen <- c(0, 0, 1, 0, 1, 0)
#' pl <- tda_mple(t, cen)
#' pl$F       # the Kaplan-Meier distribution function at each observation
#' pl$jumps   # its increase there; 0 at the censored ones
#' @export
tda_mple <- function(time, censored, dir = tempfile("tda")) {
    if (length(time) != length(censored))
        stop("`time` and `censored` must have the same length",
             call. = FALSE)
    r <- tda_mat("mple", matrix(time), matrix(censored),
                out = c("F", "D"), dir = dir)
    list(F = as.numeric(r$F), jumps = as.numeric(r$D))
}

#' Eigen decomposition of a general (possibly non-symmetric) matrix
#'
#' Eigenvalues and eigenvectors of a square matrix that need not be
#' symmetric, so eigenvalues and eigenvectors may be complex -- TDA's
#' \code{mev}. \code{t_mat.c}'s dispatch comment understates this
#' command's arity (\code{mev(A,ER,EI,EV)}, four arguments); what
#' \code{m_mev} actually parses is five (\code{mev(A,ER,EI,EVR,EVI)}),
#' splitting the eigenvector matrix into separate real and imaginary
#' parts, as \code{eigen1}'s header lists them. Column
#' \code{i} of the result, \code{EVR[, i] + 1i * EVI[, i]}, is the
#' right eigenvector for eigenvalue \code{ER[i] + 1i * EI[i]}, i.e.
#' \code{A \%*\% v == lambda * v} -- the same convention and ordering
#' as base R's \code{eigen(A)}, confirmed on two matrices, one
#' with real eigenvalues and one with an asymmetric complex
#' pair (not simply plus/minus of each other, which the first,
#' simpler test matrix happened to have and which turned out to mask
#' a construction mistake in an earlier verification attempt --
#' recorded here because it is exactly the kind of false confirmation
#' a too-simple or too-symmetric test case can produce, the same
#' lesson \code{\link{tda_mevs}} needed).
#'
#' @param A a square numeric matrix.
#' @param dir working directory.
#' @return A list: \code{values}, a complex vector of eigenvalues;
#'   \code{vectors}, a complex matrix whose columns are the
#'   corresponding right eigenvectors.
#' @details
#' TDA's \code{mev(A, ER, EI, EVR, EVI)} requires all four
#' output operands (none are optional); the wrapper supplies
#' them and folds the real and imaginary parts into complex
#' \code{$values} and \code{$vectors}.
#' @family matrix algebra
#' @examples
#' tda_mev(matrix(c(0, -1, 1, 0), 2, 2))   # eigenvalues +-i
#' @export
tda_mev <- function(A, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (nrow(A) != ncol(A))
        stop("`A` must be a square matrix", call. = FALSE)
    r <- tda_mat("mev", A, out = c("ER", "EI", "EVR", "EVI"), dir = dir)
    list(values = complex(real = as.numeric(r$ER), imaginary = as.numeric(r$EI)),
        vectors = r$EVR + 1i * r$EVI,
        er = as.numeric(r$ER), ei = as.numeric(r$EI),
        evr = r$EVR, evi = r$EVI)
}

#' Check an ordering against row/column dominance
#'
#' Flags each row/column index of a square matrix \code{A} where the
#' current index order fails a running dominance test against its
#' mirror image -- TDA's \code{mch} (undocumented beyond "check row
#' and column sum conditions" in its own header; the rule below
#' is the loop body's, not a named textbook procedure). For each
#' index \code{i}, scanning outward
#' (\code{j} from \code{i+1} to \code{n}, then separately from
#' \code{i-1} down to \code{1}), \code{A[i, j]} accumulates into a
#' running total \code{z} and its mirror \code{A[j, i]} accumulates
#' into a running total \code{s}; index \code{i} is flagged the first
#' time \code{z} falls behind \code{s} by more than a small tolerance.
#' In effect: is \code{A[i, ]}'s share of an outward pair always at
#' least as large, cumulatively, as the corresponding mirrored
#' entries in \code{A[, i]} -- a diagnostic for whether the current
#' \code{1:n} ordering already respects that dominance, without
#' finding a better ordering the way \code{\link{tda_mpz}}/
#' \code{\link{tda_mpbl}} do.
#'
#' @param A a square numeric matrix.
#' @param dir working directory.
#' @return An integer vector, length \code{nrow(A)}: 1 where that
#'   index fails the test, 0 where it passes.
#' @family matrix algebra
#' @examples
#' # A[i, j] = 5 means i strongly precedes j.  Ordered 1 < 2 < 3 the
#' # rows dominate their columns and every index passes (0); the same
#' # matrix with the order reversed fails everywhere (1).
#' A <- matrix(c(0, 5, 5,
#'               1, 0, 5,
#'               1, 1, 0), 3, byrow = TRUE)
#' tda_mch(A)
#' tda_mch(A[3:1, 3:1])
#' @export
tda_mch <- function(A, dir = tempfile("tda")) {
    A <- as.matrix(A)
    if (nrow(A) != ncol(A))
        stop("`A` must be a square matrix", call. = FALSE)
    as.integer(tda_mat("mch", A, out = "B", dir = dir))
}

#' Distribution function bounds and estimate for interval-censored data
#'
#' Given a sample of intervals \code{[lower_i, upper_i]} (each
#' interval containing an unobserved exact value), estimates the
#' distribution function of that unobserved value -- TDA's \code{midf}
#' (undocumented beyond "calculates distribution function DM, and
#' lower/upper bounds" in its own header, which also understates the
#' command's arity: it lists three outputs, \code{midf(XL,XU,DL,DU,DM)},
#' but \code{m_midf} actually returns four,
#' \code{midf(XL,XU,PT,DL,DU,DM)} -- the breakpoints themselves,
#' \code{PT}, are a result too). Evaluated at every distinct
#' interval endpoint (the
#' \dQuote{induced partition} \code{breakpoints}): \code{lower_bound}
#' counts only intervals entirely below that point (a sure lower bound
#' on the true CDF); \code{upper_bound} counts every interval that
#' could be below it (a sure upper bound); \code{cdf} is a
#' point estimate in between, assuming each interval's unobserved
#' value is uniformly distributed within it. Every interval must have
#' positive width (\code{lower < upper}, swapped automatically if
#' given in the other order). Confirmed cell-for-cell by hand on a
#' 3-interval example.
#'
#' @param lower,upper numeric vectors, the same length, with
#'   \code{lower[i] != upper[i]} for every \code{i}.
#' @param dir working directory.
#' @return A list: \code{breakpoints}, the sorted distinct interval
#'   endpoints; \code{lower_bound}, \code{upper_bound}, \code{cdf},
#'   each the same length as \code{breakpoints}.
#' @family matrix algebra
#' @examples
#' # five observations known only to intervals: [0,2], [1,3], [2,4], [3,5], [1,2]
#' lo <- c(0, 1, 2, 3, 1)
#' up <- c(2, 3, 4, 5, 2)
#' # at each endpoint: the share of intervals surely below it (lower
#' # bound), possibly below it (upper bound), and the estimate between
#' tda_midf(lo, up)
#' @export
tda_midf <- function(lower, upper, dir = tempfile("tda")) {
    if (length(lower) != length(upper))
        stop("`lower` and `upper` must have the same length", call. = FALSE)
    if (any(lower == upper))
        stop("every interval must have `lower != upper`", call. = FALSE)
    r <- tda_mat("midf", matrix(lower), matrix(upper),
                out = c("PT", "DL", "DU", "DM"), dir = dir)
    list(breakpoints = as.numeric(r$PT), lower_bound = as.numeric(r$DL),
        upper_bound = as.numeric(r$DU), cdf = as.numeric(r$DM))
}

#' Uniform-model CDF for interval-censored data, at each observation
#'
#' The same uniform-within-interval distribution function estimate as
#' \code{\link{tda_midf}}'s \code{cdf}, but evaluated at each
#' observation's \code{lower[i]} rather than at every distinct
#' breakpoint -- TDA's \code{midf1}.
#'
#' @param lower,upper numeric vectors, the same length, with
#'   \code{lower[i] != upper[i]} for every \code{i}.
#' @param dir working directory.
#' @return A numeric vector, the same length as \code{lower}, in the
#'   original (unsorted) order.
#' @family matrix algebra
#' @examples
#' lo <- c(0, 1, 2, 3, 1)
#' up <- c(2, 3, 4, 5, 2)
#' # the estimated distribution function at each interval's own position
#' tda_midf1(lo, up)
#' @export
tda_midf1 <- function(lower, upper, dir = tempfile("tda")) {
    if (length(lower) != length(upper))
        stop("`lower` and `upper` must have the same length", call. = FALSE)
    if (any(lower == upper))
        stop("every interval must have `lower != upper`", call. = FALSE)
    as.numeric(tda_mat("midf1", matrix(lower), matrix(upper), out = "F",
                       dir = dir))
}

#' Restricted mean for interval-censored data, at each observation
#'
#' For each observation, the conditional mean of the unobserved value
#' given that it exceeds that observation's \code{lower[i]} --
#' \eqn{E[X \mid X > lower_i]}{E[X | X > lower[i]]} -- estimated by
#' integrating a piecewise-uniform density built from
#' \code{\link{tda_midf}}'s CDF estimate over the induced partition --
#' TDA's \code{midf2}: it builds the CDF at every
#' breakpoint, then integrates the resulting piecewise-constant density
#' from each observation's value onward, normalizing by the
#' survival probability there.
#'
#' @param lower,upper numeric vectors, the same length, with
#'   \code{lower[i] != upper[i]} for every \code{i}.
#' @param dir working directory.
#' @return A numeric vector, the same length as \code{lower}, in the
#'   original (unsorted) order.
#' @family matrix algebra
#' @examples
#' lo <- c(0, 1, 2, 3, 1)
#' up <- c(2, 3, 4, 5, 2)
#' # the expected value of each observation given that it lies in its
#' # interval, under the estimated distribution
#' tda_midf2(lo, up)
#' @export
tda_midf2 <- function(lower, upper, dir = tempfile("tda")) {
    if (length(lower) != length(upper))
        stop("`lower` and `upper` must have the same length", call. = FALSE)
    if (any(lower == upper))
        stop("every interval must have `lower != upper`", call. = FALSE)
    as.numeric(tda_mat("midf2", matrix(lower), matrix(upper), out = "F",
                       dir = dir))
}

#' Average nearby endpoints for interval-censored data
#'
#' For each interval \code{i}, averages the lower endpoints of every
#' interval (itself included) whose lower endpoint falls within
#' \code{[lower[i], upper[i]]}, and separately averages the upper
#' endpoints falling in that same range -- TDA's \code{midf3}. A
#' local smoothing of nearby endpoints, not a distribution estimate
#' like \code{\link{tda_midf}}/\code{\link{tda_midf1}}/
#' \code{\link{tda_midf2}}.
#'
#' @param lower,upper numeric vectors, the same length.
#' @param dir working directory.
#' @return A list, both the same length as \code{lower}, in the
#'   original order: \code{lower_avg}, the average nearby lower
#'   endpoint; \code{upper_avg}, the average nearby upper endpoint.
#' @family matrix algebra
#' @examples
#' lo <- c(0, 1, 2, 3, 1)
#' up <- c(2, 3, 4, 5, 2)
#' # for each interval, the mean lower and upper endpoint of the
#' # intervals overlapping it
#' tda_midf3(lo, up)
#' @export
tda_midf3 <- function(lower, upper, dir = tempfile("tda")) {
    if (length(lower) != length(upper))
        stop("`lower` and `upper` must have the same length", call. = FALSE)
    r <- tda_mat("midf3", matrix(lower), matrix(upper),
                out = c("XL1", "XU1"), dir = dir)
    list(lower_avg = as.numeric(r$XL1), upper_avg = as.numeric(r$XU1))
}
