# Permute rows for a zero-free diagonal

Finds a row permutation that gives a square matrix a zero-free diagonal
wherever one exists (Duff's algorithm, ACM 575) – TDA's `mpz`. `p`
satisfies `A[p[i], i] != 0` for every `i` (the exact contract
`m_rperm`'s header states), and `B` is built as the gather
`B[i, ] == A[p[i], ]`, so `diag(B)[i] == A[p[i], i]` by construction. An
earlier version of TDA's C code built `B` as the opposite (a scatter,
`B[p[i], ] == A[i, ]`) despite `p` itself already satisfying the gather
contract – confirmed wrong on a 3x3 matrix with a unique, forced perfect
matching (a pure permutation matrix, one nonzero per row and column):
the scatter form returned an all-zero diagonal even though TDA's
diagnostic confirmed a complete matching was found. Fixed at the source
(`m_mpz`, `t_matc.c`), not worked around here. If `A` is structurally
singular, some diagonal entries of `B` may still be zero – check
`all(diag(B) != 0)` if that matters, rather than assuming a full match
always exists.

## Usage

``` r
tda_mpz(A, dir = tempfile("tda"))
```

## Arguments

- A:

  a square numeric matrix.

- dir:

  working directory.

## Value

A list: `B`, the permuted matrix; `p`, the permutation, with
`A[p[i], i] != 0` for every `i` where that is achievable.

## See also

Other matrix algebra: [`tda_mcel()`](tda_mcel.md),
[`tda_mcent()`](tda_mcent.md), [`tda_mch()`](tda_mch.md),
[`tda_mcross()`](tda_mcross.md), [`tda_mdiag()`](tda_mdiag.md),
[`tda_mev()`](tda_mev.md), [`tda_mevs()`](tda_mevs.md),
[`tda_midf()`](tda_midf.md), [`tda_midf1()`](tda_midf1.md),
[`tda_midf2()`](tda_midf2.md), [`tda_midf3()`](tda_midf3.md),
[`tda_minvs()`](tda_minvs.md), [`tda_mkmet()`](tda_mkmet.md),
[`tda_mkp()`](tda_mkp.md), [`tda_mmul()`](tda_mmul.md),
[`tda_mnc()`](tda_mnc.md), [`tda_mnrow()`](tda_mnrow.md),
[`tda_mpbl()`](tda_mpbl.md), [`tda_mpfit()`](tda_mpfit.md),
[`tda_mpinv()`](tda_mpinv.md), [`tda_mpit()`](tda_mpit.md),
[`tda_mple()`](tda_mple.md), [`tda_mscal1()`](tda_mscal1.md),
[`tda_msqrtd()`](tda_msqrtd.md), [`tda_msvd()`](tda_msvd.md),
[`tda_mwvec()`](tda_mwvec.md)

## Examples

``` r
A <- matrix(c(0, 3, 0, 1, 0, 5, 2, 4, 0), 3, 3)
r <- tda_mpz(A)
r$p
#> [1] 2 3 1
all(diag(r$B) != 0)
#> [1] TRUE
```
