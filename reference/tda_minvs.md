# Matrix inversion

Three different things, despite the similar names – read the one you
need. `tda_minvs` is the actual matrix inverse of a symmetric
positive-definite matrix (TDA's `minvs`, via a Cholesky-based algorithm;
refuses if not positive definite). `tda_minvd` inverts only the diagonal
elements of `A`, element by element – *not* a matrix inverse at all, and
off-diagonal entries are ignored entirely (TDA's `minvd`; a zero
diagonal entry is silently left as 0 rather than inverted). `tda_mginv`
is the Moore-Penrose pseudoinverse, requiring at least as many rows as
columns and full column rank (TDA's `mginv`; prints the pseudorank as a
side effect, discarded here – read it from `attr(x, "run")$output` if
needed).

## Usage

``` r
tda_minvs(A, dir = tempfile("tda"))

tda_minvd(A, dir = tempfile("tda"))

tda_mginv(A, dir = tempfile("tda"))
```

## Arguments

- A:

  a numeric matrix: square for `tda_minvs`; any shape for `tda_minvd`;
  at least as many rows as columns for `tda_mginv`.

- dir:

  working directory.

## Value

A numeric matrix.

## See also

Other matrix algebra: [`tda_mcel()`](tda_mcel.md),
[`tda_mcent()`](tda_mcent.md), [`tda_mch()`](tda_mch.md),
[`tda_mcross()`](tda_mcross.md), [`tda_mdiag()`](tda_mdiag.md),
[`tda_mev()`](tda_mev.md), [`tda_mevs()`](tda_mevs.md),
[`tda_midf()`](tda_midf.md), [`tda_midf1()`](tda_midf1.md),
[`tda_midf2()`](tda_midf2.md), [`tda_midf3()`](tda_midf3.md),
[`tda_mkmet()`](tda_mkmet.md), [`tda_mkp()`](tda_mkp.md),
[`tda_mmul()`](tda_mmul.md), [`tda_mnc()`](tda_mnc.md),
[`tda_mnrow()`](tda_mnrow.md), [`tda_mpbl()`](tda_mpbl.md),
[`tda_mpfit()`](tda_mpfit.md), [`tda_mpinv()`](tda_mpinv.md),
[`tda_mpit()`](tda_mpit.md), [`tda_mple()`](tda_mple.md),
[`tda_mpz()`](tda_mpz.md), [`tda_mscal1()`](tda_mscal1.md),
[`tda_msqrtd()`](tda_msqrtd.md), [`tda_msvd()`](tda_msvd.md),
[`tda_mwvec()`](tda_mwvec.md)

## Examples

``` r
A <- matrix(c(4, 2, 2, 3), 2)
tda_minvs(A)            # == solve(A) for symmetric positive-definite A
#>        [,1]  [,2]
#> [1,]  0.375 -0.25
#> [2,] -0.250  0.50
tda_minvd(matrix(c(2, 1, 4, 5, 3, 8, 9, 7, 1), 3))  # diag(1/diag(.))
#>      [,1]      [,2] [,3]
#> [1,]  0.5 0.0000000    0
#> [2,]  0.0 0.3333333    0
#> [3,]  0.0 0.0000000    1
tda_mginv(matrix(c(1, 1, 1, 1, 2, 3), 3))
#>           [,1]         [,2]       [,3]
#> [1,]  1.333333 3.333333e-01 -0.6666667
#> [2,] -0.500000 5.551115e-17  0.5000000
```
