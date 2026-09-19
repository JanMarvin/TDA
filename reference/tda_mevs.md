# Eigen decomposition of a symmetric matrix

Eigenvalues and eigenvectors of a symmetric matrix – TDA's `mevs`.
Eigenvalues come back in *descending* order, same as base R's
`eigen(A, symmetric = TRUE)`. The eigenvectors come back as the
*columns* of `vectors`, same convention as `eigen()$vectors` – confirmed
by hand on a 3x3 instance (column `i` satisfies
`A %*% vectors[, i] == values[i] * vectors[, i]`); an initial 2x2 test
matrix was checked first and satisfied this equally well as rows, which
turned out to be a coincidence of that specific matrix's eigenvector
matrix happening to be symmetric, not a real property of TDA's
convention – worth recording since it is exactly the kind of false
confirmation a too-simple test case can produce. As with any
eigenvector, the sign of each one is arbitrary; TDA's and R's choices
need not agree.

## Usage

``` r
tda_mevs(A, dir = tempfile("tda"))
```

## Arguments

- A:

  a symmetric numeric matrix (only the lower triangle is read).

- dir:

  working directory.

## Value

A list: `values`, the eigenvalues in descending order; `vectors`, a
matrix whose *columns* are the corresponding unit-norm eigenvectors.

## See also

Other matrix algebra: [`tda_mcel()`](tda_mcel.md),
[`tda_mcent()`](tda_mcent.md), [`tda_mch()`](tda_mch.md),
[`tda_mcross()`](tda_mcross.md), [`tda_mdiag()`](tda_mdiag.md),
[`tda_mev()`](tda_mev.md), [`tda_midf()`](tda_midf.md),
[`tda_midf1()`](tda_midf1.md), [`tda_midf2()`](tda_midf2.md),
[`tda_midf3()`](tda_midf3.md), [`tda_minvs()`](tda_minvs.md),
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
tda_mevs(matrix(c(2, 1, 1, 2), 2))
#> $values
#> [1] 3 1
#> 
#> $vectors
#>           [,1]       [,2]
#> [1,] 0.7071068  0.7071068
#> [2,] 0.7071068 -0.7071068
#> 
```
