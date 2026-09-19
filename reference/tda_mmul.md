# Matrix multiplication, transposition, Cholesky decomposition

`tda_mmul` multiplies any number of matrices in order
(`A1 %*% A2 %*% ... %*% Ak`) – TDA's `mmul`. `tda_mtransp` transposes –
TDA's `mtransp`. `tda_mchol` returns the lower-triangular Cholesky
factor `L` such that `L %*% t(L) == A` (i.e. `t(chol(A))` in base R
terms; base R's `chol` returns the upper factor) – TDA's `mchol`,
requiring a positive-definite `A`.

## Usage

``` r
tda_mmul(..., dir = tempfile("tda"))

tda_mtransp(A, dir = tempfile("tda"))

tda_mchol(A, dir = tempfile("tda"))
```

## Arguments

- ...:

  for `tda_mmul`, two or more numeric matrices to multiply in order.

- dir:

  working directory.

- A:

  a numeric matrix.

## Value

A numeric matrix.

## See also

Other matrix algebra: [`tda_mcel()`](tda_mcel.md),
[`tda_mcent()`](tda_mcent.md), [`tda_mch()`](tda_mch.md),
[`tda_mcross()`](tda_mcross.md), [`tda_mdiag()`](tda_mdiag.md),
[`tda_mev()`](tda_mev.md), [`tda_mevs()`](tda_mevs.md),
[`tda_midf()`](tda_midf.md), [`tda_midf1()`](tda_midf1.md),
[`tda_midf2()`](tda_midf2.md), [`tda_midf3()`](tda_midf3.md),
[`tda_minvs()`](tda_minvs.md), [`tda_mkmet()`](tda_mkmet.md),
[`tda_mkp()`](tda_mkp.md), [`tda_mnc()`](tda_mnc.md),
[`tda_mnrow()`](tda_mnrow.md), [`tda_mpbl()`](tda_mpbl.md),
[`tda_mpfit()`](tda_mpfit.md), [`tda_mpinv()`](tda_mpinv.md),
[`tda_mpit()`](tda_mpit.md), [`tda_mple()`](tda_mple.md),
[`tda_mpz()`](tda_mpz.md), [`tda_mscal1()`](tda_mscal1.md),
[`tda_msqrtd()`](tda_msqrtd.md), [`tda_msvd()`](tda_msvd.md),
[`tda_mwvec()`](tda_mwvec.md)

## Examples

``` r
tda_mmul(matrix(1:4, 2), diag(2))
#>      [,1] [,2]
#> [1,]    1    3
#> [2,]    2    4
tda_mtransp(matrix(1:6, 2))
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    3    4
#> [3,]    5    6
tda_mchol(matrix(c(4, 2, 2, 3), 2))
#>      [,1]     [,2]
#> [1,]    2 0.000000
#> [2,]    1 1.414214
```
