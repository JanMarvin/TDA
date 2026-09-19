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

Other matrix algebra:
[`tda_mcel()`](https://janmarvin.github.io/TDA/reference/tda_mcel.md),
[`tda_mcent()`](https://janmarvin.github.io/TDA/reference/tda_mcent.md),
[`tda_mch()`](https://janmarvin.github.io/TDA/reference/tda_mch.md),
[`tda_mcross()`](https://janmarvin.github.io/TDA/reference/tda_mcross.md),
[`tda_mdiag()`](https://janmarvin.github.io/TDA/reference/tda_mdiag.md),
[`tda_mev()`](https://janmarvin.github.io/TDA/reference/tda_mev.md),
[`tda_mevs()`](https://janmarvin.github.io/TDA/reference/tda_mevs.md),
[`tda_midf()`](https://janmarvin.github.io/TDA/reference/tda_midf.md),
[`tda_midf1()`](https://janmarvin.github.io/TDA/reference/tda_midf1.md),
[`tda_midf2()`](https://janmarvin.github.io/TDA/reference/tda_midf2.md),
[`tda_midf3()`](https://janmarvin.github.io/TDA/reference/tda_midf3.md),
[`tda_minvs()`](https://janmarvin.github.io/TDA/reference/tda_minvs.md),
[`tda_mkmet()`](https://janmarvin.github.io/TDA/reference/tda_mkmet.md),
[`tda_mkp()`](https://janmarvin.github.io/TDA/reference/tda_mkp.md),
[`tda_mnc()`](https://janmarvin.github.io/TDA/reference/tda_mnc.md),
[`tda_mnrow()`](https://janmarvin.github.io/TDA/reference/tda_mnrow.md),
[`tda_mpbl()`](https://janmarvin.github.io/TDA/reference/tda_mpbl.md),
[`tda_mpfit()`](https://janmarvin.github.io/TDA/reference/tda_mpfit.md),
[`tda_mpinv()`](https://janmarvin.github.io/TDA/reference/tda_mpinv.md),
[`tda_mpit()`](https://janmarvin.github.io/TDA/reference/tda_mpit.md),
[`tda_mple()`](https://janmarvin.github.io/TDA/reference/tda_mple.md),
[`tda_mpz()`](https://janmarvin.github.io/TDA/reference/tda_mpz.md),
[`tda_mscal1()`](https://janmarvin.github.io/TDA/reference/tda_mscal1.md),
[`tda_msqrtd()`](https://janmarvin.github.io/TDA/reference/tda_msqrtd.md),
[`tda_msvd()`](https://janmarvin.github.io/TDA/reference/tda_msvd.md),
[`tda_mwvec()`](https://janmarvin.github.io/TDA/reference/tda_mwvec.md)

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
