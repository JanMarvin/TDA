# Centering, standardizing, and double-centering

`tda_mcent` subtracts each column's mean – TDA's `mcent`, matching
`sweep(A, 2, colMeans(A))`. `tda_mstand` additionally divides each
column by its *population* standard deviation (denominator `n`, not base
R's [`sd()`](https://rdrr.io/r/stats/sd.html)'s `n-1`) – TDA's `mstand`;
confirmed by hand, since this is exactly the kind of denominator choice
that silently disagrees with the R-native way of doing the same thing if
assumed rather than checked. `tda_mdcent` double-centers a symmetric
matrix – TDA's `mdcent`, the classical (Torgerson) transformation
multidimensional scaling uses to turn a matrix of *squared* distances
into one whose eigendecomposition gives point coordinates.

## Usage

``` r
tda_mcent(A, dir = tempfile("tda"))

tda_mstand(A, dir = tempfile("tda"))

tda_mdcent(D, dir = tempfile("tda"))
```

## Arguments

- A, D:

  a numeric matrix (`D`, for `tda_mdcent`, should be symmetric –
  typically a distance matrix).

- dir:

  working directory.

## Value

A numeric matrix.

## See also

Other matrix algebra:
[`tda_mcel()`](https://janmarvin.github.io/TDA/reference/tda_mcel.md),
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
[`tda_mmul()`](https://janmarvin.github.io/TDA/reference/tda_mmul.md),
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
A <- matrix(c(1, 2, 6, 10, 30, 20), 3, 2)
tda_mcent(A)
#>      [,1] [,2]
#> [1,]   -2  -10
#> [2,]   -1   10
#> [3,]    3    0
tda_mstand(A)
#>            [,1]      [,2]
#> [1,] -0.9258201 -1.224745
#> [2,] -0.4629100  1.224745
#> [3,]  1.3887301  0.000000
D <- matrix(c(0, 3, 4, 3, 0, 5, 4, 5, 0), 3, 3)
tda_mdcent(D)
#>            [,1]       [,2]      [,3]
#> [1,]  2.7777778 -0.2222222 -2.555556
#> [2,] -0.2222222  5.7777778 -5.555556
#> [3,] -2.5555556 -5.5555556  8.111111
```
