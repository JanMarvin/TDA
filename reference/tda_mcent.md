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

Other matrix algebra: [`tda_mcel()`](tda_mcel.md),
[`tda_mch()`](tda_mch.md), [`tda_mcross()`](tda_mcross.md),
[`tda_mdiag()`](tda_mdiag.md), [`tda_mev()`](tda_mev.md),
[`tda_mevs()`](tda_mevs.md), [`tda_midf()`](tda_midf.md),
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
