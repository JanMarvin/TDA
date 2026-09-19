# Scale a matrix by its total sum

Divides every element of `A` by the sum of all its elements – TDA's
`mscal1`, matching `A / sum(A)`.

## Usage

``` r
tda_mscal1(A, dir = tempfile("tda"))
```

## Arguments

- A:

  a numeric matrix whose elements do not sum to zero.

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
[`tda_minvs()`](tda_minvs.md), [`tda_mkmet()`](tda_mkmet.md),
[`tda_mkp()`](tda_mkp.md), [`tda_mmul()`](tda_mmul.md),
[`tda_mnc()`](tda_mnc.md), [`tda_mnrow()`](tda_mnrow.md),
[`tda_mpbl()`](tda_mpbl.md), [`tda_mpfit()`](tda_mpfit.md),
[`tda_mpinv()`](tda_mpinv.md), [`tda_mpit()`](tda_mpit.md),
[`tda_mple()`](tda_mple.md), [`tda_mpz()`](tda_mpz.md),
[`tda_msqrtd()`](tda_msqrtd.md), [`tda_msvd()`](tda_msvd.md),
[`tda_mwvec()`](tda_mwvec.md)

## Examples

``` r
tda_mscal1(matrix(c(2, 4, 6, 8), 2))
#>      [,1] [,2]
#> [1,]  0.1  0.3
#> [2,]  0.2  0.4
```
