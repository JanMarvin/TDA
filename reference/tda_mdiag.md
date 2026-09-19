# Diagonal construction and extraction

`tda_mdiag` builds an n x n diagonal matrix from a length-n vector –
TDA's `mdiag`, matching base R's `diag(x)` when `x` has length \> 1.
`tda_mdiagd` extracts the diagonal of a matrix as a plain vector, using
`min(nrow(A), ncol(A))` entries for a non-square `A` – TDA's `mdiagd`,
matching base R's `diag(A)`.

## Usage

``` r
tda_mdiag(x, dir = tempfile("tda"))

tda_mdiagd(A, dir = tempfile("tda"))
```

## Arguments

- x:

  a numeric vector.

- dir:

  working directory.

- A:

  a numeric matrix.

## Value

`tda_mdiag` returns a matrix; `tda_mdiagd` returns a plain numeric
vector.

## See also

Other matrix algebra: [`tda_mcel()`](tda_mcel.md),
[`tda_mcent()`](tda_mcent.md), [`tda_mch()`](tda_mch.md),
[`tda_mcross()`](tda_mcross.md), [`tda_mev()`](tda_mev.md),
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
tda_mdiag(c(5, 7, 9))
#>      [,1] [,2] [,3]
#> [1,]    5    0    0
#> [2,]    0    7    0
#> [3,]    0    0    9
tda_mdiagd(matrix(1:6, 2, 3))
#> [1] 1 4
```
