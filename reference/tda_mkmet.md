# Kemeny distance between rankings

Pairwise Kemeny (Kendall-tau-with-ties) distance between every pair of
rows of `A`, each row a ranking (or score) over the same set of items –
TDA's `mkmet`. For two rankings, each pair of items contributes 0 if
both rankings order (or tie) the pair the same way, 1 if one ranking
ties the pair and the other doesn't, or 2 if the rankings strictly
disagree on the pair's order. Confirmed by hand against two 3-item
rankings, one pair with a tie.

## Usage

``` r
tda_mkmet(A, dir = tempfile("tda"))
```

## Arguments

- A:

  a numeric matrix, one ranking per row.

- dir:

  working directory.

## Value

A symmetric numeric matrix with zero diagonal, the Kemeny distance
between every pair of rows.

## See also

Other matrix algebra: [`tda_mcel()`](tda_mcel.md),
[`tda_mcent()`](tda_mcent.md), [`tda_mch()`](tda_mch.md),
[`tda_mcross()`](tda_mcross.md), [`tda_mdiag()`](tda_mdiag.md),
[`tda_mev()`](tda_mev.md), [`tda_mevs()`](tda_mevs.md),
[`tda_midf()`](tda_midf.md), [`tda_midf1()`](tda_midf1.md),
[`tda_midf2()`](tda_midf2.md), [`tda_midf3()`](tda_midf3.md),
[`tda_minvs()`](tda_minvs.md), [`tda_mkp()`](tda_mkp.md),
[`tda_mmul()`](tda_mmul.md), [`tda_mnc()`](tda_mnc.md),
[`tda_mnrow()`](tda_mnrow.md), [`tda_mpbl()`](tda_mpbl.md),
[`tda_mpfit()`](tda_mpfit.md), [`tda_mpinv()`](tda_mpinv.md),
[`tda_mpit()`](tda_mpit.md), [`tda_mple()`](tda_mple.md),
[`tda_mpz()`](tda_mpz.md), [`tda_mscal1()`](tda_mscal1.md),
[`tda_msqrtd()`](tda_msqrtd.md), [`tda_msvd()`](tda_msvd.md),
[`tda_mwvec()`](tda_mwvec.md)

## Examples

``` r
tda_mkmet(rbind(c(1, 2, 3), c(3, 2, 1), c(1, 2, 2)))
#>      [,1] [,2] [,3]
#> [1,]    0    6    1
#> [2,]    6    0    5
#> [3,]    1    5    0
```
