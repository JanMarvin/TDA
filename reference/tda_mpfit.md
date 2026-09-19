# Iterative proportional fitting

Adjusts a table `A` so its row and column sums match prescribed targets
while preserving cross-product ratios (the RAS/IPF algorithm) – TDA's
`mpfit`, confirmed against a from-scratch RAS implementation on a small
table. `A` must be non-negative with strictly positive row and column
sums, as must the target `row_sums`/`col_sums`, and their totals should
agree (TDA warns but does not refuse otherwise).

## Usage

``` r
tda_mpfit(
  A,
  row_sums,
  col_sums,
  max_iter = 100,
  eps = 1e-06,
  dir = tempfile("tda")
)
```

## Arguments

- A:

  a non-negative numeric matrix (n x m).

- row_sums:

  target row sums, length n.

- col_sums:

  target column sums, length m.

- max_iter:

  maximum number of iterations.

- eps:

  convergence tolerance on the maximum deviation between the fitted and
  target sums.

- dir:

  working directory.

## Value

A list: `B`, the fitted table; `iterations`, how many were used;
`accuracy`, the final deviation reached.

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
[`tda_mpbl()`](tda_mpbl.md), [`tda_mpinv()`](tda_mpinv.md),
[`tda_mpit()`](tda_mpit.md), [`tda_mple()`](tda_mple.md),
[`tda_mpz()`](tda_mpz.md), [`tda_mscal1()`](tda_mscal1.md),
[`tda_msqrtd()`](tda_msqrtd.md), [`tda_msvd()`](tda_msvd.md),
[`tda_mwvec()`](tda_mwvec.md)

## Examples

``` r
tda_mpfit(matrix(c(10, 20, 30, 40), 2, 2, byrow = TRUE),
         row_sums = c(45, 55), col_sums = c(40, 60))
#> $B
#>          [,1]     [,2]
#> [1,] 15.61072 29.38928
#> [2,] 24.38928 30.61072
#> 
#> $iterations
#> [1] 5
#> 
#> $accuracy
#> [1] 1.21766e-08
#> 
```
