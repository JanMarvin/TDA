# Iterative proportional fitting

Adjusts a table `A` so its row and column sums match prescribed targets
while preserving cross-product ratios (the RAS/IPF algorithm) – TDA's
`mpfit`. `A` must be non-negative with strictly positive row and column
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
[`tda_mmul()`](https://janmarvin.github.io/TDA/reference/tda_mmul.md),
[`tda_mnc()`](https://janmarvin.github.io/TDA/reference/tda_mnc.md),
[`tda_mnrow()`](https://janmarvin.github.io/TDA/reference/tda_mnrow.md),
[`tda_mpbl()`](https://janmarvin.github.io/TDA/reference/tda_mpbl.md),
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
