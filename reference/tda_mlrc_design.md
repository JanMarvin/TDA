# Design matrix for MLRC models

The within-group pairwise design matrix TDA's `mldes` builds for MLRC
models: for objects with attribute rows `z` and a grouping in
consecutive runs, the row for the ordered pair (kk, k) holds the
products `z[k, i1] * z[kk, i2]` (i2-major) plus an indicator for k == kk
when both objects are in the same group, and zeros otherwise. Semantics
established against the source and verified cell-for-cell; see
`examples/coverage/mldesops.cf`.

## Usage

``` r
tda_mlrc_design(z, groups, dir = tempfile("tda"), ...)
```

## Arguments

- z:

  an n x q matrix, one row per object.

- groups:

  length-n vector; consecutive equal values form groups.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

The n^2 x (q^2 + 1) design matrix.

## See also

Other regression: [`TDA_FAMILIES`](tda_glm.md),
[`TDA_QRMODELS`](tda_qreg.md), [`tda_freg()`](tda_freg.md),
[`tda_gdf()`](tda_gdf.md), [`tda_l1reg()`](tda_l1reg.md),
[`tda_lsreg()`](tda_lsreg.md), [`tda_mreg()`](tda_mreg.md),
[`tda_nlreg()`](tda_nlreg.md), [`tda_npreg()`](tda_npreg.md),
[`tda_zreg()`](tda_zreg.md)

## Examples

``` r
z <- rbind(c(1, 2), c(3, 4), c(5, 6))
tda_mlrc_design(z, groups = c(1, 1, 2))
#>       [,1] [,2] [,3] [,4] [,5]
#>  [1,]    1    2    2    4    1
#>  [2,]    3    4    6    8    0
#>  [3,]    0    0    0    0    0
#>  [4,]    3    6    4    8    0
#>  [5,]    9   12   12   16    1
#>  [6,]    0    0    0    0    0
#>  [7,]    0    0    0    0    0
#>  [8,]    0    0    0    0    0
#>  [9,]   25   30   30   36    1
```
