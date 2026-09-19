# Conjoint analysis

TDA's `conj` (alternating least squares) on a preference variable and
factor columns.

## Usage

``` r
tda_conj(y, x, algorithm = c("monotone", "lp"), ...)
```

## Arguments

- y:

  preference values.

- x:

  factor matrix.

- algorithm:

  "monotone" regression or "lp" (linear programming).

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

the printed output.

## Examples

``` r
x <- expand.grid(a = 1:2, b = 1:2)
r <- tda_conj(y = c(4, 3, 2, 1), x = x)
cat(tail(tda_payload(r), 8), sep = "\n")   # the part-worths
#> Final parameter change: 1
#> Final norm of least squares residuals: 2.22045e-16
#> Rank of design matrix: 3
#> Idx  Dimension  Variable    Category  Coefficient
#>   0       -      Const           -        4.0000 
#>   1       1      X1              2       -1.0000 
#>   2       2      X2              2       -2.0000 
#> Type of equivalence: 1
```
