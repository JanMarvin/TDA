# Procrustes rotation

TDA's `mproc`: rotates configuration Y to best match X. Verified in the
suite against the closed-form solution from the singular value
decomposition.

## Usage

``` r
tda_mproc(x, y, ...)
```

## Arguments

- x, y:

  n-by-m matrices.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

the rotated configuration as a matrix.

## Examples

``` r
x <- cbind(c(0, 1, 0), c(0, 0, 1))
th <- 0.7; R <- rbind(c(cos(th), -sin(th)), c(sin(th), cos(th)))
round(tda_mproc(x, x %*% R), 6)
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    1    0
#> [3,]    0    1
```
