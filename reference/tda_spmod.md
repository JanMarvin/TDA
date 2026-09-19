# Stationary solution of a Leslie matrix

TDA's `spmod`: iterates a population vector under the Leslie matrix
built from age-specific fertility (first row) and survivor rates
(subdiagonal) until the age distribution is stationary. Verified in the
suite against R's `eigen`: the growth factor equals the dominant
eigenvalue and the stationary vector the normalized dominant
eigenvector.

## Usage

``` r
tda_spmod(
  fertility,
  survival,
  start = rep(1, length(fertility)),
  tolf = 1e-08,
  mxit = 500,
  ...
)
```

## Arguments

- fertility, survival, start:

  numeric vectors of equal length: age-specific fertility, survivor
  rates (last entry unused), and the starting population.

- tolf:

  convergence tolerance.

- mxit:

  iteration cap.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

list with `growth` (the dominant eigenvalue) and `stationary` (the
stable age distribution, summing to 1).

## Examples

``` r
r <- tda_spmod(c(0, .4, .3, .1), c(.9, .8, .7, 0), rep(1, 4))
r
#> <tda_leslie>
#> growth: 0.8320931 
#> stationary: 0.2502305 0.2706518 0.2602130 0.2189047 
L <- rbind(c(0, .4, .3, .1), cbind(diag(c(.9, .8, .7)), 0))
stopifnot(abs(r$growth - Re(eigen(L)$values[1])) < 1e-5)
```
