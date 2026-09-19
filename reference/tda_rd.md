# Random values, as TDA's `rd()`/`rdn()`

The R-side equivalent of what `nvar()`'s `X = rd(0, 3)` or
`Y = sin(X) + rd` would generate inside TDA: `tda_rd(n, a, b)` is
uniform on `[a, b]` (`rd`'s two-argument form,
`a + random1() * (b - a)`), and `tda_rd(n)` alone – `rd` with no
arguments – is uniform on `[0, 1]`, not the wider, arbitrary range a
hand-rolled substitute might reach for. `tda_rdn` is the same idea for
`rdn`, a standard normal.

## Usage

``` r
tda_rd(n, a = 0, b = 1)

tda_rdn(n, mean = 0, sd = 1)
```

## Arguments

- n:

  how many values to draw.

- a, b:

  the range, for `tda_rd`; TDA's default (no arguments) is `[0, 1]`.

- mean, sd:

  the normal distribution's parameters, for `tda_rdn`; TDA's `rdn` is
  always standard (0, 1).

## Value

A numeric vector of length `n`.

## Details

These use R's random number generator, not [`tda_runif`](tda_rng.md)'s
reproduction of TDA's exact stream – for building a data frame to pass
to TDA (rather than reproducing an existing TDA example's numbers bit
for bit), R's generator is what
[`set.seed()`](https://rdrr.io/r/base/Random.html) already controls.

## Examples

``` r
set.seed(1)
x <- tda_rd(200, 0, 3)      # nvar()'s X = rd(0, 3)
y <- sin(x) + tda_rd(200)   # nvar()'s Y = sin(X) + rd
round(c(mean = mean(x), sd = sd(x)), 2)   # near 0 and 3
#> mean   sd 
#> 1.51 0.90 
```
