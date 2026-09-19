# Moving average

A weighted moving average over a set of cases, TDA's `sma`.

## Usage

``` r
tda_sma(
  x,
  width = NULL,
  weights = NULL,
  options = list(),
  dir = tempfile("tda")
)
```

## Arguments

- x:

  the variable to average.

- width:

  an odd window size for an equal-weight moving average, e.g. `3` for
  the average of each case with its immediate neighbour on each side.
  Give this or `weights`, not both.

- weights:

  TDA's half-kernel directly (see Details), for an unequal-weight window
  – `c(0.5, 0.25)` for a 3-point window weighted `0.25, 0.5, 0.25`.

- options:

  a named list of further TDA options, passed through: `opt = 1`
  (default) or `2` for how the ends are handled, `r=` a repeat factor.

- dir:

  working directory.

## Value

An object carrying a `table`: the case index, the original value, and
the smoothed value.

## Details

`gss=` is TDA's weight specification, and it takes only half a symmetric
kernel: the centre weight followed by the weights moving out to *one*
side, mirrored automatically to build the full window – `c(1, 1)`
becomes the three weights `1, 1, 1`; `c(1, 1, 1)` becomes five,
`1, 1, 1, 1, 1`. There is no `tp=` at all: a plain equal-weight window
over consecutive cases is what `sma` does, not a smooth over a separate
time axis. Weights that sum to 1 give an ordinary average, matching
[`zoo::rollmean()`](https://rdrr.io/pkg/zoo/man/rollmean.html) for the
equal-weight case aside from TDA's own end-value rule where `rollmean`
would give `NA`.

## Examples

``` r
x <- c(1, 2, 3, 10, 5, 6, 7, 8, 9, 10)
tda_sma(x, width = 3)$table       # equal-weight 3-point average
#>    sel index value smoothed
#> 1    0     1     1        1
#> 2    0     2     2        2
#> 3    0     3     3        5
#> 4    0     4    10        6
#> 5    0     5     5        7
#> 6    0     6     6        6
#> 7    0     7     7        7
#> 8    0     8     8        8
#> 9    0     9     9        9
#> 10   0    10    10       10
tda_sma(x, weights = c(0.5, 0.25))$table  # weighted 0.25, 0.5, 0.25
#>    sel index value smoothed
#> 1    0     1     1      1.0
#> 2    0     2     2      2.0
#> 3    0     3     3      4.5
#> 4    0     4    10      7.0
#> 5    0     5     5      6.5
#> 6    0     6     6      6.0
#> 7    0     7     7      7.0
#> 8    0     8     8      8.0
#> 9    0     9     9      9.0
#> 10   0    10    10     10.0
```
