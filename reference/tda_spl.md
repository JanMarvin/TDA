# Smoothing spline

A smoothing spline through a pair of variables, optionally weighted.
`smooth.spline` is the R analogue.

## Usage

``` r
tda_spl(
  x,
  y,
  weights = NULL,
  sig = 1,
  deg = NULL,
  max = 200,
  rx = NULL,
  options = list(),
  dir = tempfile("tda")
)
```

## Arguments

- x, y:

  the pair of variables.

- weights:

  optional per-point weights, the same idea as `weights` on the fitting
  functions.

- sig:

  smoothing factor; 0 interpolates.

- deg:

  degree of the spline.

- max:

  maximum number of knots.

- rx:

  range over which to evaluate, e.g. `seq(0, 10, 1)`.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

## Value

An object carrying a `table` of fitted values.

## See also

Other smoothing:
[`tda_integrate()`](https://janmarvin.github.io/TDA/reference/tda_integrate.md),
[`tda_interp()`](https://janmarvin.github.io/TDA/reference/tda_interp.md),
[`tda_isotonic()`](https://janmarvin.github.io/TDA/reference/tda_isotonic.md),
[`tda_mat()`](https://janmarvin.github.io/TDA/reference/tda_mat.md),
[`tda_sma()`](https://janmarvin.github.io/TDA/reference/tda_sma.md),
[`tda_smd()`](https://janmarvin.github.io/TDA/reference/tda_smd.md)

## Examples

``` r
set.seed(33)
x <- sort(runif(30, 0, 10))
y <- sin(x) + rnorm(30, sd = 0.15)
sp <- tda_spl(x, y, sig = 1)
head(sp$table)
#>   index         x         y    fitted          d1         d2
#> 1     1 0.1551696 0.1522848 0.2769707  0.88128812 -0.6805919
#> 2     2 0.4273416 0.3929474 0.4916244  0.69605006 -0.6805919
#> 3     3 1.1799116 0.9726723 0.8227202  0.18385699 -0.6805919
#> 4     4 1.3576507 0.9957205 0.8446484  0.06288921 -0.6805919
#> 5     5 2.2505121 0.6885517 0.6295154 -0.54478500 -0.6805919
#> 6     6 2.6048568 0.4449743 0.3937460 -0.78594920 -0.6805919
```
