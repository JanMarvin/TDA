# Statistics of an interval-valued variable

The interval-valued counterparts of the ordinary summaries: a mean, a
variance and a Gini coefficient computed from `(lower, upper)` pairs
rather than from points. Each is reported as bounds, because an interval
variable does not have a single value for these either.

## Usage

``` r
tda_imean(formula, data, options = list(), dir = tempfile("tda"), ...)

tda_ivariance(
  formula,
  data,
  max_boxes = NULL,
  max_iter = NULL,
  tol_width = NULL,
  tol_fd = NULL,
  tol_fe = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)

tda_igini(
  formula,
  data,
  x = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- formula:

  a one-sided formula naming one interval-valued variable,
  `~ iv(lo, hi)`.

- data:

  a data frame.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

- max_boxes, max_iter, tol_width, tol_fd, tol_fe:

  for `tda_ivariance` only, the branch-and-bound search's limits and
  tolerances (`nbox=`/`mxit=`/`tolbw=`/ `tolfd=`/`tolfe=`). `tda_imean`
  and `tda_igini` take no such options.

- x:

  for `tda_igini`, the points at which the Lorenz curve is evaluated.

## Value

An object carrying the bounds in `$bounds`, and the run.

## See also

Other interval-valued data:
[`tda_bounds()`](https://janmarvin.github.io/TDA/reference/tda_bounds.md),
[`tda_idf()`](https://janmarvin.github.io/TDA/reference/tda_idf.md),
[`tda_ilsreg()`](https://janmarvin.github.io/TDA/reference/tda_ilsreg.md),
[`tda_imreg()`](https://janmarvin.github.io/TDA/reference/tda_imreg.md),
[`tda_inpreg()`](https://janmarvin.github.io/TDA/reference/tda_inpreg.md),
[`tda_ivar1()`](https://janmarvin.github.io/TDA/reference/tda_ivar1.md),
[`tda_ivreg()`](https://janmarvin.github.io/TDA/reference/tda_ivreg.md),
[`tda_sddf()`](https://janmarvin.github.io/TDA/reference/tda_sddf.md)

## Examples

``` r
d <- data.frame(lo = c(1, 3, 5, 9), hi = c(2, 4, 8, 12))
tda_imean(~ iv(lo, hi), d)
#> Call: tda_imean(~iv(lo, hi), d)
#> 
#> Cases: 4 
#> Bounds: [4.5, 6.5]
tda_ivariance(~ iv(lo, hi), d)
#> Call: tda_ivariance(~iv(lo, hi), d)
#> 
#> Cases: 4 
#> Bounds: [6.5, 18.5]
tda_igini(~ iv(lo, hi), d)
#> Call: tda_igini(formula = ~iv(lo, hi), data = d)
#> 
#> Cases: 4 
#> 
#>    x lower upper
#>  1.0     0     0
#>  1.1     0     0
#>  1.2     0     0
#>  1.3     0     0
#>  1.4     0     0
#>  1.5     0     0
#>  1.6     0     0
#>  1.7     0     0
#>  1.8     0     0
#>  1.9     0     0
#>  2.0     1     1
```
