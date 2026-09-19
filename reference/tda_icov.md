# Covariance and correlation of interval-valued variables

The interval of all covariance (or correlation) values reachable when
each observation varies in its interval – the same estimand as
[`tda_ivariance`](https://janmarvin.github.io/TDA/reference/tda_imean.md),
computed by TDA's `icov` and `icorr` (added to this build; not in
Rohwer's TDA). The formula names two interval pairs:
`~ iv(xlo, xhi) + iv(ylo, yhi)`.

## Usage

``` r
tda_icov(
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

tda_icorr(
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
```

## Arguments

- formula:

  `~ iv(xlo, xhi) + iv(ylo, yhi)`.

- data:

  a data frame with the four bound columns.

- max_boxes, max_iter, tol_width, tol_fd, tol_fe:

  search limits and tolerances, as in
  [`tda_ivariance`](https://janmarvin.github.io/TDA/reference/tda_imean.md).

- options, dir, ...:

  as everywhere.

## Value

An object with `$bounds`.

## Details

Naming the same pair twice is legal but means something specific: the
two copies vary independently over the intervals, so
`tda_icov(~ iv(a, b) + iv(a, b), d)` is not the interval variance (its
lower end can be negative) and the self-correlation is not fixed at 1.
The functions warn when they see this. For the variance of one interval
variable use
[`tda_ivar`](https://janmarvin.github.io/TDA/reference/tda_ivar.md).

Covariance certifies at realistic sizes. Correlation certifies small
problems; beyond that the search stops at its limits, returns the best
values found, and warns – raise `max_iter`/`max_boxes` or accept the
honest inner values.

## Examples

``` r
# three observations, both variables interval-valued; small enough
# that both searches certify instantly (verified against a dense
# grid in the package tests)
d <- data.frame(xl = c(1, 4, 6), xu = c(2, 5, 8),
                yl = c(3, 1, 6), yu = c(5, 2, 9))
tda_icov(~ iv(xl, xu) + iv(yl, yu), d,
         max_iter = 2e5, max_boxes = 1e5)
#> Call: .iv_pairs_stat("icov", formula, data, opts, dir, "Covariance")
#> 
#> Cases: 3 
#> Bounds: [-0.6666667, 7.555556]
tda_icorr(~ iv(xl, xu) + iv(yl, yu), d,
          max_iter = 2e5, max_boxes = 1e5)
#> Call: .iv_pairs_stat("icorr", formula, data, opts, dir, "Correlation")
#> 
#> Cases: 3 
#> Bounds: [-0.1428571, 0.8934051]
# \donttest{
# a realistic size: covariance certifies, correlation warns honestly
w <- tda_interval_wages()
tda_icov(~ iv(school_lo, school_hi) + iv(wage_lo, wage_hi), w,
         max_iter = 2e5, max_boxes = 1e5)
#> Call: .iv_pairs_stat("icov", formula, data, opts, dir, "Covariance")
#> 
#> Cases: 40 
#> Bounds: [2.441875, 13.42438]
# }
```
