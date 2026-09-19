# Variance of an interval-valued variable (TDA's ivar)

Identical to [`tda_ivariance`](tda_imean.md): TDA's `ivar`
branch-and-bound with tunable limits. This is the command-style name;
[`tda_ivar1`](tda_ivar1.md) wraps the newer `ivar1` algorithm, and all
three compute the same variance bounds.

## Usage

``` r
tda_ivar(
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

- formula, data, max_boxes, max_iter, tol_width, tol_fd, tol_fe,
  options, dir, ...:

  exactly as in [`tda_ivariance`](tda_imean.md).

## Value

See [`tda_ivariance`](tda_imean.md).

## Examples

``` r
d <- data.frame(lo = c(1, 3, 5, 9), hi = c(2, 4, 8, 12))
tda_ivar(~ iv(lo, hi), d)
#> Call: tda_ivariance(formula, data, max_boxes = max_boxes, max_iter = max_iter, 
#>     tol_width = tol_width, tol_fd = tol_fd, tol_fe = tol_fe, 
#>     options = options, dir = dir, ...)
#> 
#> Cases: 4 
#> Bounds: [6.5, 18.5]
```
