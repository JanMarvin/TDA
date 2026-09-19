# Bounds from an interval regression

An interval regression identifies its coefficient only up to an
interval, so [`coef()`](https://rdrr.io/r/stats/coef.html) would be
misleading. This returns the bounds TDA reports.

## Usage

``` r
tda_bounds(x, what = c("beta", "alpha", "mean_y", "mean_x", "var_x"))
```

## Arguments

- x:

  a fit from [`tda_ivreg`](tda_ivreg.md) and friends.

- what:

  which bounds: the slope, the intercept (`tda_ilsreg` only), or the
  mean or variance summaries.

## Value

A named vector with `lower` and `upper`, or `NULL`.

## See also

Other interval-valued data: [`tda_idf()`](tda_idf.md),
[`tda_ilsreg()`](tda_ilsreg.md), [`tda_imean()`](tda_imean.md),
[`tda_imreg()`](tda_imreg.md), [`tda_inpreg()`](tda_inpreg.md),
[`tda_ivar1()`](tda_ivar1.md), [`tda_ivreg()`](tda_ivreg.md),
[`tda_sddf()`](tda_sddf.md)

## Examples

``` r
set.seed(37)
n <- 40
xlo <- round(rnorm(n), 2); xhi <- xlo + round(runif(n, 0.5, 1.5), 2)
ylo <- round(1 + 0.6 * xlo + rnorm(n, sd = 0.3), 2)
yhi <- ylo + round(runif(n, 0.5, 1.5), 2)
d <- data.frame(xlo = xlo, xhi = xhi, ylo = ylo, yhi = yhi)
ivf <- tda_ivreg(iv(ylo, yhi) ~ iv(xlo, xhi), d)
tda_bounds(ivf, "beta")
#>     lower     upper 
#> -1.168092  4.142178 
tda_bounds(ivf, "mean_x")
#>    lower    upper 
#> -0.07550  0.87275 

# ilsreg: an interval response, a point regressor -- beta and alpha (the
# intercept) both come back bounded, not just the slope
set.seed(2)
d2 <- data.frame(x = round(rnorm(30), 2))
d2$ylo <- 2 + 0.5 * d2$x + rnorm(30, sd = 0.2)
d2$yhi <- d2$ylo + round(runif(30, 0.3, 1), 2)
lsf <- tda_ilsreg(iv(ylo, yhi) ~ x, d2)
tda_bounds(lsf, "beta")
#>     lower     upper 
#> 0.2378364 0.6916871 
tda_bounds(lsf, "alpha")
#>    lower    upper 
#> 2.270285 2.367430 
#> attr(,"at_beta")
#> at_beta_min at_beta_max 
#>    2.367430    2.270285 
```
