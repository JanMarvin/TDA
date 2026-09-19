# Numerical-integration settings

TDA's `niset`: sets tolerance and subdivision options used by the `int`
operator within the same run; exposed for pipeline use.

## Usage

``` r
tda_niset(rel_error = 1e-04, method = 1, commands = character(), ...)
```

## Arguments

- rel_error:

  relative error target for the integrator.

- method:

  integration algorithm 1-5 (TDA's rhs of `niset`), default 1, QNG.

- commands:

  further command lines to run under the setting.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

the printed output, invisibly.

## Examples

``` r
out <- tda_niset(rel_error = 1e-6,
                 commands = "int(ab=0,1, fmt=12.8) = sin(x);")
grep("Approximation", out, value = TRUE)   # 1 - cos(1)
#> [1] "Approximation:   0.45969769 "
```
