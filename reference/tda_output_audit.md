# Account for every line of a run's printed output

Classifies each line of `output` as blank, boilerplate, an echo of the
inputs, an optimiser-settings line, content the fitted object stores,
content also delivered through the direct-export channel, or
*unaccounted*. The last category is the point: it lists what would be
lost if the printed output were ever discarded, and the package's tests
require it to be empty for a corpus of fits.

## Usage

``` r
tda_output_audit(x)
```

## Arguments

- x:

  a `tda_result` (from [`tda_run`](tda_run.md)) or any fitted object
  carrying one as `$run`.

## Value

a data frame with columns `line` and `class`, one row per output line,
with a `summary` attribute counting each class.

## Examples

``` r
d <- data.frame(x = 1:8, y = (1:8) * 2 + rnorm(8))
f <- tda_lsreg(y ~ x, d)
a <- tda_output_audit(f)
attr(a, "summary")
#> cls
#>       blank boilerplate        echo    settings      stored 
#>          15           9           7           3          12 
subset(a, class == "UNACCOUNTED")
#> [1] line  class
#> <0 rows> (or 0-length row.names)
```
