# General distribution functions for censored data

TDA's `gdf`: the marginal (or, with `method`, joint) distribution or
survivor function for possibly right-censored data, via a
Kaplan-Meier-style calculation – the same underlying censoring machinery
as
[`tda_lsreg1`](https://janmarvin.github.io/TDA/reference/tda_lsreg.md),
but estimating the distribution function directly rather than a
conditional expectation.

## Usage

``` r
tda_gdf(
  formula,
  data,
  censor = NULL,
  id = NULL,
  dimension = NULL,
  what = c("distribution", "survivor", "expected_values"),
  method = c("marginal", "joint1", "joint2"),
  offset = NULL,
  n_boxes = NULL,
  delta = NULL,
  control = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- formula:

  a one-sided formula naming the (possibly censored) variable, e.g.
  `~ Y`, or a list of such formulas – one per dimension – for the
  multi-dimensional, `grp=`-based form (see
  [`tda_lsreg1`](https://janmarvin.github.io/TDA/reference/tda_lsreg.md)'s
  Details for how `id` groups rows into units).

- data:

  a data frame.

- censor:

  the censoring indicator, as a name or a vector, or `NULL` (the
  default) when every observation is exact: `TRUE` or 1 marks a case
  whose value is censored (the opposite of
  [`survival::Surv`](https://rdrr.io/pkg/survival/man/Surv.html)'s event
  indicator). TDA's `cen=` is 1 for an exact observation and is
  translated here.

- id:

  for multi-dimensional data, the name of the column identifying which
  unit each row belongs to (TDA's `grp=`).

- dimension:

  with `id`, the name of the column saying which margin (dimension) each
  of a unit's observations belongs to – TDA's level variable in
  `grp=ID,L1`. Positive integers, not necessarily contiguous; TDA
  requires that no two observations of one unit share a value, so this
  is mandatory whenever an `id` value repeats. The data are in long
  format: one row per (unit, margin), and the joint methods estimate the
  joint distribution over the margins.

- what:

  the distribution function (default), the survivor function, or the
  expected value of each case's observation under the fitted
  distribution – `gdf`'s `prn=`. The third option returns a
  differently-shaped table: one row per case (not per distinct value),
  carrying the observed value, its exact/censored status (1 = exact,
  TDA's cen value), and its expected value.

- method:

  `"marginal"` (default) or, for multi-dimensional data (with `id`),
  `"joint1"`/`"joint2"` – `gdf`'s own `opt=`, the same three values
  `tda_lsreg1`'s `method` takes.

- offset, n_boxes:

  for `method = "marginal"` or `"joint1"` only, the domain offset and
  the number of grid boxes (`gdf`'s `sc=`/`n=`).

- delta:

  for `method = "joint2"` only, the delta grid spacing (`gdf`'s `d=`).

- control:

  for `method = "marginal"` or `"joint1"` only, convergence settings
  (only
  [`tda_control`](https://janmarvin.github.io/TDA/reference/tda_control.md)'s
  `maxit`/`tolf` fields apply, as `gdf`'s `mxit=`/`tolf=`).

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

An object carrying a `table` of values and the distribution (or
survivor) function at each.

## See also

Other regression:
[`TDA_FAMILIES`](https://janmarvin.github.io/TDA/reference/tda_glm.md),
[`TDA_QRMODELS`](https://janmarvin.github.io/TDA/reference/tda_qreg.md),
[`tda_freg()`](https://janmarvin.github.io/TDA/reference/tda_freg.md),
[`tda_l1reg()`](https://janmarvin.github.io/TDA/reference/tda_l1reg.md),
[`tda_lsreg()`](https://janmarvin.github.io/TDA/reference/tda_lsreg.md),
[`tda_mlrc_design()`](https://janmarvin.github.io/TDA/reference/tda_mlrc_design.md),
[`tda_mreg()`](https://janmarvin.github.io/TDA/reference/tda_mreg.md),
[`tda_nlreg()`](https://janmarvin.github.io/TDA/reference/tda_nlreg.md),
[`tda_npreg()`](https://janmarvin.github.io/TDA/reference/tda_npreg.md),
[`tda_zreg()`](https://janmarvin.github.io/TDA/reference/tda_zreg.md)

## Examples

``` r
set.seed(1)
n <- 20
d <- data.frame(y = round(rexp(n, 0.2), 1), cen = rbinom(n, 1, 0.3))
fit <- tda_gdf(~ y, d, censor = "cen")
fit$table
#>    dimension value distribution_function
#> 1          1   0.7             0.1000000
#> 2          1   1.7             0.1529412
#> 3          1   2.2             0.2058824
#> 4          1   2.7             0.2588235
#> 5          1   2.9             0.3117647
#> 6          1   3.3             0.3647059
#> 7          1   3.8             0.4176471
#> 8          1   5.2             0.4823529
#> 9          1   6.1             0.5686275
#> 10         1   6.2             0.6549020
#> 11         1   7.0             0.7411765
#> 12         1   9.4             0.8274510
#> 13         1  22.1             1.0000000

# expected values: one row per case, not per distinct value
fit2 <- tda_gdf(~ y, d, censor = "cen", what = "expected_values")
head(fit2$table)
#>   dimension    y exact  expected
#> 1         1  3.8     1  3.800000
#> 2         1  5.9     0 12.150000
#> 3         1  0.7     1  0.700000
#> 4         1  0.7     0  8.338562
#> 5         1  2.2     1  2.200000
#> 6         1 14.5     0 22.100000
```
