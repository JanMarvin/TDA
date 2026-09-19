# Distribution function of an interval-valued variable

`tda_idf`: the distribution and density function of a continuous
interval-valued variable, as bounds at each partition point – an
interval variable does not pin down a single distribution, only a range
one is consistent with. `tda_iddf` is the same idea for a discrete
interval-valued variable (a finite set of possible values, not a
continuum), and, unlike `idf`, does take `opt=` for a self-consistent
distribution; `idf` takes only `fmt=`.

## Usage

``` r
tda_idf(
  formula,
  data,
  self_consistent = FALSE,
  control = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)

tda_iddf(
  formula,
  data,
  self_consistent = FALSE,
  control = NULL,
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

- self_consistent:

  also compute the self-consistent distribution (`opt=2`), not just the
  min/max/mean bounds (`opt=1`, default) – `opt=` of both commands. For
  `tda_idf` it comes back as `$self_consistent` (partition point and
  distribution function) with `$converged`; the fixed point satisfies
  the self-consistency equation to the tolerance (checked in the tests
  against the equation iterated in R).

- control:

  convergence settings (`mxit=`/`tolf` via
  [`tda_control`](https://janmarvin.github.io/TDA/reference/tda_control.md)),
  defaults 50 and 0.001.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

An object carrying a `table` of partition points and the distribution
function's bounds and mean at each.

## See also

Other interval-valued data:
[`tda_bounds()`](https://janmarvin.github.io/TDA/reference/tda_bounds.md),
[`tda_ilsreg()`](https://janmarvin.github.io/TDA/reference/tda_ilsreg.md),
[`tda_imean()`](https://janmarvin.github.io/TDA/reference/tda_imean.md),
[`tda_imreg()`](https://janmarvin.github.io/TDA/reference/tda_imreg.md),
[`tda_inpreg()`](https://janmarvin.github.io/TDA/reference/tda_inpreg.md),
[`tda_ivar1()`](https://janmarvin.github.io/TDA/reference/tda_ivar1.md),
[`tda_ivreg()`](https://janmarvin.github.io/TDA/reference/tda_ivreg.md),
[`tda_sddf()`](https://janmarvin.github.io/TDA/reference/tda_sddf.md)

## Examples

``` r
set.seed(1)
n <- 30
d <- data.frame(lo = round(rnorm(n, 5, 2), 1))
d$hi <- d$lo + round(runif(n, 0.5, 2), 1)
fit <- tda_idf(~ iv(lo, hi), d)
head(fit$table)
#>   index partition      lower      upper    mean_df
#> 1     1       0.6 0.00000000 0.03333333 0.00000000
#> 2     2       1.0 0.00000000 0.06666667 0.01333333
#> 3     3       1.6 0.03333333 0.06666667 0.05333333
#> 4     4       2.0 0.06666667 0.06666667 0.06666667
#> 5     5       2.1 0.06666667 0.10000000 0.06666667
#> 6     6       2.8 0.10000000 0.10000000 0.10000000

# iddf: the discrete counterpart, with the self-consistent option idf
# itself does not have
fit2 <- tda_iddf(~ iv(lo, hi), d, self_consistent = TRUE)
head(fit2$table)
#>   index      lower      upper    mean_df self_consistent_df
#> 1     0 0.00000000 0.03333333 0.01666667       3.483206e-05
#> 2     1 0.00000000 0.06666667 0.03333333       4.994775e-02
#> 3     2 0.03333333 0.06666667 0.05000000       5.001742e-02
#> 4     3 0.00000000 0.13333333 0.06111111       1.392193e-05
#> 5     4 0.03333333 0.30000000 0.14444444       2.290678e-01
#> 6     5 0.00000000 0.36666667 0.15000000       3.318755e-03
```
