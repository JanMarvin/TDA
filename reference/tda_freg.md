# Non-linear regression with a user-defined function

`freg` fits a regression function you write yourself, in the same way
[`tda_fml`](https://janmarvin.github.io/TDA/reference/tda_fml.md) takes
a likelihood. The expression is TDA's language, not R's.

## Usage

``` r
tda_freg(
  definitions,
  data,
  start = NULL,
  control = NULL,
  constraints = NULL,
  residuals = FALSE,
  residual_vars = NULL,
  protocol = FALSE,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- definitions:

  a character vector of TDA assignments, the last defining `fn`, the
  quantity to minimise – or an unevaluated
  [`{ }`](https://rdrr.io/r/base/Paren.html) block of plain R
  assignments instead, translated automatically the same way
  [`tda_fml`](https://janmarvin.github.io/TDA/reference/tda_fml.md)'s
  `definitions` is (see its Details for exactly what gets translated).

- data:

  a data frame; its columns are the variables the definitions may refer
  to.

- start:

  optional named vector of starting values.

- control:

  convergence settings from
  [`tda_control`](https://janmarvin.github.io/TDA/reference/tda_control.md).

- constraints:

  optional linear constraints on the parameters – `freg`'s `con=`, the
  identical mechanism and `bN`-by-position convention as
  [`tda_fml`](https://janmarvin.github.io/TDA/reference/tda_fml.md)'s
  `constraints` (see there for the full explanation and a worked example
  of the naming trap); it works the same way here.

- residuals:

  ask TDA to also compute, per case, its contribution to `fn` at the
  converged parameters – `freg`'s own `pres=`, the identical mechanism
  as [`tda_fml`](https://janmarvin.github.io/TDA/reference/tda_fml.md)'s
  own `residuals` (see there for what it actually is: not a classical
  observed-minus-fitted residual). For `fn = r*r`, say, this is that
  squared residual's value, not `r` itself.

- residual_vars:

  with `residuals = TRUE`, extra columns to write alongside `fn`'s
  value, one per case – `freg`'s own `v=`, the same mechanism as
  [`tda_fml`](https://janmarvin.github.io/TDA/reference/tda_fml.md)'s
  `residual_vars`.

- protocol:

  ask TDA to also write its iteration-by-iteration diagnostic log –
  `freg`'s `prot=`, the same mechanism as
  [`tda_fml`](https://janmarvin.github.io/TDA/reference/tda_fml.md)'s
  `protocol`.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

An object of class `tda_fit`. With `residuals = TRUE`, also carries
`$residuals`.

## Details

[`tda_nlreg`](https://janmarvin.github.io/TDA/reference/tda_nlreg.md)'s
`expr` argument also fits a user-written nonlinear function now, using
`nlreg`'s dedicated algorithm (proper standard errors, and orthogonal
distance regression as an option) rather than a general-purpose
minimiser – prefer it for an ordinary nonlinear least-squares fit.
`tda_freg` stays the more flexible tool for the cases `nlreg` cannot
express at all: any summed objective, not just a sum of squared
residuals (a robust loss, a weighted or constrained one, anything
written by hand).

## Writing the function

`fn` is the quantity summed over cases and minimised, not the regression
function itself. Least squares is therefore the squared residual:


    tda_freg(c("r = y - a * exp(b * x)", "fn = r*r"), d,
             start = c(a = 1, b = 0.1))

## Algorithm

TDA's default minimiser is Newton, which overflows
[`exp()`](https://rdrr.io/r/base/Log.html) on the first step for
functions of this shape. `tda_freg` therefore defaults to BFGS, which is
stable and reaches the same estimates as `nls`. Pass
`control = tda_control(algorithm = 5)` for TDA's default.

## See also

Other regression:
[`TDA_FAMILIES`](https://janmarvin.github.io/TDA/reference/tda_glm.md),
[`TDA_QRMODELS`](https://janmarvin.github.io/TDA/reference/tda_qreg.md),
[`tda_gdf()`](https://janmarvin.github.io/TDA/reference/tda_gdf.md),
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
d <- data.frame(x = 1:20, y = 3 * exp(0.1 * (1:20)) + rnorm(20, 0, 0.5))
tda_freg(c("r = y - a * exp(b * x)", "fn = r*r"), d,
         start = c(a = 1, b = 0.1))
#> Call: tda_freg(definitions = c("r = y - a * exp(b * x)", "fn = r*r"), 
#>     data = d, start = c(a = 1, b = 0.1))
#> 
#> Cases: 20
#> Converged in 17 iterations
#> 
#>  Idx Parameter  Coeff  Error C/Error Signif
#>    1         a 3.0067 0.1059 28.3857 1.0000
#>    2         b 0.1005 0.0021 46.7940 1.0000

# the same model, written in R syntax instead
tda_freg({
    r = y - a * exp(b * x)
    fn = r * r
}, d, start = c(a = 1, b = 0.1))
#> Call: tda_freg(definitions = {
#>     r = y - a * exp(b * x)
#>     fn = r * r
#> }, data = d, start = c(a = 1, b = 0.1))
#> 
#> Cases: 20
#> Converged in 17 iterations
#> 
#>  Idx Parameter  Coeff  Error C/Error Signif
#>    1         a 3.0067 0.1059 28.3857 1.0000
#>    2         b 0.1005 0.0021 46.7940 1.0000
```
