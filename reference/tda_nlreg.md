# Non-linear regression

`tda_nlreg` is TDA's `nlreg` command.

## Usage

``` r
tda_nlreg(
  formula,
  data,
  expr = NULL,
  start = NULL,
  residuals = FALSE,
  control = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)

# S3 method for class 'tda_nlreg'
predict(object, newdata = NULL, ...)
```

## Arguments

- formula:

  a two-sided formula. With `expr`, only the response and the predictor
  names it uses need to be present here – the functional form itself
  comes from `expr`, not from the formula's own right-hand side.

- data:

  a data frame.

- expr:

  optional nonlinear expression in TDA's language, written against
  `formula`'s predictor names – e.g. `"a * exp(b * x)"` for \\y = a
  e^{bx}\\. Names on the right that are neither predictors nor defined
  elsewhere are the parameters to be estimated. Without this, `nlreg`
  fits an ordinary (or orthogonal distance) linear predictor instead.

- start:

  optional named list or vector of starting values, in the order `expr`
  introduces its parameters. TDA's default starting point is all `1`s,
  which will not converge for every functional form.

- residuals:

  for `tda_nlreg` and `tda_lsreg1`, ask TDA to write its per-case table
  (`pres=`) and return it as `residuals`: the case number, the response,
  the fitted value, the residual, then each predictor. TDA does not
  compute it unless asked.

- control:

  convergence settings from
  [`tda_control`](https://janmarvin.github.io/TDA/reference/tda_control.md).

- options:

  a named list of further TDA options, passed through: `opt = 1`
  (default) is ordinary least squares, `opt = 2` orthogonal distance
  regression.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

An object of class `tda_fit`.

## Non-linear fits

`nlreg` takes an optional right-hand-side expression and starting values
(`xp=`) to fit a non-linear function; without them it falls back to a
linear predictor. That right-hand-side expression is reachable through
`expr`/`start`: fitting \\y = a \\ e^{bx}\\ on data generated with
\\a=2, b=0.5\\ recovers both closely. `options = list(opt = 2)`
additionally selects orthogonal distance regression (errors in both
variables), which
[`tda_lsreg`](https://janmarvin.github.io/TDA/reference/tda_lsreg.md)
has no equivalent for.

## See also

Other regression:
[`TDA_FAMILIES`](https://janmarvin.github.io/TDA/reference/tda_glm.md),
[`TDA_QRMODELS`](https://janmarvin.github.io/TDA/reference/tda_qreg.md),
[`tda_freg()`](https://janmarvin.github.io/TDA/reference/tda_freg.md),
[`tda_gdf()`](https://janmarvin.github.io/TDA/reference/tda_gdf.md),
[`tda_l1reg()`](https://janmarvin.github.io/TDA/reference/tda_l1reg.md),
[`tda_lsreg()`](https://janmarvin.github.io/TDA/reference/tda_lsreg.md),
[`tda_mlrc_design()`](https://janmarvin.github.io/TDA/reference/tda_mlrc_design.md),
[`tda_mreg()`](https://janmarvin.github.io/TDA/reference/tda_mreg.md),
[`tda_npreg()`](https://janmarvin.github.io/TDA/reference/tda_npreg.md),
[`tda_zreg()`](https://janmarvin.github.io/TDA/reference/tda_zreg.md)

## Examples

``` r
set.seed(25)
d <- data.frame(x = round(rnorm(60, sd = 2), 3))
d$y <- round(1.5 + 0.9 * d$x + rnorm(60, sd = 1.5), 3)
reg <- tda_nlreg(y ~ x, d, options = list(opt = 1))  # ordinary least squares
coef(reg)
#> Intercept         x 
#> 1.2121419 0.7990406 
coef(tda_nlreg(y ~ x, d, options = list(opt = 2)))  # orthogonal distance
#> Intercept         x 
#>  1.328604  1.066289 
head(fitted(reg))     # no newdata needed
#> [1]  0.8733487 -0.4522596 -0.6312447  1.7259250 -1.1849799  0.5001967
predict(reg, newdata = data.frame(x = c(0, 1, 2)))
#> [1] 1.212142 2.011182 2.810223

# a nonlinear fit: y = a * exp(b * x), true a = 2, b = 0.5
set.seed(1)
d2 <- data.frame(x = round(runif(40, 0, 5), 2))
d2$y <- round(2 * exp(0.5 * d2$x) + rnorm(40, sd = 0.5), 3)
nl <- tda_nlreg(y ~ x, d2, expr = "a * exp(b * x)", start = c(1, 1))
coef(nl)
#>         a         b 
#> 2.0685764 0.4912572 
# predict()/fitted() for this case evaluate the fitted expression
# itself (in R, not TDA) rather than X %*% coefficients, which would
# be meaningless for a nonlinear model
head(fitted(nl))
#> [1]  3.975827  5.158259  8.430501 19.243404  3.397468 18.776489
predict(nl, newdata = data.frame(x = c(0, 1, 2, 5)))
#> [1]  2.068576  3.380819  5.525507 24.122540
```
