# Integrate, minimise, and evaluate a function

These take a function written in TDA's expression language, with `x` as
the formal argument, and evaluate it inside TDA. They are not given an R
function, because the expression is evaluated by TDA's interpreter.

## Usage

``` r
tda_integrate(expr, from, to, options = list(), dir = tempfile("tda"))

tda_range(
  expr,
  start = NULL,
  use_derivatives = FALSE,
  max_boxes = NULL,
  max_iter = NULL,
  tol_width = NULL,
  tol_fd = NULL,
  tol_fe = NULL,
  options = list(),
  dir = tempfile("tda")
)

tda_evalf(expr, options = list(), derivatives = 0, dir = tempfile("tda"))

tda_minimize(expr, options = list(), dir = tempfile("tda"))

tda_gmin(
  expr,
  start = NULL,
  use_derivatives = FALSE,
  max_boxes = NULL,
  max_iter = NULL,
  tol_width = NULL,
  tol_fd = NULL,
  tol_fe = NULL,
  options = list(),
  dir = tempfile("tda")
)
```

## Arguments

- expr:

  a single string, e.g. `"x*x"` or `"exp(-x)*sin(x)"`.

- from, to:

  the limits of integration.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- start:

  for `tda_range` and `tda_gmin`, a starting value for each parameter,
  as a list: `list(2)` for one parameter starting at 2 with TDA's
  default search box (`value-1` to `value+1`); `list(c(2, 0, 10))` for
  the same starting value with an explicit box, `[0, 10]`. TDA's `xp=`.
  (Beware writing the raw comma syntax by hand: TDA reads each further
  number as another parameter's starting value, not a box width for the
  one before it; the list form here avoids that.)

- use_derivatives:

  for `tda_range` and `tda_gmin`, use the function's derivatives during
  the search (`ns=1`); TDA's default does not.

- max_boxes, max_iter, tol_width, tol_fd, tol_fe:

  for `tda_range` and `tda_gmin`, the box search's limits and tolerances
  (`nbox=`/`mxit=`/`tolbw=`/`tolfd=`/`tolfe=`). Not the same search as
  `tda_minimize`'s [`tda_control`](tda_control.md)-style minimiser, so
  `tda_control` does not apply to either of these.

- derivatives:

  for `tda_evalf`, 0 for the value alone (`evalf`), 1 for the gradient
  as well (`evalf1`), 2 for the Hessian too (`evalf2`); TDA computes
  them analytically.

## Value

An object with the numeric result in `value` (`tda_integrate`,
`tda_minimize`, `tda_gmin`), the two bounds in `range` (`tda_range`), or
nothing structured (`tda_evalf` – see Details);
`tda_minimize`/`tda_gmin` also carry `estimates`, the parameter values
at the minimum, and the run is always in `run`.

## Details

`tda_evalf` does not report a value at all: `xf` is a plotting command,
drawing the function over the range `options$x` names (default 1 to 10)
to a file called `xplot.ps` in the run's directory – read it with
`tda_read_ps(x$run, which = "xplot.ps")`, the same way as any other TDA
plot.

## See also

Other smoothing: [`tda_interp()`](tda_interp.md),
[`tda_isotonic()`](tda_isotonic.md), [`tda_mat()`](tda_mat.md),
[`tda_smd()`](tda_smd.md), [`tda_spl()`](tda_spl.md)

## Examples

``` r
tda_integrate("x*x", 0, 1)$value      # 1/3
#> [1] 0.3333333
tda_range("sin(x)")$range             # [-1, 1] over its default domain
#>      lower      upper 
#> -0.8414545  0.8414545 
tda_minimize("(x - 2)^2", options = list(xp = "0"))$estimates
#>   Idx Parameter Value Error Value_E Signif
#> 1   1         x     2    NA      NA     NA

# gmin: a box-search-based *global* minimum, not tda_minimize's
# local search -- its start=/tolerances, not tda_control(). Its
# interval arithmetic does not accept ^, unlike tda_minimize's
# expression -- (x-2)*(x-2), not (x-2)^2
tda_gmin("(x-2)*(x-2)", start = list(c(2, 0, 4)))$value
#> [1] 0

# evalf: the function value at a point, 1.1342
tda_evalf("sin(x) + x^2", options = list(x = 0.7))$value
#> [1] 1.134218
```
