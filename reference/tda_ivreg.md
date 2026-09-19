# Interval regression: a coefficient identified only up to bounds

`tda_ivreg`, `tda_ivreg1` and `tda_ivreg2` are TDA's three variants of
interval regression – both response and regressor are interval-valued –
and `tda_ivls` solves the same kind of problem as a system of linear
interval equations. All four return bounds on beta rather than a point
estimate, so [`coef()`](https://rdrr.io/r/stats/coef.html) would be
misleading; use [`tda_bounds`](tda_bounds.md).

## Usage

``` r
tda_ivreg(
  formula,
  data,
  method = c("two_step", "heuristic", "exact"),
  max_boxes = NULL,
  max_iter = NULL,
  tol_width = NULL,
  tol_beta = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)

tda_ivreg1(
  formula,
  data,
  search_box = list(alpha = c(-10, 10), alpha_radius = c(0, 5), beta = c(-10, 10),
    beta_radius = c(0, 5)),
  max_boxes = 4000,
  max_iter = 4000,
  tol_width = NULL,
  tol_beta = NULL,
  tol_min = 0.001,
  options = list(),
  dir = tempfile("tda"),
  ...
)

tda_ivreg2(
  formula,
  data,
  method = c("direct", "search", "minimizer", "contour"),
  options = list(),
  dir = tempfile("tda"),
  ...
)

tda_ivls(formula, data, options = list(), dir = tempfile("tda"), ...)
```

## Arguments

- formula:

  a formula using `iv()` for interval-valued variables.

- data:

  a data frame.

- method:

  for `tda_ivreg`: `"two_step"` (default, only the first two steps),
  `"heuristic"`, or `"exact"` optimization – `ivreg`'s `opt=`. For
  `tda_ivreg2`: `"direct"` (default, direct calculation), `"search"`
  (direct search), `"minimizer"` (TDA's general minimiser), or
  `"contour"` – the same `opt=` option with a different meaning for each
  command. The choice changes the fit: on the same data, `"exact"`
  converges to a different beta than `"two_step"` and can exhaust the
  box search where the default does not.

- max_boxes, max_iter, tol_width, tol_beta:

  for `tda_ivreg` and `tda_ivreg1`, the box search's iteration limits
  and tolerances (`nbox=`/`mxit=`/`tolbw=`/`tolf=`);
  [`tda_control`](tda_control.md) does not apply here, these are
  specific to the interval-regression box search, not TDA's general
  minimiser. For `tda_ivreg1` the limits default to 4000 each – TDA's
  100/100 rarely completes a real search.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

- search_box:

  for `tda_ivreg1`, the domain the box search covers: a named list of
  `c(lower, upper)` ranges for `alpha`, `alpha_radius`, `beta` and
  `beta_radius` (equivalently a 4x2 matrix, or 8 numbers, four lower
  then four upper). The solution must lie inside it; widen it if the
  fitted bounds hug an edge. The radii cannot go below 0.

- tol_min:

  for `tda_ivreg1`, the certification tolerance (`tolfe=`): a box is
  accepted when its inclusion lower bound is within this of the best
  function value found. See Details for why TDA's default is unusable.

## Value

An object of class `tda_fit`.

## Details

Interval-valued data is a lower and an upper bound per case, a different
thing from interval *censoring*, and has no equivalent in base R or in
survival. `iv(lo, hi)` marks such a variable in a formula, the way
`Surv()` marks an episode.

`ivreg1` and `ivreg2` are global-optimization (branch-and-bound)
searches, and each reports its bounds in its own format – `ivreg2` a
centre and a half-width (“Beta (center-radius): C R”, read as
`[C - R, C + R]`), and `ivls` a bare “Final bounds” with no label at
all.

`ivreg1` needs a word of its own. As shipped, the command's search
domain is hard-coded – alpha in `[1, 2.5]`, beta in `[0, 0.1]`,
evidently the dataset it was being developed against – so for almost any
other data the answer is simply not inside the searched box and TDA
reports “0 finally accepted boxes”. The command was left mid-study, not
broken: on data whose solution lies inside the box, the search certifies
cleanly. `tda_ivreg1` therefore passes the domain in through the
package's `sbox=` option (a guarded addition to the C; the standalone
program is unchanged), always requests the `prot=` protocol file, and
reads the accepted boxes back: `boxes` holds every box the protocol
reports, `parameters` the hull of the near-optimal ones over (alpha,
alpha radius, beta, beta radius), and `beta`/`alpha` – what
[`tda_bounds`](tda_bounds.md) returns – the centre rows of that hull.
Certification follows TDA's rule: a box counts when its inclusion lower
bound is within `tol_min` of the best value found. The shipped default
of `1e-10` for that tolerance is unreachable (interval arithmetic
over-estimates a box's range), which is the second reason the command
always seemed to fail; `tol_min` defaults to `1e-3` here instead.

## See also

Other interval-valued data: [`tda_bounds()`](tda_bounds.md),
[`tda_idf()`](tda_idf.md), [`tda_ilsreg()`](tda_ilsreg.md),
[`tda_imean()`](tda_imean.md), [`tda_imreg()`](tda_imreg.md),
[`tda_inpreg()`](tda_inpreg.md), [`tda_ivar1()`](tda_ivar1.md),
[`tda_sddf()`](tda_sddf.md)

## Examples

``` r
set.seed(1)
d <- data.frame(xlo = 1:20)
d$xhi <- d$xlo + 2
d$ylo <- 2 + 0.5 * d$xlo + rnorm(20)
d$yhi <- d$ylo + 1.5
tda_ivreg(iv(ylo, yhi) ~ iv(xlo, xhi), d)
#> Call: tda_ivreg(iv(ylo, yhi) ~ iv(xlo, xhi), d)
#> 
#> Cases: 20 
#> Mean of the response: [8.040524, 8.340524] 
#> Mean of the regressor: [10.5, 12.5] 
#> Variance of the regressor: [24.225, 44.25] 
#> 
#> Bounds on the slope: [0.1868326, 0.8283934]  (outer interval)

# ivreg2: converges readily even where ivreg1 often does not
d2 <- data.frame(xlo = c(1, 2, 3, 4))
d2$xhi <- d2$xlo + 1
d2$ylo <- 2 + 0.5 * d2$xlo
d2$yhi <- d2$ylo + 1
tda_ivreg2(iv(ylo, yhi) ~ iv(xlo, xhi), d2)
#> Call: tda_ivreg2(iv(ylo, yhi) ~ iv(xlo, xhi), d2)
#> 
#> Cases: 4 
#> 
#> Bounds on the slope: [0.5, 0.5]
#> Bounds on the intercept: [0, 0] 

# ivls: the same interval-equations idea, its bound format
tda_ivls(iv(ylo, yhi) ~ iv(xlo, xhi), d2)
#> Call: tda_ivls(iv(ylo, yhi) ~ iv(xlo, xhi), d2)
#> 
#> Cases: 4 
#> 
#> Bounds on the slope: [1.25, 1.25]
```
