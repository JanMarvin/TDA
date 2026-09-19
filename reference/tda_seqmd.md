# Regression models for events

Turns sequence data into the discrete-time expanded form a logistic or
probit model needs: one row per individual per time point at risk, a
binary indicator of whether the target event (a transition from one
state to another) occurred there, dummy variables for each time period
(so period effects can be estimated freely, the way `factor(t)` would in
a regular model), and any covariates carried along. Fit the result with
[`tda_qreg`](https://janmarvin.github.io/TDA/reference/tda_qreg.md)
(`model = "logit"` or `"probit"`) or
[`tda_glm`](https://janmarvin.github.io/TDA/reference/tda_glm.md)`(family = binomial)`
on the `event` column against the period dummies and covariates – TDA's
manual (chapter 6.18) calls this approach out explicitly as the
discrete-time counterpart to
[`tda_rate`](https://janmarvin.github.io/TDA/reference/tda_rate.md).

## Usage

``` r
tda_seqmd(
  sequences,
  event,
  data = NULL,
  covariates = NULL,
  event_covariates = NULL,
  select = NULL,
  tp = NULL,
  summary = FALSE,
  options = list(),
  dir = tempfile("tda"),
  ...
)

tda_seqev(
  sequences,
  data = NULL,
  sn = 1L,
  select = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)

tda_seqevd(
  sequences,
  data = NULL,
  sn = 1L,
  select = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- sequences:

  a data frame or matrix, one row per case, one column per time point,
  in order – the same shape
  [`tda_seqm`](https://janmarvin.github.io/TDA/reference/tda_seqm.md)
  takes – or, given `data` too, a character vector naming those columns
  within it directly, rather than requiring a separately pre-subsetted
  data frame: real sequence data typically carries an ID column and
  covariates alongside the states in the same data frame (TDA's
  `examples/exam/seq.d4`, e.g.), not only the states on their own. For
  `tda_seqev`/`tda_seqevd`, also a *list* of several such specs, to
  define several independent sequence data structures at once – TDA's
  `seqdef(sn=1)`, `seqdef(sn=2)`, ... – the pattern TDA's own manual
  uses to build a second sequence structure from `seq.d4`'s `S0..S5`
  columns, entirely separate from the first (`Y0..Y5`), not a covariate
  of it. Use `sn=` to say which one a given call analyzes.

- event:

  for `tda_seqmd`, the target transition, as `c(from, to)` – the two
  states such that going from the first to the second at time `t` counts
  as the event.

- data:

  optional data frame `sequences` (and, for `tda_seqmd`, `covariates`)
  name columns within, instead of `sequences`/`covariates` being the
  actual data already.

- covariates:

  optional further time-independent covariates, carried into the output
  unchanged – a data frame or matrix, one row per case in the same order
  as `sequences`, or (with `data` given) a character vector naming
  columns within it.

- event_covariates:

  optional time-varying covariates – TDA's `xe=`, different from
  `covariates`/`v=`: a separate value for each time point, rather than
  one fixed value per case. A named list, needs `data`: each element the
  per-time-point column names for one such covariate, in the same order
  as `sequences` itself, named by what the result's column should be
  called – `list(price = paste0("Price", 0:5))`, say, for a covariate
  stored one column per time point the same way the states themselves
  are. **Exactly how TDA derives the value it puts in the output is not
  documented** – it is not simply the value at the event's time point;
  check against your data before relying on a particular reading of this
  column. say, for a covariate stored one column per time point the same
  way the states themselves are.

- select:

  optional expression (`sel=`) restricting which sequences are used –
  written with R's comparison operators (`Y0 == 1`, or `Y0 = 1`) or
  TDA's function form (`eq(Y0,1)`) interchangeably; `&`/`|` combine
  several conditions the same way in both.

- tp:

  optional restriction of the time axis, as a vector of time points or a
  single TDA range expression such as `"0(1)10"`.

- summary:

  instead of the full per-case, per-period data, just a
  `time`/`riskset`/`events` table – TDA's `seqmd(ev=...)` called with no
  output file at all, a separate mode in which TDA prints this table to
  the console instead of writing any records. `covariates`/
  `event_covariates` are not accepted together with this, since TDA's
  console summary carries no covariate columns to show.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

- sn:

  for `tda_seqev`/`tda_seqevd`, which sequence data structure to use
  when `sequences` defines more than one (see above) – TDA's `sn=`, 1 by
  default, the same as when only one structure was ever defined.

## Value

`tda_seqmd` returns a data frame: `id`, `time`, `event` (the binary
indicator), one `period` dummy column per distinct time point, and any
`covariates` – or, with `summary = TRUE`, `time`, `riskset`, `events`
instead. `tda_seqev` returns a data frame, one row per transition that
occurs in the data: `from`, `to`, `count`. `tda_seqevd` returns a data
frame: `time`, `cases`, one match-count column per event type found, and
`total`.

## Details

`tda_seqev` and `tda_seqevd` are the exploratory step before committing
to a target event: `tda_seqev` counts how often each state-to-state
transition actually occurs in the data, and `tda_seqevd` tabulates event
counts over time.

**A note on names**: the manual documents this command as `seqnd`, but
TDA 6.4 itself calls it `seqmd`; `tda_seqmd()` uses the name that works.

This wrapper always defines exactly one sequence data structure per call
(a fresh `seqdef`, every time), so TDA's `sn=` – selecting among several
already-defined structures – has nothing to select among here: `sn = 1`
is always correct, and there is no way to reach any other value through
this interface. Not exposed as an argument for that reason, rather than
included as a parameter that could never meaningfully differ from its
default.

`tda_seqev` ignores TDA's `dtda=` (no description file is written even
when asked for one); its output shape (`from`/`to`/`count`) needs no
further description to read reliably. `tda_seqmd`'s `dtda=` is used
internally to read TDA's real column structure back, since
`event_covariates`, in particular, can change how many columns come back
and where.

## See also

Other sequence analysis:
[`tda_evalfi()`](https://janmarvin.github.io/TDA/reference/tda_evalfi.md),
[`tda_fml()`](https://janmarvin.github.io/TDA/reference/tda_fml.md),
[`tda_frml()`](https://janmarvin.github.io/TDA/reference/tda_frml.md),
[`tda_seq_info()`](https://janmarvin.github.io/TDA/reference/tda_seq_info.md),
[`tda_seqgc()`](https://janmarvin.github.io/TDA/reference/tda_seqgc.md),
[`tda_seqm()`](https://janmarvin.github.io/TDA/reference/tda_seqm.md),
[`tda_seqpe()`](https://janmarvin.github.io/TDA/reference/tda_seqpe.md),
[`tda_seqpm()`](https://janmarvin.github.io/TDA/reference/tda_seqpm.md)

## Examples

``` r
# three individuals, six time points (0-5), states 1-3; the target
# event is a 1 -> 2 transition, with two time-independent covariates
s <- data.frame(Y0 = c(1, 1, 2), Y1 = c(1, 2, 1), Y2 = c(2, 1, 1),
                Y3 = c(1, 1, 2), Y4 = c(1, 1, 2), Y5 = c(3, 3, 3))
cov <- data.frame(v1 = c(1, 2, 2), v2 = c(-10, -15, -20))

tda_seqev(s)    # which transitions actually occur, and how often
#>   from to count
#> 1    1  2     3
#> 2    1  3     2
#> 3    2  1     3
#> 4    2  3     1
tda_seqevd(s)   # event counts over time
#>   time cases ev1_2 ev1_3 ev2_1 ev2_3 total
#> 1    1     3     1     0     1     0     2
#> 2    2     3     1     0     1     0     2
#> 3    3     3     1     0     1     0     2
#> 4    4     3     0     0     0     0     0
#> 5    5     3     0     2     0     1     3

d <- tda_seqmd(s, event = c(1, 2), covariates = cov)
d
#>    id time event period1 period2 period3 period4 period5 v1  v2
#> 1   1    1     0       1       0       0       0       0  1 -10
#> 2   1    2     1       0       1       0       0       0  1 -10
#> 3   1    4     0       0       0       0       1       0  1 -10
#> 4   1    5     0       0       0       0       0       1  1 -10
#> 5   2    1     1       1       0       0       0       0  2 -15
#> 6   2    3     0       0       0       1       0       0  2 -15
#> 7   2    4     0       0       0       0       1       0  2 -15
#> 8   2    5     0       0       0       0       0       1  2 -15
#> 9   3    2     0       0       1       0       0       0  2 -20
#> 10  3    3     1       0       0       1       0       0  2 -20
# fit it: a period dummy for every distinct time point, the way
# factor(time) would in an ordinary logistic regression
pd <- grep("^period", names(d), value = TRUE)
tda_qreg(stats::as.formula(paste("event ~", paste(c(pd, "v1", "v2"),
                                                   collapse = "+"))),
        d, model = "logit", intercept = FALSE)
#> Call: tda_qreg(formula = stats::as.formula(paste("event ~", paste(c(pd, 
#>     "v1", "v2"), collapse = "+"))), data = d, model = "logit", 
#>     intercept = FALSE)
#> 
#> Cases: 10
#> Model:   logit (TDA 1)
#> logLik (starting values): -6.931472
#> logLik:  -4.158883
#> Converged in 18 iterations
#> 
#>  Idx    Cat Term Variable    Coeff     Error C/Error Signif
#>    1 1.0000    X  period1  -0.0000    3.3665 -0.0000 0.0000
#>    2 1.0000    X  period2  -0.0000    3.7417 -0.0000 0.0000
#>    3 1.0000    X  period3  -0.0000    4.3970 -0.0000 0.0000
#>    4 1.0000    X  period4 -18.2029 6341.5209 -0.0029 0.0023
#>    5 1.0000    X  period5 -18.2029 6341.5209 -0.0029 0.0023
#>    6 1.0000    X       v1  -0.0000    4.0000 -0.0000 0.0000
#>    7 1.0000    X       v2  -0.0000    0.4619 -0.0000 0.0000

# the same data, but with an ID column and the states/covariates all
# together in one data frame, the more realistic shape -- sequences=
# and covariates= just name columns in data= directly
full <- cbind(id = 1:3, s, cov)
tda_seqmd(paste0("Y", 0:5), event = c(1, 2), data = full,
         covariates = c("v1", "v2"))
#>    id time event period1 period2 period3 period4 period5 v1  v2
#> 1   1    1     0       1       0       0       0       0  1 -10
#> 2   1    2     1       0       1       0       0       0  1 -10
#> 3   1    4     0       0       0       0       1       0  1 -10
#> 4   1    5     0       0       0       0       0       1  1 -10
#> 5   2    1     1       1       0       0       0       0  2 -15
#> 6   2    3     0       0       0       1       0       0  2 -15
#> 7   2    4     0       0       0       0       1       0  2 -15
#> 8   2    5     0       0       0       0       0       1  2 -15
#> 9   3    2     0       0       1       0       0       0  2 -20
#> 10  3    3     1       0       0       1       0       0  2 -20
```
