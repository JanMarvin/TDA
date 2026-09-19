# tdaR: Transition Data Analysis

An R interface to TDA 6.4, Goetz Rohwer and Ulrich Poetter's program for
the analysis of event history and transition data.

## Details

Rohwer began TDA in 1989; it grew inside Hans-Peter Blossfeld's
*Household Dynamics and Social Inequality* project at the European
University Institute in Florence and then at the University of Bremen,
reached maturity with version 5.7 (1994), and was redesigned as the 6.x
line with Ulrich Poetter as co-author, the last releases coming from
Rohwer's later years at the Ruhr-Universitaet Bochum (the account TDA's
manual preface gives). It became one of the standard tools of
quantitative life course research: transition rate models, product-limit
and life table estimation, sequence analysis, and a long tail of
supporting methods, all driven by plain-text command files. Blossfeld
and Rohwer's *Techniques of Event History Modeling* (1995, 2002) is
written around it. This package embeds TDA's original C sources – the
models are computed by TDA itself, unchanged – behind ordinary R
functions, formulas and data frames, so that its methods, and results
computed with it decades ago, stay reachable from a current environment.

**Where to start**, by task:

- Event histories:
  [`tda_rate`](https://janmarvin.github.io/TDA/reference/tda_rate.md)
  fits the transition rate models (exponential through Gompertz,
  Weibull, log-logistic, Cox and the rest of `TDA_MODELS`), with
  `Surv()`-style formulas;
  [`tda_ple`](https://janmarvin.github.io/TDA/reference/tda_ltb.md),
  [`tda_ltb`](https://janmarvin.github.io/TDA/reference/tda_ltb.md) and
  [`tda_km`](https://janmarvin.github.io/TDA/reference/tda_ltb.md)
  estimate survivor functions and life tables.

- Regression:
  [`tda_lsreg`](https://janmarvin.github.io/TDA/reference/tda_lsreg.md)
  (least squares, also censored and grouped variants),
  [`tda_glm`](https://janmarvin.github.io/TDA/reference/tda_glm.md),
  [`tda_qreg`](https://janmarvin.github.io/TDA/reference/tda_qreg.md)
  (logit/probit and their ordered and multinomial relatives),
  [`tda_loglin`](https://janmarvin.github.io/TDA/reference/tda_loglin.md).

- Sequences:
  [`tda_seqm`](https://janmarvin.github.io/TDA/reference/tda_seqm.md)
  computes optimal matching distances (the result is a `dist`, ready for
  [`hclust()`](https://rdrr.io/r/stats/hclust.html) or
  [`tda_cluster`](https://janmarvin.github.io/TDA/reference/tda_cluster.md));
  [`tda_seqmd`](https://janmarvin.github.io/TDA/reference/tda_seqmd.md)
  turns sequences into episode data.

- Description:
  [`tda_dstat`](https://janmarvin.github.io/TDA/reference/tda_dstat.md),
  [`tda_freq`](https://janmarvin.github.io/TDA/reference/tda_freq1.md),
  [`tda_corr`](https://janmarvin.github.io/TDA/reference/tda_dstat.md)
  and relatives.

- Anything else:
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md) runs
  any TDA command file against an R data frame, and
  [`tda_help`](https://janmarvin.github.io/TDA/reference/tda_help.md)
  shows TDA's manual entry for any command, by its TDA name –
  `tda_help("lsreg")`, not `"tda_lsreg"`.

Fitted models answer to the usual verbs –
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`predict()`](https://rdrr.io/r/stats/predict.html) where TDA computes
predictions – and
[`tda_estimates()`](https://janmarvin.github.io/TDA/reference/tda_estimates.md)
returns TDA's coefficient table as printed.

## The shape of the interface

The package's model is base R's statistical interface –
[`lm()`](https://rdrr.io/r/stats/lm.html),
[`glm()`](https://rdrr.io/r/stats/glm.html), survival – not a new
grammar. Every export is one of five shapes, and knowing which shape a
function is tells you how to call it before reading its page:

- Model fits:

  `tda_<cmd>(formula, data, ...)` – exactly the
  [`lm()`](https://rdrr.io/r/stats/lm.html) idiom, with `Surv()` on the
  left for episodes and `iv()` for interval-valued variables. They
  return a `tda_fit` answering to the standard verbs above.

- Descriptives:

  `tda_<cmd>(...)` – data first, like
  [`cor()`](https://rdrr.io/r/stats/cor.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html): a data frame, a
  matrix, or loose vectors. They return an object carrying `table`.

- Builder sessions:

  a session object built up by a family of verbs and run at the end, the
  `ggplot()` idea:
  [`tda_ps()`](https://janmarvin.github.io/TDA/reference/tda_ps.md) then
  `tda_pl_*()` for plots,
  [`tda_graph()`](https://janmarvin.github.io/TDA/reference/tda_graph.md)
  then `tda_g_*()` for graphs,
  [`tda_spatial()`](https://janmarvin.github.io/TDA/reference/tda_spatial.md)
  then `tda_sd_*()` for spatial data. Each verb takes the session as its
  first argument and returns it, so the steps chain with `|>`.

- Runners:

  [`tda_run()`](https://janmarvin.github.io/TDA/reference/tda_run.md)
  executes any TDA command file against an R data frame, with
  [`tda_nvar()`](https://janmarvin.github.io/TDA/reference/tda_write_data.md)
  and
  [`tda_block()`](https://janmarvin.github.io/TDA/reference/tda_write_data.md)
  to build the text; everything above is ultimately this.

- Readers:

  `tda_read_*()` bring TDA's file formats back into R – PostScript
  plots, spatial files, data archives.

Common argument grammar across all of them: `formula` and `data` come
first where they exist; `weights`, `censor` and `select` always mean the
same thing wherever they appear; `control` takes
[`tda_control()`](https://janmarvin.github.io/TDA/reference/tda_control.md)
for iteration settings; `options = list(...)` is the escape hatch that
reaches *any* TDA option under TDA's name, so nothing the program can do
is out of reach even where no named argument exists; and `dir` pins the
working directory when you want to keep the run's files.

Two namespaces coexist deliberately. Argument names are R's (`weights`,
not `cwt`; `max_iter`, not `mxit`), but each function keeps its TDA
*command* name as the suffix: `tda_ple`, not `tda_kaplan_meier`. That is
what keeps thirty years of TDA literature usable – any command file in
Blossfeld and Rohwer, and any entry
[`tda_help()`](https://janmarvin.github.io/TDA/reference/tda_help.md)
shows, maps 1:1 onto the function that wraps it.

**Missing values.** TDA itself has no NA: it stores a numeric missing
value (`msys`, -5 by default) and computes with it like any other
number. Its manual is direct about this (section 6.1.1): once the data
matrix is created, “TDA no longer makes any distinction between valid
and missing values”, and excluding them “by appropriate case selection
commands” is the user's responsibility. The model-fitting functions here
do that selection for you: incomplete cases are dropped before TDA sees
the data, like [`lm()`](https://rdrr.io/r/stats/lm.html), with a message
saying how many. The descriptive functions warn instead, since without a
formula there is no saying which cases a computation uses. Clean your
data first if -5 appearing as a value would be a surprise.

`inst/examples/ehhnew.R` is a full worked showcase, meant to be read
beside *Techniques of Event History Modeling*: every model from the
book's command files, written out in full.

## References

Rohwer, G. and Poetter, U., *TDA User's Manual*. Ruhr-Universitaet
Bochum. TDA's homepage: <https://www.stat.rub.de/tda.html>.

Blossfeld, H.-P. and Rohwer, G. (2002). *Techniques of Event History
Modeling: New Approaches to Causal Analysis*, 2nd ed. Mahwah, NJ:
Lawrence Erlbaum.

## See also

Useful links:

- <https://github.com/JanMarvin/TDA>

- <https://janmarvin.github.io/TDA)>

- <https://www.stat.rub.de/tda.html>

- Report bugs at <https://github.com/JanMarvin/TDA/issues>

## Author

**Maintainer**: Jan Marvin Garbuszus <jan.garbuszus@ruhr-uni-bochum.de>

Authors:

- Jan Marvin Garbuszus <jan.garbuszus@ruhr-uni-bochum.de>

- Goetz Rohwer \[copyright holder\]

- Ulrich Poetter \[copyright holder\]

## Examples

``` r
set.seed(1)
d <- data.frame(x = 1:40, z = rnorm(40))
d$y <- 2 + 0.5 * d$x + d$z + rnorm(40)
fit <- tda_lsreg(y ~ x + z, d)
coef(fit)
#> Intercept         x         z 
#> 2.2338633 0.4934024 1.2353156 

tda_help("lsreg")  # the command's entry in TDA's manual
#> lsreg: The lsreg command performs least squares regression, optionally
#> with constraints. The syntax is:
#> 
#>     lsreg(
#>         w=...,          case weight variable, def. no weigths
#>         ni=1,           if without intercept, def. ni=0
#>         lsecon=...,     equality constraints
#>         lsicon=...,     inequality constraints
#>         dgrp=...,       estimate dummy variables with constraints
#>         tfmt=...,       print format for results, def. tfmt=10.4
#>         s=1,            use robust covariance matrix, def. s=0
#>         df=...,         print data to an output file
#>         fmt=...,        print format for df option, def. 10.4
#>         dtda=...,       TDA description file for df option
#>         ppar=...,       print estimated coefficients to output file
#>         pcov=...,       print covariance matrix to output file
#>         pres=...,       print residuals to output file
#>         mfmt=...,       print format for pcov and pres option
#>         mplog=...,      write norm of residuals into matrix
#>         mppar=...,      write parameters into matrix
#>         mpcov=...,      write covariance matrix into matrix
#>     ) = varlist;
#> 
#> varlist is a comma-separated list of variables. The first variable in
#> varlist is interpreted as the dependent variable.
#> 
#> See also: l1reg, mpcov, mplog, mppar
```
