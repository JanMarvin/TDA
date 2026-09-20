# Diagnostics TDA printed during a run

Every error, warning, note or "cannot" line TDA printed, on stdout or
stderr, for a run or a fitted model. TDA reports problems by printing
them and carrying on, so a fit can come back with a full table and a
message about what it could not do; nothing here is dropped from
[`tda_output()`](https://janmarvin.github.io/TDA/reference/tda_estimates.md),
this is the short list.

## Usage

``` r
tda_diagnostics(x)
```

## Arguments

- x:

  a `tda_result` from
  [`tda_run()`](https://janmarvin.github.io/TDA/reference/tda_run.md) or
  any fitted model.

## Value

a character vector, empty when TDA reported nothing.

## Examples

``` r
# a clean run reports nothing
tda_diagnostics(tda_run("nvar(noc = 5, X = case);"))
#> character(0)

# dstat asks for a variable that was never defined: TDA prints the
# error and continues, so the run returns and the message is here
r <- tda_run(c("nvar(noc = 5, X = case);", "dstat = Y;"))
tda_diagnostics(r)
#> [1] "Syntax error or undefined variables."
#> [2] "Error: =Y"                           
#> [3] "Syntax error."                       

# a fit stopped after one iteration still returns a coefficient
# table; the problem is reported alongside it
d <- data.frame(x = c(-1, -0.5, 0, 0.5, 1, 1.5), y = c(0, 0, 1, 0, 1, 1))
f <- suppressWarnings(tda_qreg(y ~ x, d, control = tda_control(maxit = 1)))
tda_diagnostics(f)
#> [1] "Problem: reached max number of iterations."
```
