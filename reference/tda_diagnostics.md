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
