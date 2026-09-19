# Declare a block structure

TDA's standalone `dblock` on the current-run data: blocks are maximal
runs of identical values. Mostly useful inside longer
[`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md)
pipelines; for block-mode variable creation put `dblock=` inside `nvar`
(see the package tests).

## Usage

``` r
tda_dblock(data, by, commands = character(), ...)
```

## Arguments

- data:

  data frame.

- by:

  column name defining blocks.

- commands:

  further command lines run under the block structure. Column names are
  TDA-sanitised (upper-cased) first, so refer to them in that form.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

the printed output, invisibly.

## Examples

``` r
out <- tda_dblock(data.frame(G = c(1, 1, 2), X = 1:3), by = "G",
                  commands = "dstat = X;")
# blocks are defined by G; dstat runs under the block structure
grep("^X", out, value = TRUE)
#> [1] "X            1.0000     3.0000     2.0000     1.0000           6.0000"
```
