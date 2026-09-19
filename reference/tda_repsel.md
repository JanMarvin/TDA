# Repeat or select block cases

TDA's `repsel` with a matrix expression giving, per block, the
repetition count.

## Usage

``` r
tda_repsel(data, by, times, commands = character(), ...)
```

## Arguments

- data:

  data frame.

- by:

  block-defining column.

- times:

  matrix expression (TDA syntax) of dimension blocks x 1.

- commands:

  further command lines run under the replicated case structure. Column
  names are TDA-sanitised (upper-cased) first, so refer to them in that
  form.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

the printed output, invisibly.

## Examples

``` r
# blocks by G, block 1 twice and block 2 once: cases become
# X = 1,2,1,2,3, and the statistics confirm it (mean 1.8, sum 9)
out <- tda_repsel(data.frame(G = c(1, 1, 2), X = 1:3), by = "G",
                  times = "<2, 1>", commands = "dstat = X;")
grep("^X", out, value = TRUE)
#> [1] "X            1.0000     3.0000     1.8000     0.8367           9.0000"
```
