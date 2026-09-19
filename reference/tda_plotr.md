# Grey-scale relief of a matrix

TDA's `plotr`.

## Usage

``` r
tda_plotr(z, grey_range = c(0, 1), ...)
```

## Arguments

- z:

  value matrix.

- grey_range:

  two greys mapped to the value extremes.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

path of the PostScript file, invisibly.

## Examples

``` r
r <- tda_plotr(outer(1:5, 1:5, "+"))
file.exists(grep("[.]ps$", r, value = TRUE))   # the finished plot
#> [1] TRUE
```
