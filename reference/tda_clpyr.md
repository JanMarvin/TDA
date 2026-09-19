# Pyramidal clustering

TDA's `clpyr` on a distance matrix.

## Usage

``` r
tda_clpyr(d, linkage = c("single", "complete"), ...)
```

## Arguments

- d:

  symmetric distance matrix.

- linkage:

  "single" or "complete" link.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

the printed output.

## Examples

``` r
r <- tda_clpyr(as.matrix(dist(c(0, 1, 5, 6))))
# the top of the pyramid joins everything: 4 objects spanning 1..4
tail(r$pyramid, 1)
#>    cluster left right size min max next_left next_right index
#> 10      10    9     8    4   4   1         0          0     6
```
