# Union of two graphs

TDA's `gdu`: graph 1 is loaded from `edges1` (with the `perm` option of
`gdd`), then the union with the edge list `edges2` is written to the
output file.

## Usage

``` r
tda_gdu(edges1, edges2, missing_value = -1, ...)
```

## Arguments

- edges1, edges2:

  data frames: from, to, value.

- missing_value:

  substitute written for missing values.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

the output file's lines.

## Examples

``` r
e1 <- data.frame(i = c(1, 2), j = c(2, 3), v = c(1, 1))
e2 <- data.frame(i = c(1, 3), j = c(3, 1), v = c(2, 2))
tda_gdu(e1, e2)   # value1/value2 side by side, -1 where absent
#>   from to value1 value2
#> 1    1  2      1     -1
#> 2    1  3     -1      2
#> 3    2  3      1     -1
#> 4    3  1     -1      2
```
