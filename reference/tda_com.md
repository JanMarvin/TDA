# Combinatorial patterns

TDA's `com`: enumerate tuples, sets, permutations or partitions into an
output file.

## Usage

``` r
tda_com(
  pattern = c("tuples", "subsets", "m_tuples", "permutations", "partitions",
    "partitions_m", "bipartitions"),
  n = 1,
  m = 1,
  ...
)
```

## Arguments

- pattern:

  one of "tuples", "subsets", "m_tuples", "permutations", "partitions",
  "partitions_m", "bipartitions".

- n, m:

  dimensions (n the base size, m where the pattern needs it).

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

the output file's lines.

## Examples

``` r
tda_com("permutations", n = 3)   # the 6 permutations of 0,1,2
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    1    2
#> [2,]    2    1    0    2
#> [3,]    3    1    2    0
#> [4,]    4    2    1    0
#> [5,]    5    2    0    1
#> [6,]    6    0    2    1
```
