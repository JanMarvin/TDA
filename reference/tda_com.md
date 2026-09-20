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

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

the output file's lines.

## Examples

``` r
# the six 2-element subsets of {0, 1, 2, 3}: a row number, then the members
tda_com("subsets", n = 4, m = 2)
#>      [,1] [,2] [,3]
#> [1,]    1    0    1
#> [2,]    2    0    2
#> [3,]    3    0    3
#> [4,]    4    1    2
#> [5,]    5    1    3
#> [6,]    6    2    3
# the five partitions of 4: 4, 3+1, 2+2, 2+1+1, 1+1+1+1
tda_com("partitions", n = 4)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    1    4    0    0    0
#> [2,]    2    2    3    1    0    0
#> [3,]    3    2    2    2    0    0
#> [4,]    4    3    2    1    1    0
#> [5,]    5    4    1    1    1    1
```
