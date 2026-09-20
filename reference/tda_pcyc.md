# Cycle structure of permutations

TDA's `pcyc`: each record of the given integer variables is a
permutation; the command reports its cycle decomposition. Verified in
the suite against a direct R cycle count.

## Usage

``` r
tda_pcyc(perm, ...)
```

## Arguments

- perm:

  a matrix or data frame, one permutation per row.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

the output file's lines (one cycle report per record).

## Examples

``` r
# each row is a permutation of 1..4, written as where each element
# goes; the last column of the output is its cycle decomposition
tda_pcyc(rbind(c(2, 3, 1, 4),    # (1,2,3)(4)
               c(1, 2, 3, 4),    # the identity, four fixed points
               c(2, 1, 4, 3)))   # (1,2)(3,4)
#> [1] "   1    4    2    0    2    0    2    3    1    4   (1,2,3)(4)"  
#> [2] "   2    4    4    0    0    0    1    2    3    4   (1)(2)(3)(4)"
#> [3] "   3    4    2    0    2    0    2    1    4    3   (1,2)(3,4)"  
```
