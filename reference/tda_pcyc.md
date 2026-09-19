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

  passed to [`tda_run`](tda_run.md).

## Value

the output file's lines (one cycle report per record).

## Examples

``` r
tda_pcyc(rbind(c(2, 1, 3), c(2, 3, 1)))
#> [1] "   1    3    2    0    1    0    2    1    3   (1,2)(3)"
#> [2] "   2    3    1    0    2    0    2    3    1   (1,2,3)" 
```
