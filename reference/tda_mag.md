# Aggregating a matrix by row and column groups

Sums `A`'s entries within each (row-group, column-group) block – TDA's
`mag`. `row_groups`/`col_groups` assign each row/column of `A` to a
group (any positive integers; groups need not be contiguous or start at
1). The result has one row per distinct value in `row_groups` and one
column per distinct value in `col_groups`, in ascending order of group
number.

## Usage

``` r
tda_mag(A, row_groups, col_groups, dir = tempfile("tda"))
```

## Arguments

- A:

  a numeric matrix.

- row_groups:

  integer group index for each row, length `nrow(A)`.

- col_groups:

  integer group index for each column, length `ncol(A)`.

- dir:

  working directory.

## Value

A numeric matrix.

## See also

Other matrix reshaping:
[`tda_mcath()`](https://janmarvin.github.io/TDA/reference/tda_mcath.md),
[`tda_mcvec()`](https://janmarvin.github.io/TDA/reference/tda_mcvec.md),
[`tda_mrsum()`](https://janmarvin.github.io/TDA/reference/tda_mrsum.md),
[`tda_msort()`](https://janmarvin.github.io/TDA/reference/tda_msort.md),
[`tda_msrow()`](https://janmarvin.github.io/TDA/reference/tda_msrow.md),
[`tda_mtrim()`](https://janmarvin.github.io/TDA/reference/tda_mtrim.md)

## Examples

``` r
A <- matrix(1:12, 3, 4, byrow = TRUE)
tda_mag(A, row_groups = c(1, 2, 1), col_groups = c(1, 1, 2, 2))
#>      [,1] [,2]
#> [1,]   22   30
#> [2,]   11   15
```
