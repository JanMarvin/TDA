# Trim or pad a matrix's edges

Deletes or adds whole rows/columns at the edges of `A` – TDA's `mtrim`.
A positive count deletes that many rows or columns from that edge; a
negative count adds that many zero rows/columns instead.
`leading_cols`/`trailing_cols` act on the first and last columns,
`leading_rows`/ `trailing_rows` on the first and last rows.

## Usage

``` r
tda_mtrim(
  A,
  leading_cols = 0,
  leading_rows = 0,
  trailing_cols = 0,
  trailing_rows = 0,
  dir = tempfile("tda")
)
```

## Arguments

- A:

  a numeric matrix.

- leading_cols, trailing_cols, leading_rows, trailing_rows:

  integer counts (default 0); positive deletes, negative pads with
  zeros.

- dir:

  working directory.

## Value

A numeric matrix.

## See also

Other matrix reshaping:
[`tda_mag()`](https://janmarvin.github.io/TDA/reference/tda_mag.md),
[`tda_mcath()`](https://janmarvin.github.io/TDA/reference/tda_mcath.md),
[`tda_mcvec()`](https://janmarvin.github.io/TDA/reference/tda_mcvec.md),
[`tda_mrsum()`](https://janmarvin.github.io/TDA/reference/tda_mrsum.md),
[`tda_msort()`](https://janmarvin.github.io/TDA/reference/tda_msort.md),
[`tda_msrow()`](https://janmarvin.github.io/TDA/reference/tda_msrow.md)

## Examples

``` r
A <- matrix(1:4, 2)
tda_mtrim(A, leading_cols = 1)    # drop the first column
#>      [,1]
#> [1,]    3
#> [2,]    4
tda_mtrim(A, leading_cols = -1)   # add a zero column in front
#>      [,1] [,2] [,3]
#> [1,]    0    1    3
#> [2,]    0    2    4
```
