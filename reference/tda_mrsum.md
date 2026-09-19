# Row and column sums, as vectors or diagonal matrices

`tda_mrsum`/`tda_mcsum` are TDA's `mrsum`/ `mcsum`, matching base R's
`rowSums`/`colSums`. `tda_mdrow`/`tda_mdcol` return the same sums as a
diagonal matrix instead of a vector – TDA's `mdrow`/ `mdcol`, matching
`diag(rowSums(A))`/ `diag(colSums(A))`; useful for building a graph's
degree matrix from its adjacency matrix.

## Usage

``` r
tda_mrsum(A, dir = tempfile("tda"))

tda_mcsum(A, dir = tempfile("tda"))

tda_mdrow(A, dir = tempfile("tda"))

tda_mdcol(A, dir = tempfile("tda"))
```

## Arguments

- A:

  a numeric matrix.

- dir:

  working directory.

## Value

`tda_mrsum`/`tda_mcsum` return a plain numeric vector;
`tda_mdrow`/`tda_mdcol` return a diagonal matrix.

## See also

Other matrix reshaping:
[`tda_mag()`](https://janmarvin.github.io/TDA/reference/tda_mag.md),
[`tda_mcath()`](https://janmarvin.github.io/TDA/reference/tda_mcath.md),
[`tda_mcvec()`](https://janmarvin.github.io/TDA/reference/tda_mcvec.md),
[`tda_msort()`](https://janmarvin.github.io/TDA/reference/tda_msort.md),
[`tda_msrow()`](https://janmarvin.github.io/TDA/reference/tda_msrow.md),
[`tda_mtrim()`](https://janmarvin.github.io/TDA/reference/tda_mtrim.md)

## Examples

``` r
A <- matrix(1:6, 2, 3)
tda_mrsum(A); tda_mcsum(A)
#> [1]  9 12
#> [1]  3  7 11
tda_mdrow(A); tda_mdcol(A)
#>      [,1] [,2]
#> [1,]    9    0
#> [2,]    0   12
#>      [,1] [,2] [,3]
#> [1,]    3    0    0
#> [2,]    0    7    0
#> [3,]    0    0   11
```
