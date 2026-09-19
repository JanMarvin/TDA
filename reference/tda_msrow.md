# Selecting or permuting rows and columns

`tda_msrow`/`tda_mscol` select (and may reorder or repeat) rows/columns
by index – TDA's `msrow`/`mscol`, matching `X[rows, ]`/`X[, cols]`.
`tda_mprow`/ `tda_mpcol`/`tda_mpsym` permute rows, columns, or both by
the same permutation – TDA's `mprow`/`mpcol`/ `mpsym`, matching
`A[p, ]`/`A[, p]`/`A[p, p]` – confirmed against all five by hand
(`mprow`/`mpcol`/ `mpsym` additionally have their exact convention
stated in `m_mperm`'s header, `b(i,j)=a(p(i),j)` etc.). `mprow`/`mpcol`
require a permutation of `1:nrow(A)`/`1:ncol(A)`, so unlike `msrow`/
`mscol` they cannot repeat or drop indices; `mpsym` additionally
requires `A` to be square.

## Usage

``` r
tda_msrow(X, rows, dir = tempfile("tda"))

tda_mscol(X, cols, dir = tempfile("tda"))

tda_mprow(A, p, dir = tempfile("tda"))

tda_mpcol(A, p, dir = tempfile("tda"))

tda_mpsym(A, p, dir = tempfile("tda"))
```

## Arguments

- X, A:

  a numeric matrix.

- rows, cols:

  integer indices (1-based; `tda_msrow`/ `tda_mscol` allow repeats and
  omissions).

- dir:

  working directory.

- p:

  an integer permutation of `1:nrow(A)` (`tda_mprow`, `tda_mpsym`) or
  `1:ncol(A)` (`tda_mpcol`).

## Value

A numeric matrix.

## See also

Other matrix reshaping: [`tda_mag()`](tda_mag.md),
[`tda_mcath()`](tda_mcath.md), [`tda_mcvec()`](tda_mcvec.md),
[`tda_mrsum()`](tda_mrsum.md), [`tda_msort()`](tda_msort.md),
[`tda_mtrim()`](tda_mtrim.md)

## Examples

``` r
A <- matrix(1:9, 3, 3)
tda_msrow(A, c(2, 1))       # == A[c(2, 1), ]
#>      [,1] [,2] [,3]
#> [1,]    2    5    8
#> [2,]    1    4    7
tda_mprow(A, c(2, 3, 1))    # == A[c(2, 3, 1), ]
#>      [,1] [,2] [,3]
#> [1,]    2    5    8
#> [2,]    3    6    9
#> [3,]    1    4    7
tda_mpsym(A, c(2, 3, 1))    # == A[c(2, 3, 1), c(2, 3, 1)]
#>      [,1] [,2] [,3]
#> [1,]    5    8    2
#> [2,]    6    9    3
#> [3,]    4    7    1
```
