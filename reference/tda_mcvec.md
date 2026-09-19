# Vectorizing a matrix and its inverse

`tda_mcvec` stacks the columns of `A` into a single vector – TDA's
`mcvec`, the classic matrix-algebra `vec(A)` operator, matching base R's
`as.vector(A)` (R matrices are column-major already). `tda_mrvec` stacks
the rows instead – TDA's `mrvec`, matching `as.vector(t(A))`.
`tda_mivec` is the inverse of `tda_mcvec`: given a vector and a row
count `n`, it rebuilds an n-row matrix by filling column by column –
TDA's `mivec`, matching base R's `matrix(x, nrow = n)`.

## Usage

``` r
tda_mcvec(A, dir = tempfile("tda"))

tda_mrvec(A, dir = tempfile("tda"))

tda_mivec(x, n, dir = tempfile("tda"))
```

## Arguments

- A:

  a numeric matrix.

- dir:

  working directory.

- x:

  a numeric vector, whose length must divide evenly by `n`.

- n:

  the number of rows for `tda_mivec`'s result.

## Value

`tda_mcvec`/`tda_mrvec` return a plain numeric vector; `tda_mivec`
returns a matrix.

## See also

Other matrix reshaping: [`tda_mag()`](tda_mag.md),
[`tda_mcath()`](tda_mcath.md), [`tda_mrsum()`](tda_mrsum.md),
[`tda_msort()`](tda_msort.md), [`tda_msrow()`](tda_msrow.md),
[`tda_mtrim()`](tda_mtrim.md)

## Examples

``` r
A <- matrix(1:6, 2, 3)
tda_mcvec(A)               # == as.vector(A)
#> [1] 1 2 3 4 5 6
tda_mrvec(A)               # == as.vector(t(A))
#> [1] 1 3 5 2 4 6
tda_mivec(1:6, 2)          # == matrix(1:6, nrow = 2)
#>      [,1] [,2] [,3]
#> [1,]    1    3    5
#> [2,]    2    4    6
```
