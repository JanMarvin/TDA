# Concatenating matrices

`tda_mcath` concatenates horizontally – TDA's `mcath`, matching `cbind`.
`tda_mcatv` concatenates vertically – TDA's `mcatv`, matching `rbind`.
`tda_mcathv` takes the matrix direct sum: the inputs placed as
block-diagonal blocks, zero elsewhere – TDA's `mcathv` (all three take
two or more matrices).

## Usage

``` r
tda_mcath(..., dir = tempfile("tda"))

tda_mcatv(..., dir = tempfile("tda"))

tda_mcathv(..., dir = tempfile("tda"))
```

## Arguments

- ...:

  two or more numeric matrices.

- dir:

  working directory.

## Value

A numeric matrix.

## See also

Other matrix reshaping: [`tda_mag()`](tda_mag.md),
[`tda_mcvec()`](tda_mcvec.md), [`tda_mrsum()`](tda_mrsum.md),
[`tda_msort()`](tda_msort.md), [`tda_msrow()`](tda_msrow.md),
[`tda_mtrim()`](tda_mtrim.md)

## Examples

``` r
A <- matrix(1:4, 2); B <- matrix(5:8, 2)
tda_mcath(A, B)    # == cbind(A, B)
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    3    5    7
#> [2,]    2    4    6    8
tda_mcatv(A, B)    # == rbind(A, B)
#>      [,1] [,2]
#> [1,]    1    3
#> [2,]    2    4
#> [3,]    5    7
#> [4,]    6    8
tda_mcathv(A, B)   # block-diagonal direct sum
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    3    0    0
#> [2,]    2    4    0    0
#> [3,]    0    0    5    7
#> [4,]    0    0    6    8
```
