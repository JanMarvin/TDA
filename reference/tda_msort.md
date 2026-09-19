# Sort or rank the rows of a matrix

`tda_msort` sorts the rows of `X` ascending by the columns named in
`by`, first column first – TDA's `msort`, matching
`X[do.call(order, as.data.frame(X[, by])), ]`. `tda_msort1` does the
same and then drops consecutive duplicate rows (comparing every column,
not just `by`) – TDA's `msort1`, matching `unique(sorted_X)`.
`tda_mrank` – despite the name – is *not* a linear-algebra matrix rank
(no TDA command returns one as a matrix value;
[`tda_mginv`](tda_minvs.md) prints a pseudorank as a side effect, but
never returns it). It returns the sort *permutation* itself: an integer
vector `p` such that `X[p, ]` is sorted by `by` – TDA's `mrank`,
confirmed against a hand-checked instance to equal base R's
`order(X[, by[1]], X[, by[2]], ...)`, not
[`rank()`](https://rdrr.io/r/base/rank.html).

## Usage

``` r
tda_msort(X, by, dir = tempfile("tda"))

tda_msort1(X, by, dir = tempfile("tda"))

tda_mrank(X, by, dir = tempfile("tda"))
```

## Arguments

- X:

  a numeric matrix.

- by:

  which columns to sort by, in order (1-based, at least one).

- dir:

  working directory.

## Value

`tda_msort`/`tda_msort1` return a matrix; `tda_mrank` returns an integer
vector, the same length as `nrow(X)`.

## See also

Other matrix reshaping: [`tda_mag()`](tda_mag.md),
[`tda_mcath()`](tda_mcath.md), [`tda_mcvec()`](tda_mcvec.md),
[`tda_mrsum()`](tda_mrsum.md), [`tda_msrow()`](tda_msrow.md),
[`tda_mtrim()`](tda_mtrim.md)

## Examples

``` r
X <- rbind(c(30, 1), c(10, 1), c(20, 1), c(40, 1))
tda_msort(X, 1)
#>      [,1] [,2]
#> [1,]   10    1
#> [2,]   20    1
#> [3,]   30    1
#> [4,]   40    1
tda_mrank(X, 1)   # == order(X[, 1]), i.e. c(2, 3, 1, 4)
#> [1] 2 3 1 4
```
