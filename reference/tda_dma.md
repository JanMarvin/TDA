# Principal components, factor analysis, correspondence analysis

`dma`: seven related methods for reducing a data matrix to a handful of
derived scores, selected with `alg` – principal components (1, from the
raw data; 2, from a covariance or correlation matrix already computed),
factor analysis (3, 4, the same two starting points but with
eigenvectors rescaled as loadings), dual scaling and correspondence
analysis (5, 6, for a frequency table rather than a data matrix), and a
direct SVD-based projection (7). Verified against
[`prcomp`](https://rdrr.io/r/stats/prcomp.html) for algorithm 1:
eigenvalues match exactly, eigenvectors match up to sign (an ordinary,
harmless ambiguity in any eigendecomposition).

## Usage

``` r
tda_dma(
  data,
  alg = 1,
  preprocess = NULL,
  ns = NULL,
  variables = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- data:

  a data frame or matrix. For algorithms 1–4 and 7, a standard data
  matrix (rows are cases); for 2 and 4, a covariance or correlation
  matrix (square, symmetric); for 5 and 6, a frequency table.

- alg:

  which of the seven algorithms; see Description. Default 1, ordinary
  principal components.

- preprocess:

  preprocessing, for algorithms 1, 5 and 7 only – for `alg = 1`: `1`
  nothing (TDA's default), `2` mean centre, `3` standardize; for `alg` 5
  or 7: `2`/`3` convert to row or table relative frequencies.

- ns:

  keep at most this many components; all of them by default.

- variables:

  optional subset of `data`'s column names to use, instead of all of
  them.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).
  Further `dma` options can be given the same way, notably `opt=`
  (algorithm variant) and `pcf=` (print classification frequencies).

## Value

A list with `values` (the eigenvalues), `percent` (each eigenvalue's
share of the total, as TDA prints it), `vectors` (the eigenvectors or
factor loadings, one column per component, row names taken from `data`),
and `scores` (the derived score for each case, one column per component)
– `scores` is empty for algorithms where TDA does not calculate it (2
and 4).

## Details

TDA's `v=` (choosing a subset of variables from a larger data matrix)
needs variable names that start uppercase, the same rule as anywhere
else in TDA – lowercase R column names are displayed correctly when
declared but cannot be *referenced* again later in the same script,
including by `v=` itself, failing with a “Syntax error or undefined
variables” that gives no hint the name's case is the problem.
`variables` handles this automatically; there is no need to rename
columns yourself first.

## See also

Other clustering:
[`TDA_CLUSTER`](https://janmarvin.github.io/TDA/reference/tda_cluster.md),
[`TDA_MDS`](https://janmarvin.github.io/TDA/reference/tda_mds.md),
[`tda_conjoint()`](https://janmarvin.github.io/TDA/reference/tda_conjoint.md),
[`tda_cutree()`](https://janmarvin.github.io/TDA/reference/tda_cutree.md),
[`tda_pdatd()`](https://janmarvin.github.io/TDA/reference/tda_pdatd.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x1 = rnorm(30), x2 = rnorm(30))
d$x3 <- d$x1 + d$x2 + rnorm(30, sd = 0.3)   # correlated with both
pc <- tda_dma(d)
pc$values                 # most of the variance is in the first PC
#> [1] 2.70895046 0.70467635 0.01930123
pc$vectors                # loadings: x3 dominates it
#>          PC1         PC2        PC3
#> x1 0.4244330  0.72193373 -0.5465056
#> x2 0.3387373 -0.68633924 -0.6435802
#> x3 0.8397105 -0.08803484  0.5358509
head(pc$scores)
#>             PC1        PC2          PC3
#> [1,]  1.4142045 -1.5126624  0.246377700
#> [2,]  0.1011365  0.1970438  0.002808551
#> [3,] -0.4257493 -0.8481232  0.078019005
#> [4,]  1.9603115  1.0521724 -0.006699393
#> [5,] -1.3934880  1.2948640  0.025353851
#> [6,] -1.4786799 -0.2037191  0.083798313

# a subset of the variables, instead of all of them
tda_dma(d, variables = c("x1", "x2"))$values
#> [1] 0.8711004 0.6407060
```
