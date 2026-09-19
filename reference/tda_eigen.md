# Eigenvalues of a matrix, via TDA

Runs TDA's `etest`, which computes the eigenvalues and eigenvectors of a
square matrix and reports the residual of each. It is a test of TDA's
eigen routines rather than a general facility –
[`eigen`](https://rdrr.io/r/base/eigen.html) is what you want for
ordinary work – but it is the way to check what TDA itself computes for
a matrix a model is about to be fitted on.

## Usage

``` r
tda_eigen(
  x,
  algorithm = 1,
  max_iter = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- x:

  a square numeric matrix.

- algorithm:

  1 (default), 2 (`eigen1`) or 3 (`eigen2`).

- max_iter:

  for `algorithm = 3`: iteration limit, TDA's default is 100.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

A data frame with one row per eigenvalue: `re` and `im`, its real and
imaginary parts. `attr(x, "run")` carries the run, whose output holds
the eigenvectors and residuals.

## Details

The eigenvalues come from an exporter rather than from the printed
table, so they arrive at full precision instead of rounded to the print
format.

## Examples

``` r
m <- diag(c(2, 3, 5))
tda_eigen(m)
#>   re im
#> 1  2  0
#> 2  3  0
#> 3  5  0
```
