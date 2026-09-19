# Evaluate an inclusion function over an interval

`evalfi`: given an expression and an interval domain for each of its
arguments, returns an interval guaranteed to contain the expression's
true range over that domain – interval arithmetic, not a numeric
evaluation at a point. The bound is conservative rather than tight where
a variable appears more than once (the standard "dependency problem" of
interval arithmetic: each occurrence is widened independently, as if
they were unrelated variables that happened to share a domain) –
`evalfi(x = c(0, 3), "x*x-2*x+1")` returns `[-5, 10]`, not the true
range `[0, 4]` that `(x-1)^2` takes over `[0, 3]`. `evalfi1` is the same
but for the first derivative.

## Usage

``` r
tda_evalfi(
  expr,
  ...,
  derivative = FALSE,
  options = list(),
  dir = tempfile("tda")
)
```

## Arguments

- expr:

  the expression, as a string, in TDA's language.

- ...:

  one interval per argument the expression uses, named for the argument,
  each as `c(lower, upper)` – `evalfi(x = c(0, 3), expr = "...")`.

- derivative:

  if `TRUE`, use `evalfi1` (the first derivative's inclusion function)
  instead of `evalfi`.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

## Value

A named numeric vector, `c(lower, upper)`.

## Details

`^` is not a valid operator in an interval expression – write a power as
repeated multiplication (`x*x`, not `x^2`).

## See also

Other sequence analysis: [`tda_fml()`](tda_fml.md),
[`tda_frml()`](tda_frml.md), [`tda_seq_info()`](tda_seq_info.md),
[`tda_seqgc()`](tda_seqgc.md), [`tda_seqm()`](tda_seqm.md),
[`tda_seqmd()`](tda_seqmd.md), [`tda_seqpe()`](tda_seqpe.md),
[`tda_seqpm()`](tda_seqpm.md)

## Examples

``` r
# true range of (x-1)^2 over [0,3] is [0,4] -- this is wider, since x
# appears twice and interval arithmetic cannot see they are the same x
tda_evalfi("x*x-2*x+1", x = c(0, 3))
#> lower upper 
#>    -5    10 
tda_evalfi("x*x-2*x+1", x = c(0, 3), derivative = TRUE)  # 2x-2, exact
#> lower upper 
#>    -2     4 
                                                         # here: [-2,4]
```
