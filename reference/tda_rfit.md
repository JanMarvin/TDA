# Fit binary relations to a valued graph

TDA's `rfit`: finds the relations with the requested properties closest
to the graph, via linear programming.

## Usage

``` r
tda_rfit(
  d,
  properties = c("symmetric", "transitive"),
  max_solutions = 10,
  variant = c("standard", "experimental"),
  ...
)
```

## Arguments

- d:

  symmetric matrix (proximities).

- properties:

  character subset of "reflexive", "symmetric", "antisymmetric",
  "transitive", "complete".

- max_solutions:

  cap on the number of equally good relations TDA keeps (`max=`, 10).
  When more exist the fit is refused with an error naming this argument;
  a linear-order request (antisymmetric, transitive, complete) on a
  symmetric input has 17 optimal relations for 4 nodes and 491 for 6.
  What comes back is correct but not always complete: against every
  relation on 4 nodes enumerated in R (test-gap-wrappers.R, and 30
  directed cases in the session-52 notes) `rfit`'s value is always the
  optimum and every relation it lists is optimal with the properties,
  but in about one case in five with ties it lists a subset of the tied
  optima (2 of 4, 4 of 5) – a limit of the CACM 449 enumeration it uses.

- variant:

  "standard" runs rfit; "experimental" runs the rfit1 sibling (same
  syntax, alternative implementation).

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

`relations`, one fitted relation matrix per solution, and `output`.

## Examples

``` r
d <- rbind(c(0, 2, 1), c(2, 0, 2), c(1, 2, 0))
r <- tda_rfit(d)
r$relations[[1]]
#>      [,1] [,2] [,3]
#> [1,]    1    1    1
#> [2,]    1    1    1
#> [3,]    1    1    1
```
