# Boolean function minimization

Which combinations of binary conditions produce the outcome? TDA's `bfa`
builds the truth table of `y` over the condition variables, finds all
prime implicants (Quine-McCluskey), and selects minimal covers. This is
the machinery behind crisp-set QCA-style analyses.

## Usage

``` r
tda_boolean_min(
  y,
  conditions,
  data = NULL,
  undefined = c("dont_care", "true", "false"),
  algorithm = c("lawler", "petrick", "lawler2", "none"),
  max_implicants = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- y:

  the outcome: a 0/1 vector, or a column name in `data`.

- conditions:

  the condition variables: a 0/1 matrix or data frame, or column names
  in `data` (at most 15).

- data:

  optional data frame supplying the columns.

- undefined:

  how rows never observed are treated: `"dont_care"` (default),
  `"true"`, or `"false"`.

- algorithm:

  cover selection: `"lawler"` (default), `"petrick"`, `"lawler2"`, or
  `"none"` for prime implicants only.

- max_implicants:

  storage cap.

- options:

  a named list of further TDA options for `bfa` (`ptab`, `prot`, ...).

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

A list: `$selections`, a data frame with one row per term of each
minimal cover – `selection` number and one column per condition holding
1, 0, or NA for "does not matter" – and `$expression`, the printed
symbolic forms.

## See also

Other optimization:
[`tda_mlp()`](https://janmarvin.github.io/TDA/reference/tda_mlp.md),
[`tda_mlpi()`](https://janmarvin.github.io/TDA/reference/tda_mlpi.md),
[`tda_mls()`](https://janmarvin.github.io/TDA/reference/tda_mls.md),
[`tda_mqap()`](https://janmarvin.github.io/TDA/reference/tda_mqap.md),
[`tda_mqp()`](https://janmarvin.github.io/TDA/reference/tda_mqp.md)

## Examples

``` r
d <- expand.grid(X1 = 0:1, X2 = 0:1, X3 = 0:1)
d$Y <- as.integer(d$X1 & d$X2 | d$X3)
tda_boolean_min("Y", c("X1", "X2", "X3"), data = d)$selections
#>   selection X1 X2 X3
#> 1         1  1  1 NA
#> 2         1 NA NA  1
```
