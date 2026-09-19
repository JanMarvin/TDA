# Sequence pattern matching

Counts, per sequence, how many times each of up to 20 patterns occurs –
TDA's `seqpm`. A pattern is a vector of states and wildcards, one
element per position: `c(3, 3)` looks for state 3 immediately followed
by state 3; `"?"` matches any one state, `"*"` any run of states
(including none), `"+"` any repeat of the state before it, and `"-"` any
run of *identical* states. Several patterns are given as a list,
`list(c(3, 3), c(3, "*", 3))`.

## Usage

``` r
tda_seqpm(
  sequences,
  patterns,
  variables = NULL,
  id = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- sequences:

  a data frame or matrix of one column per time point, in order – the
  same shape
  [`tda_seqm`](https://janmarvin.github.io/TDA/reference/tda_seqm.md)
  takes – optionally with further, non-sequence columns (an id, say)
  alongside them; use `variables`/`id` to say which is which.

- patterns:

  a list of patterns, each a vector of states (numbers) and wildcards
  (`"?"`, `"*"`, `"+"`, `"-"`, as strings), as described above. At most
  20.

- variables:

  which columns of `sequences` are the actual time-point columns, in
  order – TDA's `seqdef=Y0,,Y7`, which likewise names the sequence out
  of a data matrix holding other variables too. A character vector of
  names in `sequences`; defaults to every column (so an
  already-subsetted `sequences` still works unchanged).

- id:

  further columns of `sequences` to carry through to the output as-is,
  one value per case – `seqpm`'s `v=`. A character vector of names in
  `sequences`.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

An object carrying a `table` with one row per case: its sequence id, its
length, one match-count column per pattern, and any `id` columns
requested.

## See also

Other sequence analysis:
[`tda_evalfi()`](https://janmarvin.github.io/TDA/reference/tda_evalfi.md),
[`tda_fml()`](https://janmarvin.github.io/TDA/reference/tda_fml.md),
[`tda_frml()`](https://janmarvin.github.io/TDA/reference/tda_frml.md),
[`tda_seq_info()`](https://janmarvin.github.io/TDA/reference/tda_seq_info.md),
[`tda_seqgc()`](https://janmarvin.github.io/TDA/reference/tda_seqgc.md),
[`tda_seqm()`](https://janmarvin.github.io/TDA/reference/tda_seqm.md),
[`tda_seqmd()`](https://janmarvin.github.io/TDA/reference/tda_seqmd.md),
[`tda_seqpe()`](https://janmarvin.github.io/TDA/reference/tda_seqpe.md)

## Examples

``` r
# six short sequences over three states; does "3,3" (state 3 twice in a
# row) occur, and does the sequence ever return to 1 after leaving it?
s <- data.frame(Y1 = c(1, 1, 3, 3, 3, 1), Y2 = c(1, 3, 3, 3, 1, 2),
                Y3 = c(1, 3, 3, 1, 2, 2), Y4 = c(2, 3, 1, 2, 2, 2),
                Y5 = c(2, 1, 2, 2, 2, 3))
r <- tda_seqpm(s, list(c(3, 3), c(1, "*", 1)))
r$table
#>   case length matches1 matches2
#> 1    1      5        0        1
#> 2    2      5        1        1
#> 3    3      5        1        0
#> 4    4      5        1        0
#> 5    5      5        0        0
#> 6    6      5        0        0

# an id column alongside the sequence, positively selected
s2 <- cbind(person = 101:106, s)
r2 <- tda_seqpm(s2, list(c(3, 3)), variables = names(s), id = "person")
r2$table
#>   case length matches1 Vperson
#> 1    1      5        0     101
#> 2    2      5        1     102
#> 3    3      5        1     103
#> 4    4      5        1     104
#> 5    5      5        0     105
#> 6    6      5        0     106
```
