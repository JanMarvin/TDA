# Sequence characteristics, lengths, gaps and state distributions

Descriptive summaries of a sequence data structure: `tda_seqgc` basic
characteristics per case (length, number of states visited, number of
transitions, time spent in and number of episodes in each state),
`tda_seqlg` sequence lengths and gaps, `tda_seqsd` the distribution of
states over time (how many cases are in each state at each time point),
`tda_seqsi` a state indicator matrix (one row per case-time combination,
one column per state, 1 where that case is in that state at that time),
and `tda_seqen` entropy of the state distribution at each time point.
Column names come from TDA's own description of each output, not a fixed
layout guessed in advance – `tda_seqgc`'s columns, for instance, depend
on how many distinct states the data actually has.

## Usage

``` r
tda_seqgc(
  sequences,
  select = NULL,
  variables = NULL,
  id = NULL,
  options = list(),
  type = 1,
  dir = tempfile("tda"),
  ...
)

tda_seqlg(
  sequences,
  select = NULL,
  variables = NULL,
  id = NULL,
  options = list(),
  type = 1,
  dir = tempfile("tda"),
  ...
)

tda_seqsd(
  sequences,
  select = NULL,
  variables = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)

tda_seqsi(
  sequences,
  tp,
  select = NULL,
  variables = NULL,
  id = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)

tda_seqen(
  sequences,
  select = NULL,
  variables = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- sequences:

  a data frame or matrix of one column per time point, in order, the
  same shape
  [`tda_seqm`](https://janmarvin.github.io/TDA/reference/tda_seqm.md)
  takes – optionally with further, non-sequence columns (an id, say)
  alongside them; use `variables`/`id` to say which is which rather than
  pre-subsetting `sequences` yourself.

- select:

  optional expression (`sel=`) restricting which sequences are used –
  written with R's comparison operators (`Y0 == 1`, or `Y0 = 1`) or
  TDA's function form (`eq(Y0,1)`) interchangeably; `&`/`|` combine
  several conditions the same way in both.

- variables:

  which columns of `sequences` are the actual time-point columns, in
  order – TDA's `seqdef=Y0,,Y7`, which likewise names the sequence out
  of a data matrix holding other variables too rather than requiring
  them to be dropped first. A character vector of names in `sequences`;
  defaults to every column (so passing an already-subsetted `sequences`
  still works exactly as before).

- id:

  further columns of `sequences` to carry through to the output as-is,
  one value per case – `seqgc`/`seqlg`/ `seqpm`/`seqsi`'s `v=`, which
  writes the case's raw value onto its output row, no more. A character
  vector of names in `sequences`. `seqsd` and `seqen` have no per-case
  row to attach a value to (one row per time point, over all cases) and
  error if given one – TDA accepts `v=` there without complaint but
  never actually writes the column.

- options:

  a named list of further TDA options, passed through.

- type:

  the sequence data type, TDA's `seqdef` `m=`: `1` (default) one column
  per time point, `2` columns read as state/time pairs, one pair per
  spell, with the time axis taken from the times rather than from the
  column count.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

- tp:

  for `tda_seqsi`, the time points to build the indicator matrix at –
  required, as a vector or a single TDA range expression such as
  `"0(1)10"`.

## Value

A data frame, one row per case (`tda_seqgc`, `tda_seqlg`) or per time
point (`tda_seqsd`, `tda_seqen`), or per case-time combination
(`tda_seqsi`).

## See also

Other sequence analysis:
[`tda_evalfi()`](https://janmarvin.github.io/TDA/reference/tda_evalfi.md),
[`tda_fml()`](https://janmarvin.github.io/TDA/reference/tda_fml.md),
[`tda_frml()`](https://janmarvin.github.io/TDA/reference/tda_frml.md),
[`tda_seq_info()`](https://janmarvin.github.io/TDA/reference/tda_seq_info.md),
[`tda_seqm()`](https://janmarvin.github.io/TDA/reference/tda_seqm.md),
[`tda_seqmd()`](https://janmarvin.github.io/TDA/reference/tda_seqmd.md),
[`tda_seqpe()`](https://janmarvin.github.io/TDA/reference/tda_seqpe.md),
[`tda_seqpm()`](https://janmarvin.github.io/TDA/reference/tda_seqpm.md)

## Examples

``` r
s <- data.frame(Y0 = c(1, 1, 2), Y1 = c(1, 2, 1), Y2 = c(2, 1, 1),
                Y3 = c(1, 1, 2), Y4 = c(1, 1, 2), Y5 = c(3, 3, 3))
tda_seqgc(s)   # length, states visited, transitions, time in each state
#>   case slen nds nev dur1 dur2 dur3 durm nep1 nep2 nep3 nepm
#> 1    1    6   3   3    4    1    1    0    2    1    1    0
#> 2    2    6   3   3    4    1    1    0    2    1    1    0
#> 3    3    6   3   3    2    3    1    0    1    2    1    0
tda_seqlg(s)   # lengths and gaps
#>   case ts tf slen glen gmin gmax ngap
#> 1    1  0  5    6    0    0    0    0
#> 2    2  0  5    6    0    0    0    0
#> 3    3  0  5    6    0    0    0    0
tda_seqsd(s)   # how many cases are in each state, at each time point
#>   time nst1 nst2 nst3 valid nmiss total
#> 1    0    2    1    0     3     0     3
#> 2    1    2    1    0     3     0     3
#> 3    2    2    1    0     3     0     3
#> 4    3    2    1    0     3     0     3
#> 5    4    2    1    0     3     0     3
#> 6    5    0    0    3     3     0     3
tda_seqsi(s, tp = "0(1)5")   # state indicator matrix
#>   case s1_0 s2_0 s3_0 s1_1 s2_1 s3_1 s1_2 s2_2 s3_2 s1_3 s2_3 s3_3 s1_4 s2_4
#> 1    1    1    0    0    1    0    0    0    1    0    1    0    0    1    0
#> 2    2    1    0    0    0    1    0    1    0    0    1    0    0    1    0
#> 3    3    0    1    0    1    0    0    1    0    0    0    1    0    0    1
#>   s3_4 s1_5 s2_5 s3_5
#> 1    0    0    0    1
#> 2    0    0    0    1
#> 3    0    0    0    1
tda_seqen(s)   # entropy of the state distribution at each time point
#>   time n       ent
#> 1    0 3 0.6365142
#> 2    1 3 0.6365142
#> 3    2 3 0.6365142
#> 4    3 3 0.6365142
#> 5    4 3 0.6365142
#> 6    5 3 0.0000000

# an id column alongside the sequence, positively selected rather
# than requiring s2[, -1]
s2 <- cbind(id = 101:103, s)
tda_seqgc(s2, variables = names(s), id = "id")
#>   case slen nds nev dur1 dur2 dur3 durm nep1 nep2 nep3 nepm  id
#> 1    1    6   3   3    4    1    1    0    2    1    1    0 101
#> 2    2    6   3   3    4    1    1    0    2    1    1    0 102
#> 3    3    6   3   3    2    3    1    0    1    2    1    0 103
```
