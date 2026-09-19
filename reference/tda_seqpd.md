# Write sequence data back out

`seqpd` prints the sequence data in one of four layouts, chosen by
`layout` (TDA's `m=`):

- 1:

  every sequence structure side by side, one row per case

- 2:

  one row per case and structure

- 3:

  one row per case, structure and position

- 4:

  one row per spell: its origin and destination state, and the times it
  spans

## Usage

``` r
tda_seqpd(
  sequences,
  layout = 1,
  select = NULL,
  variables = NULL,
  id = NULL,
  options = list(),
  type = 1,
  second = NULL,
  second_type = 2,
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

- layout:

  1 to 4, the shape of the result, as listed above.

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

- second:

  names of the columns of a second sequence structure, declared beside
  the first as the manual's `seq5.cf` does; the tables then carry one
  column per structure.

- second_type:

  the sequence data type of the second structure, as `type`; default 2,
  state/time pairs.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

A data frame in the requested layout.

## Examples

``` r
# four cases over three states, one sequence structure
s <- data.frame(ID = 1:4,
                Y1 = c(1, 2, 3, 1), Y2 = c(1, 2, 3, 3),
                Y3 = c(2, 2, 1, 3), Y4 = c(2, 3, 1, 1))
tda_seqpd(s, id = "ID")
#>   case y0_0 y0_1 y0_2 y0_3 id
#> 1    1    1    1    2    2  1
#> 2    2    2    2    2    3  2
#> 3    3    3    3    1    1  3
#> 4    4    1    3    3    1  4

# one row per spell instead, with the times each spans
tda_seqpd(s, layout = 4, id = "ID")
#>   case sn org des ts tf id
#> 1    1  1   1   2  0  2  1
#> 2    1  2   2   2  2  4  1
#> 3    2  1   2   3  0  3  2
#> 4    2  2   3   3  3  4  2
#> 5    3  1   3   1  0  2  3
#> 6    3  2   1   1  2  4  3
#> 7    4  1   1   3  0  1  4
#> 8    4  2   3   1  1  3  4
#> 9    4  3   1   1  3  4  4
```
