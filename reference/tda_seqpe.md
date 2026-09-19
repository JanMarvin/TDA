# Derive a sequence data structure from episode data

`seqpe`: builds a sequence – a state at every time point – from episode
data (a start time, an end time, and origin/destination states for each
spell). The result can be read back with a fresh call to
[`tda_seqm`](https://janmarvin.github.io/TDA/reference/tda_seqm.md) or
any of the other sequence functions once written out and re-read; this
function only performs the conversion, since TDA's output is a plain
data file rather than something kept in memory across calls the way
[`tda_seqm`](https://janmarvin.github.io/TDA/reference/tda_seqm.md)
manages internally.

## Usage

``` r
tda_seqpe(
  data,
  id,
  origin,
  destination,
  start,
  end,
  tp,
  missing = NULL,
  missing_start = NULL,
  missing_end = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- data:

  a data frame of episodes, one row per spell.

- id:

  name of the case identifier column in `data`.

- origin, destination:

  names of the origin- and destination-state columns.

- start, end:

  names of the start- and end-time columns.

- tp:

  definition of the time points the sequence is built over, as a vector
  or a single TDA range expression such as `"0(1)10"`.

- missing:

  value to use for a gap at the beginning (`missing_start`), at the end
  (`missing_end`), or anywhere else (`missing`); TDA's default is `-1`
  for all three.

- missing_start, missing_end:

  see `missing`.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

A data frame: the case id, then one column per time point in `tp`, each
holding the state the case was in at that time.

## See also

Other sequence analysis:
[`tda_evalfi()`](https://janmarvin.github.io/TDA/reference/tda_evalfi.md),
[`tda_fml()`](https://janmarvin.github.io/TDA/reference/tda_fml.md),
[`tda_frml()`](https://janmarvin.github.io/TDA/reference/tda_frml.md),
[`tda_seq_info()`](https://janmarvin.github.io/TDA/reference/tda_seq_info.md),
[`tda_seqgc()`](https://janmarvin.github.io/TDA/reference/tda_seqgc.md),
[`tda_seqm()`](https://janmarvin.github.io/TDA/reference/tda_seqm.md),
[`tda_seqmd()`](https://janmarvin.github.io/TDA/reference/tda_seqmd.md),
[`tda_seqpm()`](https://janmarvin.github.io/TDA/reference/tda_seqpm.md)

## Examples

``` r
# two cases, two spells each
e <- data.frame(id = c(1, 1, 2, 2), org = c(1, 2, 1, 3),
                des = c(2, 3, 3, 1), ts = c(0, 3, 0, 2),
                tf = c(3, 5, 2, 5))
tda_seqpe(e, id = "id", origin = "org", destination = "des",
         start = "ts", end = "tf", tp = "0(1)5")
#>   id t0 t1 t2 t3 t4 t5
#> 1  1  1  1  1  2  2 -1
#> 2  2  1  1  3  3  3 -1
```
