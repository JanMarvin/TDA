# Information about currently defined sequence data structures

Runs TDA's `seq;` against the sequence structure(s) given – the same
table TDA itself prints after defining them, listing each structure's
number of state variables, its time axis, and the distinct states it
actually contains.

## Usage

``` r
tda_seq_info(
  sequences,
  data = NULL,
  options = list(),
  type = 1,
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- sequences:

  a data frame or matrix, one row per case, one column per time point,
  in order – the same shape [`tda_seqm`](tda_seqm.md) takes – or, given
  `data` too, a character vector naming those columns within it
  directly, rather than requiring a separately pre-subsetted data frame:
  real sequence data typically carries an ID column and covariates
  alongside the states in the same data frame (TDA's
  `examples/exam/seq.d4`, e.g.), not only the states on their own. For
  `tda_seqev`/`tda_seqevd`, also a *list* of several such specs, to
  define several independent sequence data structures at once – TDA's
  `seqdef(sn=1)`, `seqdef(sn=2)`, ... – the pattern TDA's own manual
  uses to build a second sequence structure from `seq.d4`'s `S0..S5`
  columns, entirely separate from the first (`Y0..Y5`), not a covariate
  of it. Use `sn=` to say which one a given call analyzes.

- data:

  optional data frame `sequences` (and, for `tda_seqmd`, `covariates`)
  name columns within, instead of `sequences`/`covariates` being the
  actual data already.

- options:

  a named list of further TDA options, passed through.

- type:

  the sequence data type, TDA's `seqdef` `m=`: `1` (default) one column
  per time point, `2` columns read as state/time pairs, one pair per
  spell.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

A data frame, one row per sequence structure: `sn`, `type`, `variables`
(number of time points), `tmin`/`tmax` (the time axis), `nstates`, and
`states` (a list column, the actual state values, since their count
varies by structure).

## See also

Other sequence analysis: [`tda_evalfi()`](tda_evalfi.md),
[`tda_fml()`](tda_fml.md), [`tda_frml()`](tda_frml.md),
[`tda_seqgc()`](tda_seqgc.md), [`tda_seqm()`](tda_seqm.md),
[`tda_seqmd()`](tda_seqmd.md), [`tda_seqpe()`](tda_seqpe.md),
[`tda_seqpm()`](tda_seqpm.md)

## Examples

``` r
s <- data.frame(Y0 = c(1, 1, 2), Y1 = c(1, 2, 1), Y2 = c(2, 1, 1))
tda_seq_info(s)
#>   sn type variables tmin tmax nstates states
#> 1  1    1         3    0    2       2   1, 2

# several structures at once, the same way tda_seqev's sn= does
d <- data.frame(Y0 = c(1, 1, 2), Y1 = c(1, 2, 1), S0 = c(1, 3, 1),
                S1 = c(3, 1, 3))
tda_seq_info(list(c("Y0", "Y1"), c("S0", "S1")), data = d)
#>   sn type variables tmin tmax nstates states
#> 1  1    1         2    0    1       2   1, 2
#> 2  2    1         2    0    1       2   1, 3
```
