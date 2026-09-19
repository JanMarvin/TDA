# Sequence analysis by optimal matching

Computes proximities between sequences and returns a `dist`, so
`hclust`, `cmdscale` and the rest work on the result directly.

## Usage

``` r
tda_seqm(
  sequences,
  method = 1,
  indel = NULL,
  subcost = NULL,
  common_length = FALSE,
  max_restrict = NULL,
  compare_with = NULL,
  preprocess = NULL,
  n_random = NULL,
  print = NULL,
  dp_matrix = FALSE,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- sequences:

  one row per case, one column per time point – a data frame or matrix,
  or a list of vectors of differing lengths, one per case, padded with
  `NA` to a common width internally (TDA's `seqdef` accepts `NA` for a
  case's unused trailing positions). States may be numeric, character or
  factor; the alphabet is taken over the whole frame.

- method:

  optimal matching method, TDA's `m=`.

- indel:

  indel cost, TDA's `icost=`.

- subcost:

  substitution costs, TDA's `scost=`: either a single number selecting
  one of TDA's built-in cost schemes, or a full `k x k` matrix of costs
  (`k` = number of distinct states in `sequences`), which gets written
  out as an `mdef()` block and referenced by name (TDA itself only
  accepts a matrix *name* for `scost=`, never a literal list of
  numbers). Rows and columns follow the states in ascending order, TDA's
  internal ordering – see `states` in the return value for exactly what
  that order is for a given call.

- common_length:

  require all sequences to share one length before aligning them
  (`seqm`'s `rr=1`); TDA's manual text for this is brief ("use common
  sequence length"), so this is documented no further than that.

- max_restrict:

  restrict the alignment search (`seqm`'s `max=`); TDA's manual text for
  this is brief ("alignment restriction"), so this is documented no
  further than that.

- compare_with:

  compare every sequence against specific sequences rather than every
  pair (`seqm`'s `cn=`, one or more case numbers); TDA's manual text for
  this is brief ("compare with specified sequences"), so this is
  documented no further than that.

- preprocess:

  one or both of `"skip_gaps"` (ignore internal gaps – missing values –
  during alignment) and `"skip_identical"` (collapse runs of identical
  states before aligning) – `seqm`'s `sm=`.

- n_random:

  request an approximate random sample of this many sequence pairs
  rather than comparing every pair – `seqm`'s own `r=`.

- print:

  add further columns to `attr(x, "pairs")` beyond the default
  index/length/distance: `"sequential"` replaces the single `distance`
  column with a running, cumulative-so-far cost at each time point
  (`dist_t1..dist_tn`, `n` = number of time points in `sequences`), also
  built into `attr(x, "dist_by_time")`, an `n x n x n_time` array –
  `[, , t]` is the full pairwise distance matrix using only the
  alignment cost accumulated up to time point `t`, the same symmetric,
  zero-diagonal shape `x` itself has. `distance` itself becomes the
  final time point, `dist_tn`, kept for compatibility with the default
  shape. `"lcs"` adds, alongside the usual single `distance`, the
  longest common subsequence's own length (`lcs_length`) and its state
  sequence, one column per time point (`lcs_t1..lcs_tn`, `-1`-padded
  past the LCS's length). `"sequential"` errors if `compare_with` is
  also given, its shape there not being one this wrapper parses; `"lcs"`
  with `compare_with` runs but leaves its extra columns unnamed for the
  same reason (one distance/LCS group per reference sequence, not the
  single group above) – read `attr(x, "pairs")` (or `table` for the
  `compare_with` shape) directly for these.

- dp_matrix:

  capture TDA's step-by-step dump of the alignment algorithm's internal
  dynamic-programming matrix (`tst=2,3` + `df=`), one block per pair
  actually compared – an illustration of how the distance was arrived
  at, not part of TDA's regular output at all. When `TRUE`, `dp_matrix`
  in the return value is a list keyed `"caseA_caseB"`, each holding
  `seqA`/`seqB` (the two sequences, states in order) and a `D` matrix
  (plus `E`/`F` for `method`s with affine costs) – `D[i, j]` is the
  alignment cost using only the first `i` states of `seqA` and first `j`
  of `seqB`, so `D`'s bottom-right corner is that pair's final distance.
  Off by default: turning it on can produce a large file, one block per
  pair compared.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

Ordinarily a `dist`, with the pairwise table kept in its `pairs`
attribute, and, with `print = "sequential"`, the per-time-point
breakdown also in `dist_by_time` (see `print`). Either shape also
carries `states`, the distinct states `sequences` actually has, in the
ascending order `subcost` (if given as a matrix) must follow, and, with
`dp_matrix = TRUE`, the per-pair alignment breakdown in `dp_matrix` (see
`dp_matrix`). Also carries `run`, the underlying
[`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md) result
– `run$dir` is where the generated TDA commands (`$run$commands`) and,
with `dp_matrix = TRUE`, the raw, unparsed debug file TDA itself wrote
(`dp.tst`) both sit, if what `dp_matrix` parsed out of it needs checking
directly. With `compare_with` given, there is no full pairwise matrix to
return – `cn=` writes each sequence's distance to the reference set
instead – so the return value is instead an object carrying that as
`table` (`case`, `length`, `distance`), and `states` the same way.

## Costs

With the default costs a substitution costs twice an indel, so it never
beats a delete plus an insert, and the result is the indel distance
`|a| + |b| - 2 * LCS(a, b)`. Setting `indel = 2` makes substitution the
cheaper move and the result becomes twice the Levenshtein distance.

## See also

Other sequence analysis:
[`tda_evalfi()`](https://janmarvin.github.io/TDA/reference/tda_evalfi.md),
[`tda_fml()`](https://janmarvin.github.io/TDA/reference/tda_fml.md),
[`tda_frml()`](https://janmarvin.github.io/TDA/reference/tda_frml.md),
[`tda_seq_info()`](https://janmarvin.github.io/TDA/reference/tda_seq_info.md),
[`tda_seqgc()`](https://janmarvin.github.io/TDA/reference/tda_seqgc.md),
[`tda_seqmd()`](https://janmarvin.github.io/TDA/reference/tda_seqmd.md),
[`tda_seqpe()`](https://janmarvin.github.io/TDA/reference/tda_seqpe.md),
[`tda_seqpm()`](https://janmarvin.github.io/TDA/reference/tda_seqpm.md)

## Examples

``` r
s <- data.frame(t1 = c(1, 1, 2), t2 = c(1, 2, 2), t3 = c(2, 2, 3))
d <- tda_seqm(s)
hclust(d)
#> 
#> Call:
#> hclust(d = d)
#> 
#> Cluster method   : complete 
#> Number of objects: 3 
#> 
```
