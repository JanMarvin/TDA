# The episode data behind a model

`tda_episodes` returns the episode data `edef` builds from a formula,
which is what every model in this package is actually fitted to. It is
the way to check that a `Surv()` specification means what you intended
before trusting a fit.

## Usage

``` r
tda_episodes(
  formula,
  data,
  vars = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)

tda_state_dist(
  formula,
  data,
  times,
  id = NULL,
  spell = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- formula:

  the same `Surv()` formula the model functions take.

- data:

  a data frame.

- vars:

  additional variables to carry into the output, as names.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

- times:

  time points at which to evaluate the state distribution.

- id, spell:

  for `tda_state_dist`, columns identifying the case and numbering its
  spells, for multi-episode data (`edef`'s `id=` and `sn=`).

## Value

A data frame of episodes, or of state proportions by time.

## Details

`tda_state_dist` gives the distribution over states at a set of time
points – how many cases occupy each state when – which is the natural
summary of multi-state data.

## See also

Other episodes: [`tda_esort()`](tda_esort.md)

## Examples

``` r
d <- data.frame(t = c(4, 3, 1, 5), s = c(1, 0, 1, 1), x = c(2, 1, 3, 1))
tda_episodes(Surv(t, s) ~ x, d)
#>   episode case subsample transition org des ts tf x
#> 1       1    1         1          1   0   1  0  4 2
#> 2       2    1         1          1   0   0  0  3 1
#> 3       3    1         1          1   0   1  0  1 3
#> 4       4    1         1          1   0   1  0  5 1
tda_state_dist(Surv(t, s) ~ x, d, times = c(1, 3, 5))
#>   time state1 total missing
#> 1    1      3     3       1
#> 2    3      2     2       2
#> 3    5      0     0       4
```
