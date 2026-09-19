# One-dimensional multifacility location

Place k new (variable) points among m existing (fixed) points on a line
so that the weighted interaction costs are minimal – TDA's `gloc`. Each
variable point ends up at one of the fixed points; only the order of the
fixed points on the line matters, not their spacing (verified against
exhaustive search in the test suite).

## Usage

``` r
tda_locate_line(
  fixed_weights,
  variable_weights = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- fixed_weights:

  a k x m matrix: interaction weight between each variable point (row)
  and each fixed point (column, in line order).

- variable_weights:

  optional symmetric k x k matrix of weights between the variable points
  themselves (default: none).

- options:

  a named list of further TDA options for `gloc`.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md). The graph-setup options `opt=`,
  `gt=`, `perm=` and `sc=` pass through to the underlying edge load, as
  in [`tda_g`](tda_g.md).

## Value

An integer vector of length k: for each variable point, the index of the
fixed point it is placed at.

## See also

Other graph analysis: [`plot.tda_graph()`](plot.tda_graph.md),
[`tda_dmet()`](tda_dmet.md), [`tda_g()`](tda_g.md),
[`tda_g_analyses`](tda_g_analyses.md), [`tda_graph()`](tda_graph.md),
[`tda_ptree()`](tda_ptree.md)

## Examples

``` r
# two facilities, four sites; facility 1 mostly serves site 2,
# facility 2 site 3, and the two facilities interact
vf <- rbind(c(8, 1, 1, 2), c(1, 1, 7, 3))
tda_locate_line(vf, variable_weights = rbind(c(0, 5), c(5, 0)))
#> [1] 2 3
```
