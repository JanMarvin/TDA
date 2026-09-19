# Run a graph command

`tda_g` reaches any of TDA's forty graph commands; the named functions
are the ones with an obvious meaning.

## Usage

``` r
tda_g(g, cmd, options = list(), dir = tempfile("tda"), ...)
```

## Arguments

- g:

  a
  [`tda_graph`](https://janmarvin.github.io/TDA/reference/tda_graph.md).

- cmd:

  the command name, e.g. `"gsp"`, `"gmst"`.

- options:

  a named list of options for the command.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

An object carrying a `table` and the run.

## See also

Other graph analysis:
[`plot.tda_graph()`](https://janmarvin.github.io/TDA/reference/plot.tda_graph.md),
[`tda_dmet()`](https://janmarvin.github.io/TDA/reference/tda_dmet.md),
[`tda_g_analyses`](https://janmarvin.github.io/TDA/reference/tda_g_analyses.md),
[`tda_graph()`](https://janmarvin.github.io/TDA/reference/tda_graph.md),
[`tda_locate_line()`](https://janmarvin.github.io/TDA/reference/tda_locate_line.md),
[`tda_ptree()`](https://janmarvin.github.io/TDA/reference/tda_ptree.md)

## Examples

``` r
e <- data.frame(from = c(1, 1, 2, 3), to = c(2, 3, 4, 4), value = 1)
g <- tda_graph(e, directed = FALSE)
tda_g(g, "gmst")$table  # minimum spanning tree
#>   tree node_i node_j value
#> 1    1      1      2     1
#> 2    1      2      4     1
#> 3    1      3      4     1
```
