# Combine several PostScript plots into one file

`tda_combine_ps` reaches TDA's `dplot` command, which lays out a grid of
previously-created PostScript files (each its
[`tda_ps()`](https://janmarvin.github.io/TDA/reference/tda_ps.md)
session, already run) into a single output file – one row per element of
`rows`, side by side within a row. Unlike every other `tda_pl_*`
function here, this is not a command added to an in-progress session:
`dplot` combines files that already exist on disk into a new one of its
own, so this takes file paths (or `tda_ps` objects that have already
been run, via
[`tda_ps_file`](https://janmarvin.github.io/TDA/reference/tda_pl.md))
and runs standalone, the way
[`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md) does.

## Usage

``` r
tda_combine_ps(
  rows,
  width = 120,
  height = 80,
  origin = c(10, 10),
  file = "combined.ps",
  ...
)
```

## Arguments

- rows:

  A list of PostScript files to lay out in a grid, one row per list
  element; each element is a character vector of file paths (or `tda_ps`
  run objects) placed side by side in that row.

- width, height:

  Size of the combined output, in mm (TDA's defaults: 120 and 80).

- origin:

  Physical origin of the combined plot, `c(x, y)` in mm (TDA's default:
  `c(10, 10)`).

- file:

  Name of the combined output file.

- ...:

  Further options passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

The result of
[`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md), whose
own run directory holds `file`.

## Examples

``` r
# two small plots, then both side by side in one file
p1 <- tda_ps(xlim = c(0, 10), ylim = c(-1, 1), width = 55, height = 40)
p1 <- tda_pl_axes(tda_pl_function(p1, "sin(x)"))
p2 <- tda_ps(xlim = c(0, 10), ylim = c(0, 100), width = 55, height = 40)
p2 <- tda_pl_axes(tda_pl_function(p2, "x * x"))
r <- tda_combine_ps(list(c(tda_ps_file(p1), tda_ps_file(p2))))
file.exists(file.path(r$dir, "combined.ps"))
#> [1] TRUE
```
