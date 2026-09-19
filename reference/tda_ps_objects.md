# What is inside a TDA PostScript file

Reopens a PostScript file that TDA wrote (`xopen`) and lists the plot
objects in it (`xlog1`). TDA numbers each thing it draws – an axis, a
series, a convex hull – and records the command that drew it as a
`%\#N: command` comment. This reports that list, which is how you find
the number to pass to a later `xdelete`.

## Usage

``` r
tda_ps_objects(file, options = list(), dir = tempfile("tda"), ...)
```

## Arguments

- file:

  a PostScript file written by TDA.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

A data frame with `object`, the number, and `command`, the command that
drew it – zero rows if the file holds none. `attr(x, "xlim")`,
`attr(x, "ylim")` and `attr(x, "bbox")` carry the coordinate system;
`attr(x, "run")` carries the run.

## Details

The coordinate system and bounding box come back as attributes, because
they are what decide whether a second plot can be added to this file at
all.

## See also

[`tda_read_ps`](tda_read_ps.md), [`tda_combine_ps`](tda_combine_ps.md)

## Examples

``` r
d <- data.frame(X = 1:10, Y = c(2, 4, 3, 6, 5, 8, 7, 10, 9, 12))
p <- tda_ps(d, file = "p.ps", xlim = c(0, 11), ylim = c(0, 13))
tda_ps_objects(tda_ps_file(p))
#> [1] object  command
#> <0 rows> (or 0-length row.names)
```
