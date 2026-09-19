# Build polygons from line segments

`sdcpol` takes a set of line segments, finds where they meet and
assembles the polygons they enclose. It works on a plain table of
segments rather than on a
[`tda_spatial`](https://janmarvin.github.io/TDA/reference/tda_spatial.md)
structure, and with `option = 3` it writes a spatial data file the rest
of the family can read.

## Usage

``` r
tda_polygons(segments, option = 3, dir = tempfile("tda"), ...)
```

## Arguments

- segments:

  a data frame of four columns: the two ends of each segment, as x1, y1,
  x2, y2.

- option:

  1 prints the node list of the graph the segments form, 2 its edge
  list, 3 constructs the polygons.

- dir:

  working directory for the run.

- ...:

  further options for `sdcpol`; `tol` is the tolerance for treating two
  points as the same, default 1e-4.

## Value

For `option = 3` a `tda_spatial` structure; otherwise the run, with the
table in `$table`.

## See also

Other spatial analysis:
[`plot.tda_spatial()`](https://janmarvin.github.io/TDA/reference/plot.tda_spatial.md),
[`tda_map()`](https://janmarvin.github.io/TDA/reference/tda_map.md),
[`tda_read_dbf()`](https://janmarvin.github.io/TDA/reference/tda_read_dbf.md),
[`tda_read_shapefile()`](https://janmarvin.github.io/TDA/reference/tda_read_shapefile.md),
[`tda_read_spatial`](https://janmarvin.github.io/TDA/reference/tda_read_spatial.md),
[`tda_sd()`](https://janmarvin.github.io/TDA/reference/tda_sd.md),
[`tda_sd_analyses`](https://janmarvin.github.io/TDA/reference/tda_sd_analyses.md),
[`tda_sd_ps3()`](https://janmarvin.github.io/TDA/reference/tda_sd_ps3.md),
[`tda_spatial()`](https://janmarvin.github.io/TDA/reference/tda_spatial.md)

## Examples

``` r
sq <- data.frame(x1 = c(0, 1, 1, 0), y1 = c(0, 0, 1, 1),
                 x2 = c(1, 1, 0, 0), y2 = c(0, 1, 1, 0))
tda_polygons(sq)
#> TDA spatial data (poly.sd)
#>   (no summary in the reader's output; see $run$output)
```
