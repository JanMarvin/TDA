# Draw spatial data

`sdplot` needs a PostScript coordinate system as well as the spatial
data, so the ranges have to be given.

## Usage

``` r
# S3 method for class 'tda_spatial'
plot(
  x,
  xlim = NULL,
  ylim = NULL,
  width = 100,
  height = NULL,
  newpage = TRUE,
  what = c("objects", "hull", "hull_per_polygon"),
  symbol = NULL,
  size = NULL,
  lty = NULL,
  lw = NULL,
  gray = NULL,
  ...
)
```

## Arguments

- x:

  a [`tda_spatial`](tda_spatial.md).

- xlim, ylim:

  ranges of the coordinate system.

- width, height:

  size of the plotting area in millimetres; the height follows the
  data's aspect ratio unless given.

- newpage:

  start a new page before drawing.

- what:

  what to draw: `"objects"` (default), the convex hull over every point,
  or the convex hull of each polygon separately; `sdplot`'s `opt=`.

- symbol, size, lty, lw, gray:

  marker symbol and size, line type, line width in mm, and grey level
  from 0 (black) to 1 (white); `sdplot`'s own
  `s=`/`fs=`/`lt=`/`lw=`/`gs=`.

- ...:

  further options for `sdplot`.

## Value

The parsed drawing operations, invisibly.

## Details

Points need a marker or nothing is drawn: give `symbol` (or another TDA
symbol number). This is the same trap noted for
[`tda_sd_plot3`](tda_sd_ps3.md) – without a marker or a `zvar`
attribute, `sdplot` draws nothing at all rather than erroring.

## See also

Other spatial analysis: [`tda_map()`](tda_map.md),
[`tda_polygons()`](tda_polygons.md),
[`tda_read_dbf()`](tda_read_dbf.md),
[`tda_read_shapefile()`](tda_read_shapefile.md),
[`tda_read_spatial`](tda_read_spatial.md), [`tda_sd()`](tda_sd.md),
[`tda_sd_analyses`](tda_sd_analyses.md),
[`tda_sd_ps3()`](tda_sd_ps3.md), [`tda_spatial()`](tda_spatial.md)

## Examples

``` r
d <- data.frame(id = 1:4, x = c(0, 1, 1, 0), y = c(0, 0, 1, 1))
s <- tda_spatial(d)
plot(s, symbol = 1)  # a marker is needed; points draw nothing without one
```
