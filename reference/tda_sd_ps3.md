# Three-dimensional spatial plots

`sdplot31`, `sdplot32` and `sdplot33` draw spatial objects in a
three-dimensional coordinate system rather than the flat one
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) uses. That
system comes from `psetup3`, so the ranges of all three axes and the
direction the plot is viewed from are set up once with `tda_sd_ps3()`,
the 3-D counterpart of
[`tda_ps`](https://janmarvin.github.io/TDA/reference/tda_ps.md) –
`tda_sd_pl3_points()`, `tda_sd_pl3_lines()` and `tda_sd_pl3_polygons()`
then add drawing commands to it the way
[`tda_pl_lines`](https://janmarvin.github.io/TDA/reference/tda_pl.md)
and its siblings add to a `tda_ps` session, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
result.
[`tda_pl_text3`](https://janmarvin.github.io/TDA/reference/tda_ps3.md)
and
[`tda_pl_points3`](https://janmarvin.github.io/TDA/reference/tda_ps3.md)
(TDA's 3-D-aware commands, not tied to a spatial file) work on a
`tda_sd_ps3` session too, which is how a title or an extra marker gets
added – TDA's ordinary `plabel()`
([`tda_pl_labels`](https://janmarvin.github.io/TDA/reference/tda_pl_hist.md))
is 2-D only and refuses outright with “current coordinate system is
3-dimensional” if tried here.
[`tda_pl_lines3`](https://janmarvin.github.io/TDA/reference/tda_ps3.md),
by contrast, does **not** work on a `tda_sd_ps3` session: it folds its
coordinates into a regular `nvar` variable, which needs an ordinary
`nvar` declaration that a spatial session, declared with `sdnvar`
instead, does not have – “Syntax error or undefined variables”. Use
`tda_pl_points3` instead, which takes literal coordinates and needs no
such declaration.

## Usage

``` r
tda_sd_ps3(
  s,
  xlim = NULL,
  ylim = NULL,
  zlim,
  view = c(30, 30),
  width = 100,
  ...
)

tda_sd_pl3(p, cmd, ..., zvar = NULL)

tda_sd_pl3_points(p, symbol = 1, size = NULL, zvar = NULL, zval = NULL, ...)

tda_sd_pl3_lines(p, symbol = NULL, size = NULL, zvar = NULL, zval = NULL, ...)

tda_sd_pl3_polygons(
  p,
  symbol = NULL,
  size = NULL,
  zvar = NULL,
  zval = NULL,
  ...
)

tda_sd_plot3(
  s,
  z,
  xlim = NULL,
  ylim = NULL,
  view = c(30, 30),
  width = 100,
  newpage = TRUE,
  symbol = 1,
  size = NULL,
  zvar = NULL,
  zval = NULL,
  ...
)

tda_sd_plot3_lines(
  s,
  z,
  xlim = NULL,
  ylim = NULL,
  view = c(30, 30),
  width = 100,
  newpage = TRUE,
  symbol = NULL,
  size = NULL,
  zvar = NULL,
  zval = NULL,
  ...
)

tda_sd_plot3_polygons(
  s,
  z,
  xlim = NULL,
  ylim = NULL,
  view = c(30, 30),
  width = 100,
  newpage = TRUE,
  symbol = NULL,
  size = NULL,
  zvar = NULL,
  zval = NULL,
  ...
)
```

## Arguments

- s:

  a
  [`tda_spatial`](https://janmarvin.github.io/TDA/reference/tda_spatial.md)
  structure.

- xlim, ylim:

  ranges of the x and y axes, taken from the structure when not given.

- view:

  direction of the projection, as longitude and latitude in degrees;
  TDA's default is 30,30.

- width:

  horizontal size of the plot in mm.

- ...:

  further options for the command.

- p:

  a `tda_sd_ps3` session.

- cmd:

  the TDA command name.

- zvar:

  an attribute carried in from the source file, used as each point's or
  polygon's height instead of a flat `zval` – for the polygon commands,
  each polygon is raised into a prism of that height rather than drawn
  flat. Name it as it appears in the object record – for a shapefile, a
  column of the `.dbf`.

- symbol:

  marker symbol to draw at each point, as TDA's `s=`. It has its
  argument because `s` is the spatial object.

- size:

  marker size in mm.

- zval:

  a single height for every point or polygon, when they do not have
  their via `zvar`; TDA's default is 0.

- z, zlim:

  the range of the z axis, as two numbers. Unlike x and y there is
  nothing in the spatial file to take it from, so it has no default.

- newpage:

  start a new grid page.

## Value

`tda_sd_ps3()` and the `tda_sd_pl3_*` functions return the session, with
the command added.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a session,
and the one-shot functions, draw the plot and return the parsed drawing
operations invisibly.

## Details

**A north-up view of geographic data** (longitude/latitude, the way an
ordinary map looks) is `view = c(-90, 90)`, not TDA's own default of
`c(30, 30)`. This is a property of TDA's azimuth/elevation convention
itself, not of any particular dataset, so it is the same for any
longitude/latitude data.

**Each of the three drawing commands only draws one kind of object**,
silently drawing nothing at all for any other kind: `tda_sd_pl3_points`
draws points, `tda_sd_pl3_lines` draws lines, `tda_sd_pl3_polygons`
draws polygons. A shapefile of counties, say, is polygons – points on it
draws nothing, not because anything is wrong, but because it is looking
for points that are not there. Check
[`tda_sd_info`](https://janmarvin.github.io/TDA/reference/tda_sd_analyses.md)
first if it is not obvious which kind a structure holds.

`tda_sd_plot3`, `tda_sd_plot3_lines` and `tda_sd_plot3_polygons` are
one-shot equivalents – set up a session, add one drawing command, plot
it – kept for the simple case where nothing else needs adding to the
same plot.

## See also

Other spatial analysis:
[`plot.tda_spatial()`](https://janmarvin.github.io/TDA/reference/plot.tda_spatial.md),
[`tda_map()`](https://janmarvin.github.io/TDA/reference/tda_map.md),
[`tda_polygons()`](https://janmarvin.github.io/TDA/reference/tda_polygons.md),
[`tda_read_dbf()`](https://janmarvin.github.io/TDA/reference/tda_read_dbf.md),
[`tda_read_shapefile()`](https://janmarvin.github.io/TDA/reference/tda_read_shapefile.md),
[`tda_read_spatial`](https://janmarvin.github.io/TDA/reference/tda_read_spatial.md),
[`tda_sd()`](https://janmarvin.github.io/TDA/reference/tda_sd.md),
[`tda_sd_analyses`](https://janmarvin.github.io/TDA/reference/tda_sd_analyses.md),
[`tda_spatial()`](https://janmarvin.github.io/TDA/reference/tda_spatial.md)

## Examples

``` r
set.seed(1)
d <- data.frame(id = 1:8, x = runif(8), y = runif(8), z = runif(8))
s <- tda_spatial(d, attributes = "z")
pdf(NULL)

# composable form: build a session, add commands, plot it -- a title
# (via TDA's 3-D-aware text command, not the 2-D-only plabel) and an
# extra marker layered onto the same plot as the data
p <- tda_sd_ps3(s, zlim = c(0, 1))
p <- tda_sd_pl3_points(p, symbol = 1, size = 3)
p <- tda_pl_text3(p, "eight points", at = c(0.5, 0.5, 0.9))
plot(p)

# the one-shot form, for when nothing else needs adding
tda_sd_plot3(s, z = c(0, 1), symbol = 1, size = 3)

# plot3_lines wants line-type spatial data, not scattered points --
# a single connected line here, four vertices sharing one id
ln <- data.frame(id = 1, x = c(0, 0.3, 0.6, 1), y = c(0, 0.4, 0.5, 1),
                 z = c(0, 0.3, 0.6, 1))
sl <- tda_spatial(ln, type = 2, attributes = "z")
tda_sd_plot3_lines(sl, z = c(0, 1), symbol = 1)

# plot3_polygons draws flat polygons at a single height in the z range,
# or, given zvar, raises each one into a prism of that height -- h here
poly <- data.frame(id = 1, x = c(0, 2, 2, 0), y = c(0, 0, 2, 2),
                   h = c(2, 2, 2, 2))
sp <- tda_spatial(poly, type = 3, attributes = "h")
tda_sd_plot3_polygons(sp, z = c(0, 3), zvar = "h")

# a real shapefile, if one is available: sf's example data, North
# Carolina counties, in a proper north-up view, with city markers and
# labels layered on via TDA's 3-D-aware point/text commands
if (requireNamespace("sf", quietly = TRUE)) {
  nc <- tda_read_shapefile(system.file("shape/nc.shp", package = "sf"))
  cities <- data.frame(
    name = c("Raleigh", "Charlotte"),
    lon = c(-78.6382, -80.8431), lat = c(35.7796, 35.2271))

  pnc <- tda_sd_ps3(nc, zlim = c(0, 1), view = c(-90, 90))
  pnc <- tda_sd_pl3_polygons(pnc)
  # fs= is the marker/font size in mm, independent for each command --
  # the default marker size (2mm) reads as a large blob at this map's
  # own small default width (100mm), so it is turned down here, and
  # the default text size turned up, rather than using both defaults
  pnc <- tda_pl_points3(pnc, cities$lon, cities$lat, rep(0.5, 2),
                        symbol = 5, size = 0.8, lty = 0)
  for (i in seq_len(nrow(cities)))
    pnc <- tda_pl_text3(pnc, cities$name[i],
                        at = c(cities$lon[i] + 0.35, cities$lat[i], 0.5),
                        fs = 4)
  plot(pnc)
}
dev.off()
#> agg_record_1ec1dd1ad22 
#>                      2 
```
