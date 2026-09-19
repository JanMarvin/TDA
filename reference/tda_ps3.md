# Three-dimensional plots

`tda_ps3` opens a plot session using `psetup3`, TDA's three-dimensional
coordinate setup, and the rest add to it: points, text, a polyline and a
surface, each taking three coordinates where their two-dimensional
counterparts take two. Note that the 3-D commands are the only ones that
emit `rotate` and `scale`, which
[`tda_read_ps`](https://janmarvin.github.io/TDA/reference/tda_read_ps.md)
does not yet interpret, so the replayed picture may be positioned
wrongly even where the PostScript itself is correct.

## Usage

``` r
tda_ps3(
  data = NULL,
  width = 100,
  xlim = c(0, 100),
  ylim = c(0, 100),
  zlim = c(0, 100),
  file = "plot.ps",
  ...
)

tda_pl_points3(
  p,
  x,
  y,
  z,
  symbol = NULL,
  size = NULL,
  lty = 0,
  lw = NULL,
  gray = NULL,
  arrow = NULL,
  ...
)

tda_pl_text3(
  p,
  label,
  at,
  fs = NULL,
  rotate = NULL,
  center = NULL,
  symbol = NULL,
  white = NULL,
  ...
)

tda_pl_surface3(
  p,
  ru,
  rv,
  f1,
  f2,
  f3,
  contour = NULL,
  lty = NULL,
  lw = NULL,
  gray = NULL,
  ...
)

tda_pl_curve3(p, rx, f1, f2, f3, lty = NULL, lw = NULL, gray = NULL, ...)

tda_pl_lines3(
  p,
  x,
  y,
  z,
  symbol = NULL,
  size = NULL,
  lty = NULL,
  lw = NULL,
  gray = NULL,
  arrow = NULL,
  ...
)
```

## Arguments

- data:

  a data frame, or nothing. Coordinates can also be handed to the
  drawing commands as plain vectors, which are folded into the session's
  data as they arrive.

- width:

  size of the plotting area in millimetres.

- xlim, ylim, zlim:

  ranges of the three axes.

- file:

  name for the PostScript file.

- ...:

  further options for the command.

- p:

  a `tda_ps` session opened by `tda_ps3`.

- x, y, z:

  coordinates, either as numbers or as names of columns in the session's
  data. `plot3` beneath `tda_pl_lines3` names three variables rather
  than taking coordinates, so loose vectors are folded into the
  session's data for it; either form works. At least two points are
  needed in both, since these draw a polyline and one triple is not a
  line.

- symbol, size, lty, lw, gray, arrow:

  for `tda_pl_points3` and `tda_pl_lines3`: marker symbol and size, line
  type, line width, grey level, and an arrowhead at the end, as
  `c(length, width)` in mm – `plotp3`/`plot3`'s `s=`/`fs=`/
  `lt=`/`lw=`/`gs=`/`a=`. `tda_pl_text3` also takes `symbol`, drawn at
  the same position as the text.

- label:

  text to draw.

- at:

  where to draw it, as a triple of coordinates.

- fs:

  font size in mm, for `tda_pl_text3`.

- rotate:

  rotation in degrees, for `tda_pl_text3`.

- center:

  for `tda_pl_text3`, centre the string on `at` rather than starting it
  there; `pltext3`'s `sc=`.

- white:

  for `tda_pl_text3`, draw the label on a white background rather than
  transparently; `pltext3`'s `wf=1`.

- ru, rv:

  for `tda_pl_surface3`, the two parameters' ranges and grid resolution,
  each as `"a,b,n,m"`: from `a` to `b`, evaluated at `n` points, drawn
  with `m` grid lines. Keep `n` modest – 10 or so – if `gray` is also
  given: TDA shades every individual grid cell of the mesh, and the
  PostScript output size (and the time to parse and render it) grows
  with `n` squared, not linearly. Checked: an otherwise ordinary
  `n = 30` surface produced over a million lines of PostScript once
  `gray` was added, and took long enough to parse that it looked hung
  rather than merely slow; `n = 5` renders in about a second.

- f1, f2, f3:

  the x, y and z coordinate as expressions in the curve's or surface's
  parameter(s) (`x` for `tda_pl_curve3`; `u` and `v` for
  `tda_pl_surface3`) – parametric, not a single `z = f(x, y)`: for a
  plain surface height field, use `f1 = "u"`, `f2 = "v"`, and the height
  expression itself for `f3`.

- contour:

  for `tda_pl_surface3`, draw internal grid lines over the surface
  (`cont = 1`) in addition to its outline.

- rx:

  for `tda_pl_curve3`, the range of its one parameter, as TDA's step
  notation `"a(step)b"` – not the comma form `ru`/`rv` take.

## Value

A `tda_ps` session, with the command added.

## See also

Other plotting:
[`plot.tda_ple()`](https://janmarvin.github.io/TDA/reference/plot.tda_ple.md),
[`tda_check_ps()`](https://janmarvin.github.io/TDA/reference/tda_check_ps.md),
[`tda_pl()`](https://janmarvin.github.io/TDA/reference/tda_pl.md),
[`tda_pl_arc()`](https://janmarvin.github.io/TDA/reference/tda_pl_arc.md),
[`tda_pl_axis()`](https://janmarvin.github.io/TDA/reference/tda_pl_axis.md),
[`tda_pl_graph()`](https://janmarvin.github.io/TDA/reference/tda_pl_graph.md),
[`tda_pl_hist`](https://janmarvin.github.io/TDA/reference/tda_pl_hist.md),
[`tda_pl_panel()`](https://janmarvin.github.io/TDA/reference/tda_pl_panel.md),
[`tda_pl_regression()`](https://janmarvin.github.io/TDA/reference/tda_pl_regression.md),
[`tda_pl_scatter()`](https://janmarvin.github.io/TDA/reference/tda_pl_scatter.md),
[`tda_plot_ps()`](https://janmarvin.github.io/TDA/reference/tda_plot_ps.md),
[`tda_ps()`](https://janmarvin.github.io/TDA/reference/tda_ps.md),
[`tda_read_ps()`](https://janmarvin.github.io/TDA/reference/tda_read_ps.md)

## Examples

``` r
x <- c(1, 3, 5, 7)
y <- c(2, 6, 3, 8)
z <- c(1, 4, 2, 6)
# no data frame needed: the vectors are folded in as they are used
p <- tda_ps3(xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
p <- tda_pl_points3(p, x, y, z, symbol = 4)
p <- tda_pl_lines3(p, x, y, z)
p <- tda_pl_text3(p, "a corner", at = c(1, 2, 1))

# column names work too, resolved against the session's data
d <- data.frame(x = x, y = y, z = z)
p2 <- tda_ps3(d, xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
p2 <- tda_pl_points3(p2, "x", "y", "z", symbol = 4)

# a surface, defined parametrically over a u,v grid rather than from
# data -- plsurf3 needs all three coordinate functions and both
# ranges explicitly, not a single z = f(x,y) expression
p3 <- tda_ps3(xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
p3 <- tda_pl_surface3(p3, ru = "0,8,30,10", rv = "0,10,30,10",
                      f1 = "u", f2 = "v", f3 = "u*v/8")

# shaded, with internal contour lines -- n kept small (5, not 30) since
# gray= shades every grid cell individually; see the ru/rv argument
p3b <- tda_ps3(xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
p3b <- tda_pl_surface3(p3b, ru = "0,8,5,5", rv = "0,10,5,5",
                       f1 = "u", f2 = "v", f3 = "u*v/8",
                       gray = 0.9, contour = TRUE)

# a parametric curve -- rx takes TDA's step notation (a(step)b),
# not the comma form the surface's ru/rv take
p4 <- tda_ps3(xlim = c(-2, 2), ylim = c(-2, 2), zlim = c(-2, 2))
p4 <- tda_pl_curve3(p4, rx = "0(0.1)6.28",
                    f1 = "cos(x)", f2 = "sin(x)", f3 = "x/3-1",
                    lty = "dotted", lw = 0.6)

# a single point: plotp3 (what tda_pl_points3 calls) rejects one triple
# outright, so this duplicates it -- the same marker drawn twice at the
# same place, joined by an invisible zero-length line, looks identical
# to drawing it once
p4 <- tda_ps3(xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
p4 <- tda_pl_points3(p4, 4, 5, 3, symbol = 5, size = 3)

# a rotated, centred label on a white background, next to a bounding
# line -- a lone label with nothing else in the plot has no coordinate
# range to size the plot from and would come out unreadably small
p5 <- tda_ps3(xlim = c(0, 8), ylim = c(0, 10), zlim = c(0, 8))
p5 <- tda_pl_points3(p5, c(0, 8), c(0, 10), c(0, 8), symbol = 1,
                     size = 0.1)
p5 <- tda_pl_text3(p5, "hello", at = c(4, 5, 4), fs = 4, rotate = 30,
                   center = TRUE, white = TRUE)
```
