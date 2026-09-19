# Start a new panel within the same session

TDA can lay out more than one plot on a single page by calling `psetup`
again, at a different physical position (`psorg=`) – this is that,
reachable with the same `xlim`/`ylim`/ `width`/`height` vocabulary
[`tda_ps`](https://janmarvin.github.io/TDA/reference/tda_ps.md) itself
uses to start the first one, rather than a raw, hand-built `psetup`
block. Every `tda_pl_*` command after this one draws into the new panel;
the commands already added stay in the first (or previous) panel.

## Usage

``` r
tda_pl_panel(
  p,
  origin,
  width = NULL,
  height = NULL,
  xlim = NULL,
  ylim = NULL,
  ...
)
```

## Arguments

- p:

  a [`tda_ps`](https://janmarvin.github.io/TDA/reference/tda_ps.md)
  session.

- origin:

  the new panel's position **relative to the previous panel's own** (the
  session's `origin=` for the first panel, TDA's default if that was not
  given), `c(dx, dy)` in mm – `c(width, 0)` of the previous panel places
  this one immediately to its right, for instance. Not an absolute page
  position: TDA's `psorg=` is itself relative to whatever the previous
  one was (PostScript's cumulative `translate`), and this mirrors that
  rather than asking the caller for TDA's undocumented default starting
  point.

  `origin` is the gap between each panel's *physical position*, not
  between their visible content: an axis tick's label extends a few mm
  past the plot's logical boundary (the width the tick mark itself sits
  at), so an `origin` exactly equal to the previous panel's `width` (no
  gap at all) can let the two panels' own tick labels collide right at
  the seam. Adding a real gap on top of `width` – 10-15mm is normally
  enough – keeps each panel's labels clear of its neighbour's.

- width, height:

  size of the new panel, in mm (TDA's `pxlen=`/ `pylen=`); default to
  the same size the session itself started with.

- xlim, ylim:

  the new panel's coordinate range; default to the same range the
  session itself started with.

- ...:

  further `psetup` options.

## Value

The session, with the new panel started.

## See also

Other plotting:
[`plot.tda_ple()`](https://janmarvin.github.io/TDA/reference/plot.tda_ple.md),
[`tda_check_ps()`](https://janmarvin.github.io/TDA/reference/tda_check_ps.md),
[`tda_pl()`](https://janmarvin.github.io/TDA/reference/tda_pl.md),
[`tda_pl_arc()`](https://janmarvin.github.io/TDA/reference/tda_pl_arc.md),
[`tda_pl_axis()`](https://janmarvin.github.io/TDA/reference/tda_pl_axis.md),
[`tda_pl_graph()`](https://janmarvin.github.io/TDA/reference/tda_pl_graph.md),
[`tda_pl_hist`](https://janmarvin.github.io/TDA/reference/tda_pl_hist.md),
[`tda_pl_regression()`](https://janmarvin.github.io/TDA/reference/tda_pl_regression.md),
[`tda_pl_scatter()`](https://janmarvin.github.io/TDA/reference/tda_pl_scatter.md),
[`tda_plot_ps()`](https://janmarvin.github.io/TDA/reference/tda_plot_ps.md),
[`tda_ps()`](https://janmarvin.github.io/TDA/reference/tda_ps.md),
[`tda_ps3()`](https://janmarvin.github.io/TDA/reference/tda_ps3.md),
[`tda_read_ps()`](https://janmarvin.github.io/TDA/reference/tda_read_ps.md)

## Examples

``` r
p <- tda_ps(xlim = c(0, 6), ylim = c(-1, 1), width = 50, height = 40)
p <- tda_pl_axes(p, sc = 1)
p <- tda_pl_function(p, "sin(x1)", range = c(0, 6))
p <- tda_pl_text(p, "Plot 1", at = c(4, 0.2))
# 50 (the previous panel's width) plus a real gap, or its own
# tick labels collide with this new panel's own.
p <- tda_pl_panel(p, origin = c(65, 0), width = 50, height = 30)
p <- tda_pl_axes(p, sc = 1)
p <- tda_pl_function(p, "sin(x1)", range = c(0, 6))
p <- tda_pl_text(p, "Plot 2", at = c(4, 0.2))
```
