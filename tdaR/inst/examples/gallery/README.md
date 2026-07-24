# tdaR examples

Each `.R` script here replicates one of TDA's shipped example plots
(from `examples/exam/*.cf`) using only documented, exported `tda_pl_*` /
`tda_ps()` / `tda_graph()` functions -- no internal, undocumented, or
raw-option-string calls. `tools/check_gallery.R` compares the PostScript
TDA writes for each script with TDA's `examples/exam/<stem>.ps`, comment
lines aside; as of session 52, 46 scripts match byte for byte, 5 differ
by design (listed in the tool) and 5 still differ (listed in
handover.md) -- an earlier version of this paragraph claimed every one
had been checked, and plot3/plot4 had their grid arguments swapped,
plot5 the wrong page size and plot10 the wrong marker positions.

Run any of them directly:

    Rscript plot1.R

Each writes its `..._r.png` into `out/` beside the scripts, which is
git-ignored.

## What each one demonstrates

- **plot1.R** -- a function and its derivative
  (`tda_pl_function(deriv=)`)
- **plot2.R** -- a shaded frame and an explicit grid
  (`tda_pl_frame(gray=)`, `tda_pl_grid(at_x=, at_y=)`)
- **plot5.R** -- a second axis at an explicit position
  (`tda_pl_axis(at=, dir=)`)
- **plot6.R** -- a second panel on the same page
  (`tda_pl_panel(origin=)`)
- **plot7.R** -- combining several already-created plots into a grid
  (`tda_combine_ps()`)
- **plot8.R** -- symbols connected by a line (`tda_pl_points()`)
- **plot9.R** -- all nine of TDA's line types, distinguishable side by
  side
- **gd9.R** -- two small graphs, straight and curved edges, arrowheads
  (`tda_pl_graph()`)
- **gd12.R** -- a larger graph with a non-default node shape and several
  curved edges
- **gd15.R** -- a tree layout computed by TDA itself (`tda_graph()` /
  `plot.tda_graph(layout = "tree")`)
- **ple5p.R** -- two survivor curves from a real, hand-computed
  Kaplan-Meier estimate (step functions, with Greenwood's-formula
  confidence bands that widen as fewer subjects remain at risk), each
  shaded TDA's way (`tda_pl_lines(band=)`)
- **scplot1.R** -- a basic scatterplot (`tda_pl_scatter()`)

## Things worth knowing if you write your own

**Draw order matters, always.** PostScript paints strictly in the order
commands are given -- whatever is drawn last sits on top of everything
before it. A shaded band's white "clear" fill (see `band=`) can reach
down to the x axis; if the frame is drawn *before* the band, that white
fill paints right over the axis line. `ple5p.R` draws its frame *after*
the data for exactly this reason.

**`tda_pl_panel()`'s `origin=` is a gap between panels' physical
position, not their visible content.** An axis tick's label extends a
few mm past the plot's logical boundary, so placing a new panel's origin
at exactly the previous panel's width (no gap at all) lets the two
panels' own tick labels collide at the seam. Add a real gap (10-15mm is
usually enough) on top of the width. `plot6.R` shows this.

**`tda_combine_ps()`'s inputs need distinct content, not necessarily
distinct file names** -- it generates its internal file names when
copying inputs together, so two sessions that both used `tda_ps()`'s
default `file="plot.ps"` combine correctly regardless.

## Full file list

All ~53 of TDA's shipped example plots are now replicated here: every
`plot*.cf`, `plot-*.cf`, `gd*.cf`, `ple5p.cf`/`ple3p.cf`, and
`scplot*.cf` that produces a PostScript output. Two exceptions, both
documented in their own script:

- **gdf3.R** does not pixel-match the original: `gdf3.cf` itself fails
  to run to completion through the real TDA binary (confirmed directly),
  so there is no successful reference output to match. The script
  instead demonstrates the same real capability (`tda_gdf()` fitting a
  distribution to right-censored data) on data structured the same way.
- **plot14.R**, **plots.R** substitute a plain-text label for TDA's own
  special-character codes (`@141` etc.), which are not documented
  anywhere in TDA's help text.

Randomly-generated examples (`plot11`, `plot17`, `plot18`, `scplot1-3`,
etc.) use R's RNG, not TDA's `rd()` -- the specific values differ from
the original by construction, but the plotting technique is the same.
Where the original's `rd` (no arguments) is TDA's uniform(0,1)
generator, the R side uses `runif()` to match -- confirmed against TDA's
documented behaviour for `rd`, not assumed.

Examples built from TDA's real, shipped data (`ple5p.R`, `ple3p.R`,
using `tda_rrdat()`/`tda_ple()`) are checked against the real TDA
binary's output, not just visually similar.

## Known, disclosed rendering limitation

`plot10.R`'s markers 6, 7, 14, and 17 (the four symbol types that
combine an outline shape with a cross/x mark drawn through it) render
differently from TDA's real output for this one file. Investigated in
depth: the real cause is that TDA's PostScript generator leaves an
unbalanced `gsave`/`grestore` count by the time it draws these later
symbols in a long, 17-symbol sequence -- checked by counting: the
balance is -2 by marker 17, -5 by marker 14. This is a quirk in TDA's
generated PostScript (confirmed by comparing an isolated test of just
the symbol's procedure calls, which renders cleanly, against the same
sequence run with the real file's accumulated prior state, which does
not). This package's own renderer treats each symbol as an independent
drawing operation and correctly reproduces what each symbol's procedure
describes in isolation; it does not replicate this accumulated,
file-position-dependent PostScript state from TDA's generator. No other
example in this folder uses these four symbol types, so the scope of
this is limited to `plot10.R` alone.
