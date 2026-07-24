# Adversarial fuzz tests

Deliberately weird, extreme, or degenerate `.cf` files -- not TDA's
own shipped examples -- built to stress-test the parser/renderer in
ways ordinary, well-behaved plots never do. Each has a matching
`_mine.R` script using the documented tdaR API, meant to be compared
against the real TDA binary's own output for the `.cf` file of the
same number.

To regenerate a real reference and compare:

    tda cf=fuzzN.cf
    gs -sDEVICE=png16m -r150 -dEPSCrop -o fuzzN_real.png fuzzN.ps
    Rscript fuzzN_mine.R

## What each one targets, and what it found

- **fuzz1** -- a circle far larger than the plot's own declared
  range, a near-zero-scale y axis, long/special-character text.
  Found a real bug: circle ops carried no clip region at all, so an
  oversized circle inflated the whole viewport instead of being
  clipped to the plot's own frame the way TDA's real output is.
  Fixed in `tda_plot_ps()` / `.draw_arc()`.
- **fuzz2** -- two nodes at the identical position, a self-loop, a
  zero-length edge, edges with curvature too extreme for TDA itself
  to draw (silently skipped in the real output). Clean -- matched.
- **fuzz3** -- a strongly anisotropic plot (150mm wide, 15mm tall)
  with overlapping rotated text. Found a real bug: the page-wide
  text-scaling factor compared the largest axis span against the
  smallest available device dimension, mixing two different
  directions' own numbers and shrinking text to an illegible
  fraction of a point. Fixed in `tda_plot_ps()`.
- **fuzz4** -- a near-degenerate, extremely small coordinate range
  (0.0002 units wide). Clean -- matched.
- **fuzz5** -- extreme line widths (5mm and 0.01mm), overlapping
  semi-opaque-looking ellipses, wildly different text sizes. Clean --
  matched.
- **fuzz6** -- four overlapping edges of different curvature between
  the same node pair, plus coincident dashed/dot-dash edges. Clean --
  matched.
