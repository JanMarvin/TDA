# Parametric surface in a 3-d plot

TDA's `plsurf3`: x, y, z as expressions in u and v.

## Usage

``` r
tda_plsurf3(
  fx,
  fy,
  fz,
  u_range = c(-1, 1, 5, 10),
  v_range = c(-1, 1, 5, 10),
  ...
)
```

## Arguments

- fx, fy, fz:

  coordinate expressions in `u` and `v` (TDA syntax).

- u_range, v_range:

  c(from, to, grid_lines, points_per_line).

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

path of the PostScript file, invisibly.

## Examples

``` r
# z = sin(u) cos(v) over a 5 x 5 grid of lines, 12 points each
f <- tda_plsurf3(fx = "u", fy = "v", fz = "sin(u) * cos(v)",
                 u_range = c(-2, 2, 5, 12),
                 v_range = c(-2, 2, 5, 12))
tda_plot_ps(tda_read_ps(f))

# tda_pl_surface3() adds the same surface to a tda_ps3() session
```
