# Circle in a 3-d plot

TDA's `plcirc3` on a `psetup3` system.

## Usage

``` r
tda_plcirc3(radius, center = c(0, 0, 0), normal = c(0, 0, 1), ...)
```

## Arguments

- radius:

  circle radius.

- center:

  xyz center.

- normal:

  direction vector of the circle plane.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

path of the PostScript file, invisibly.

## Examples

``` r
# a unit circle in the xy plane, drawn on its own psetup3 system
f <- tda_plcirc3(1, center = c(0, 0, 0), normal = c(0, 0, 1))
tda_plot_ps(tda_read_ps(f))


# the same command added to a session, so it can share a plot with
# other 3-d commands: TDA's own option names, xyz= and dvec=
p <- tda_ps3(xlim = c(-2, 2), ylim = c(-2, 2), zlim = c(-2, 2))
p <- tda_pl(p, "plcirc3", xyz = "0,0,0", dvec = "0,0,1", rhs = 1)
p <- tda_pl(p, "plcirc3", xyz = "0,0,0", dvec = "0,1,0", rhs = 1)
plot(p)
```
