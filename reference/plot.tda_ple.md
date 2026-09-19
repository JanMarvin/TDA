# Plot survivor curves and fitted rates

Base graphics, drawn from the tables the estimators already return, so
no plotting dependency is needed. Censored times are marked, a band is
drawn where a standard error is available, and groups get a legend.

## Usage

``` r
# S3 method for class 'tda_ple'
plot(
  x,
  conf.int = TRUE,
  mark.censored = TRUE,
  col = NULL,
  lty = 1,
  lwd = 1.5,
  xlab = "Time",
  ylab = "Survivor function",
  main = NULL,
  legend = TRUE,
  ylim = c(0, 1),
  ...
)

# S3 method for class 'tda_ltb'
plot(x, ...)

# S3 method for class 'tda_rate'
plot(
  x,
  what = c("survivor", "rate", "density"),
  xlab = "Time",
  ylab = NULL,
  col = "black",
  lwd = 1.5,
  main = NULL,
  ...
)
```

## Arguments

- x:

  an estimate from
  [`tda_km`](https://janmarvin.github.io/TDA/reference/tda_ltb.md),
  `tda_ple`, `tda_ltb`, or a fit from
  [`tda_rate`](https://janmarvin.github.io/TDA/reference/tda_rate.md).

- conf.int:

  draw a band at plus or minus 1.96 standard errors.

- mark.censored:

  mark censoring times with a cross.

- col, lty, lwd:

  colour, line type and width, recycled over groups.

- xlab, ylab, main, ylim:

  as usual.

- legend:

  draw a legend when there is more than one group.

- ...:

  passed to `plot`.

- what:

  for a fitted rate model, which of the survivor, rate or density curves
  to draw. The fit must have been made with `prate`.

## Value

The object, invisibly.

## See also

Other plotting:
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
[`tda_ps3()`](https://janmarvin.github.io/TDA/reference/tda_ps3.md),
[`tda_read_ps()`](https://janmarvin.github.io/TDA/reference/tda_read_ps.md)

## Examples

``` r
set.seed(40)
d <- data.frame(t = round(rexp(60, 0.1), 1) + 0.5, s = rbinom(60, 1, 0.8))
km <- tda_km(Surv(t, s) ~ 1, d)
plot(km)
```
