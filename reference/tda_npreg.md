# Non-parametric regression

`tda_npreg` is a non-parametric regression evaluated at points you
choose: TDA's `npreg`, a kernel-smoothed mean, quantile or frequency.

## Usage

``` r
tda_npreg(
  formula,
  data,
  x,
  method = c("mean", "quantile", "frequency", "lowess", "midmeans"),
  kernel = c("uniform", "triangle", "quartic", "epanechnikov"),
  bandwidth = NULL,
  select = NULL,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- formula:

  a two-sided formula.

- data:

  a data frame.

- x:

  for `tda_npreg`, the points at which to evaluate the fit.

- method:

  see Details; `npreg`'s `opt=`. Either a name (partially matched) or
  TDA's number directly (`opt=1` is `"mean"`, and so on) – useful when
  translating a `.cf` file, which only ever writes the number.

- kernel:

  the smoothing kernel, only meaningful for `method = "mean"` –
  `npreg`'s `k=`. A name or TDA's number, the same as `method`.

- bandwidth:

  the kernel bandwidth (`npreg`'s `d=`, default 1).

- select:

  a case-selection expression, TDA's `sel=`.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

An object carrying the fitted values as `table`.

## Details

`method` picks which: `"mean"` (the default), a kernel-smoothed mean at
each point in `x`; `"quantile"`, smoothed quantiles; `"frequency"`,
smoothed frequencies. Only these three exist. `npreg` names the columns
for `method` `"mean"` and `"quantile"`, which `tda_npreg` reproduces;
`"frequency"` writes columns that depend on the data (the distinct
response categories) and are returned as TDA numbered them, `V1`, `V2`,
...

## See also

Other regression: [`TDA_FAMILIES`](tda_glm.md),
[`TDA_QRMODELS`](tda_qreg.md), [`tda_freg()`](tda_freg.md),
[`tda_gdf()`](tda_gdf.md), [`tda_l1reg()`](tda_l1reg.md),
[`tda_lsreg()`](tda_lsreg.md),
[`tda_mlrc_design()`](tda_mlrc_design.md), [`tda_mreg()`](tda_mreg.md),
[`tda_nlreg()`](tda_nlreg.md), [`tda_zreg()`](tda_zreg.md)

## Examples

``` r
set.seed(19)
d <- data.frame(x = round(rnorm(60), 3))
d$y <- round(1 + sin(d$x) + rnorm(60, sd = 0.2), 3)
np <- tda_npreg(y ~ x, d, x = seq(-2, 2, 0.5))
np$table
#>   RECN    X NX          XM        YM       YSD
#> 1    1 -2.0  3 -1.90300000 0.3163333 0.2787406
#> 2    2 -1.5  7 -1.27557143 0.1440000 0.1282433
#> 3    3 -1.0 15 -0.84100000 0.3114000 0.1956454
#> 4    4 -0.5 22 -0.45181818 0.6211364 0.3012182
#> 5    5  0.0 19 -0.07663158 0.9491053 0.3062905
#> 6    6  0.5 20  0.53270000 1.4328500 0.2978702
#> 7    7  1.0 17  0.84864706 1.6395882 0.2886109
#> 8    8  1.5  8  1.53475000 1.9345000 0.2688335
#> 9    9  2.0  5  1.83320000 1.9724000 0.2669191

# a triangle kernel, wider bandwidth
np2 <- tda_npreg(y ~ x, d, x = seq(-2, 2, 0.5), kernel = "triangle",
                 bandwidth = 0.5)
np2$table
#>   RECN    X NX          XM        YM       YSD
#> 1    1 -2.0  1 -1.82700000 0.1890000 0.0000000
#> 2    2 -1.5  1 -1.51300000 0.1240000 0.0000000
#> 3    3 -1.0  7 -1.07957143 0.1782073 0.1528711
#> 4    4 -0.5 14 -0.49385714 0.5344353 0.1960580
#> 5    5  0.0 11 -0.02218182 1.0768303 0.2161748
#> 6    6  0.5  9  0.51211111 1.4221763 0.2066423
#> 7    7  1.0  7  0.90300000 1.6029630 0.2058259
#> 8    8  1.5  5  1.47000000 1.9501246 0.1462970
#> 9    9  2.0  3  1.98566667 2.0862431 0.3072334
```
