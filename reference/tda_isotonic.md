# Isotonic regression, with or without tie groups

The non-decreasing sequence closest (least squares) to `y` – TDA's
`mmp`, and with tie groups its `mmp1`/`mmp2`. Without `groups` this is
plain pool-adjacent-violators and agrees with
[`stats::isoreg`](https://rdrr.io/r/stats/isoreg.html) exactly.

## Usage

``` r
tda_isotonic(
  y,
  groups = NULL,
  ties = c("primary", "secondary"),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- y:

  a numeric vector, in the order to be respected.

- groups:

  optional vector of the same length marking tie groups by consecutive
  runs of equal values.

- ties:

  how tie groups are treated: `"primary"` or `"secondary"`; ignored
  without `groups`.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

The fitted values, a numeric vector like `y`.

## Details

`groups` marks observations whose order among themselves is not
meaningful (ties in the ordering variable). Consecutive equal values
form one group – the grouping follows runs, exactly as TDA reads its tie
matrix, so `c(1, 1, 2)` is two groups but `c(1, 2, 1)` is three. How a
group is treated is `ties`: `"primary"` (Kruskal's primary approach,
`mmp1`) lets the fitted values untie freely – the group's values are
sorted, fitted, and each observation receives its rank slot's value, so
the result is monotone only up to reordering within groups;
`"secondary"` (`mmp2`) forces one common value per group – the group
mean is fitted and broadcast back.

## See also

Other smoothing:
[`tda_integrate()`](https://janmarvin.github.io/TDA/reference/tda_integrate.md),
[`tda_interp()`](https://janmarvin.github.io/TDA/reference/tda_interp.md),
[`tda_mat()`](https://janmarvin.github.io/TDA/reference/tda_mat.md),
[`tda_sma()`](https://janmarvin.github.io/TDA/reference/tda_sma.md),
[`tda_smd()`](https://janmarvin.github.io/TDA/reference/tda_smd.md),
[`tda_spl()`](https://janmarvin.github.io/TDA/reference/tda_spl.md)

## Examples

``` r
tda_isotonic(c(5, 3, 4, 6, 2, 7))
#> [1] 4 4 4 4 4 7
# rows 1-2 tied, rows 4-6 tied:
tda_isotonic(c(5, 3, 4, 6, 2, 7), groups = c(1, 1, 2, 3, 3, 3),
             ties = "secondary")
#> [1] 4 4 4 5 5 5
```
