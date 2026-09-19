# Read TDA's output files

`tda_read_table` reads one table, taking column names from the comment
header where they can be matched to the columns. `tda_read_blocks`
splits a file holding several tables, which is what a life table writes
when there are groups. `tda_file` and `tda_blocks` do the same relative
to a run's directory, and report what TDA said if the file is missing.

## Usage

``` r
tda_read_table(file, col_names = NULL)

tda_read_blocks(file, col_names = NULL)

tda_blocks(x, file, ...)

tda_file(x, file, ...)
```

## Arguments

- file:

  a file name.

- col_names:

  column names to apply; a list for `tda_read_blocks`, one entry per
  block.

- x:

  a run or fitted model.

- ...:

  passed on.

## Value

A data frame, or a list of them.

## See also

Other TDA infrastructure:
[`tda_ccnt()`](https://janmarvin.github.io/TDA/reference/tda_ccnt.md),
[`tda_help()`](https://janmarvin.github.io/TDA/reference/tda_help.md),
[`tda_output()`](https://janmarvin.github.io/TDA/reference/tda_estimates.md),
[`tda_run()`](https://janmarvin.github.io/TDA/reference/tda_run.md),
[`tda_write_data()`](https://janmarvin.github.io/TDA/reference/tda_write_data.md)

## Examples

``` r
set.seed(39)
d <- data.frame(income = round(rlnorm(30, 8, 0.5)))
fit <- tda_ineq(d)
# tda_file() reads a named output file relative to the run's directory
tda_file(fit, "out.txt", col_names = c("index", "cases", "min", "max",
                                       "mean", "sd", "vcoeff", "gini"))
#>   index cases  min  max     mean      sd    vcoeff      gini
#> 1     1    30 1256 5459 2929.233 1099.43 0.3753304 0.2099765

# tda_blocks(): a life table with groups writes several tables into one
# file -- counts, then a survivor/density/rate table, once per group
set.seed(1)
d2 <- data.frame(t = c(4, 3, 1, 5, 8, 2), s = c(1, 1, 0, 1, 1, 1),
                 g = c(1, 1, 1, 2, 2, 2))
lt <- tda_ltb(Surv(t, s) ~ as.factor(g), d2, tp = seq(0, 10, 2))
length(tda_blocks(lt, "out.ltb"))  # 4: two groups, two tables each
#> [1] 4

# tda_read_blocks(): the same thing from a raw file path, when there is
# no fit object to resolve one from
path <- file.path(lt$run$dir, "out.ltb")
length(tda_read_blocks(path))
#> [1] 4
```
