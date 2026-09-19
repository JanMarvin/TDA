# Rates in age-period form

TDA's `rap`: builds age-period tables from starting time, cohort, ending
time and destination variables.

## Usage

``` r
tda_rap(cohort, entry_year, exit_year, destination, years, ages, ...)
```

## Arguments

- cohort, entry_year, exit_year, destination:

  integer vectors, one entry per case: birth year, year of entry into
  observation, year of exit, and the destination code (0 censored). The
  engine requires cohort \<= entry_year \<= exit_year.

- years, ages:

  length-2 table ranges.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

`table`, the rate table TDA wrote, and `output`.

## Examples

``` r
# three subjects born 1960/1961, observed 1980-1986: the risk
# table fills along the Lexis diagonals
r <- tda_rap(cohort = c(60, 61, 60), entry_year = c(80, 82, 83),
             exit_year = c(85, 86, 86), destination = c(1, 1, 0),
             years = c(80, 86), ages = c(18, 26))
r$risk
#>       age y80 y81 y82 y83 y84 y85 y86
#>  [1,]  26   0   0   0   0   0   0   1
#>  [2,]  25   0   0   0   0   0   2   1
#>  [3,]  24   0   0   0   0   2   1   0
#>  [4,]  23   0   0   0   2   1   0   0
#>  [5,]  22   0   0   1   1   0   0   0
#>  [6,]  21   0   1   1   0   0   0   0
#>  [7,]  20   1   0   0   0   0   0   0
#>  [8,]  19   0   0   0   0   0   0   0
#>  [9,]  18   0   0   0   0   0   0   0
```
