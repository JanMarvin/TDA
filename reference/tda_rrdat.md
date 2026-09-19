# The job-episode data from Blossfeld and Rohwer

`rrdat.1` is the data file used by nearly every example in
`examples/ehhnew`, the command files accompanying Blossfeld and Rohwer,
*Techniques of Event History Modeling*. This reads it and names its
columns as the examples do.

## Usage

``` r
tda_rrdat(file = NULL, states = c(2, 4), upward = c("gt", "ge"))
```

## Arguments

- file:

  path to `rrdat.1`. Defaults to the copy shipped with the package.

- states:

  2 for a single destination state, 4 to distinguish upward, lateral and
  downward job changes.

- upward:

  with `states = 4`, how an upward move is coded: `"gt"` is `ed2.cf`'s
  `gt(PRESN/PRES - 1, 0.2)`, `"ge"` is `rt1m.cf`'s `ge(...)`; one
  episode lies on the boundary, so the two give 84/219 and 85/218
  upward/lateral moves. Both use TDA's comparison tolerance.

## Value

A data frame of 600 episodes.

## Source

The data are 600 job episodes of 201 randomly selected respondents from
the German Life History Study (GLHS), collected at the
Max-Planck-Institut für Bildungsforschung, Berlin, in 1981-1983 for the
birth cohorts 1929-31, 1939-41 and 1949-51 (Mayer and Brückner 1989).
TDA's manual (section 3.3.3) thanks Karl Ulrich Mayer and Hans-Peter
Blossfeld for providing the data set; it has shipped with TDA since 1994
as `rrdat.1` and is the example throughout the manual and the book.

## References

Mayer, K. U., Brückner, E. (1989). *Lebensverläufe und
Wohlfahrtsentwicklung. Konzeption, Design und Methodik der Erhebung von
Lebensverläufen der Geburtsjahrgänge 1929-1931, 1939-1941, 1949-1951.*
Berlin: Max-Planck-Institut für Bildungsforschung.

Blossfeld, H.-P., Rohwer, G. (2002). *Techniques of Event History
Modeling. New Approaches to Causal Analysis*, 2nd ed. Mahwah, NJ:
Lawrence Erlbaum.

The derived variables are the ones the examples define in their `nvar`
blocks:

- TFP:

  duration of the episode, `TFin - TStart + 1`.

- DES:

  destination state. With `states = 2` it is 1 for a job change and 0
  for a censored episode, which is what most of the examples use. With
  `states = 4` the job changes are split by prestige: 1 upward, 2
  lateral, 3 downward.

- LFX:

  labour force experience before the episode, `TStart - TE`.

- PNOJ:

  number of previous jobs, `NOJ - 1`.

- COHO2, COHO3:

  birth cohort indicators, 1939-41 and 1949-51.

## Examples

``` r
d <- tda_rrdat()
tda_ltb(Surv(TFP, DES) ~ 1, d, tp = seq(0, 500, by = 30))
#> Call: tda_ltb(formula = Surv(TFP, DES) ~ 1, data = d, tp = seq(0, 500, 
#>     by = 30))
#> 
#> Episodes: 600 
#> 
#>  start midpoint entering censored exposed events       prob
#>      0       15      600       28   586.0    223 0.38054608
#>     30       45      349       23   337.5    113 0.33481481
#>     60       75      213       15   205.5     51 0.24817518
#>     90      105      147       16   139.0     25 0.17985612
#>    120      135      106       15    98.5     24 0.24365482
#>    150      165       67        5    64.5      9 0.13953488
#>    180      195       53        9    48.5      4 0.08247423
#>    210      225       40        5    37.5      3 0.08000000
#>    240      255       32        5    29.5      0 0.00000000
#>    270      285       27        7    23.5      2 0.08510638
#> ... 5 more rows
```
