# Estimates and raw output

`tda_estimates` returns the coefficient table, reading the column layout
from TDA's header so that it works for the rate models (SN/Org/Des/MT),
`qreg` (Cat/Term) and `fml` (Parameter/Value) alike. A parameter TDA
held fixed, printed as `---`, comes back `NA`. `tda_output` prints
everything TDA wrote.

## Usage

``` r
tda_output(x)

tda_estimates(x)
```

## Arguments

- x:

  a run or fitted model.

## Value

A data frame, a list of data frames, or `NULL`.

## See also

Other TDA infrastructure:
[`tda_ccnt()`](https://janmarvin.github.io/TDA/reference/tda_ccnt.md),
[`tda_help()`](https://janmarvin.github.io/TDA/reference/tda_help.md),
[`tda_read_table()`](https://janmarvin.github.io/TDA/reference/tda_read_table.md),
[`tda_run()`](https://janmarvin.github.io/TDA/reference/tda_run.md),
[`tda_write_data()`](https://janmarvin.github.io/TDA/reference/tda_write_data.md)

## Examples

``` r
set.seed(38)
d <- data.frame(x = rnorm(100))
d$y <- rbinom(100, 1, plogis(0.2 + 0.7 * d$x))
f <- tda_qreg(y ~ x, d)
tda_estimates(f)
#>   Idx Cat Term  Variable     Coeff     Error  C/Error    Signif
#> 1   1   1    I Intercept 0.4058624 0.2195944 1.848236 0.9354318
#> 2   2   1    X         x 0.9012704 0.2672257 3.372694 0.9992556
tda_output(f)  # everything TDA printed, not just the coefficient table
#> TDA. Analysis of Transition Data (6.4q). Sun Sep 20 03:44:26 2026
#> Current memory: 390032 bytes.
#> 
#> Reading command file: commands
#> ============================================================================
#> Idx Variable  T   S  PFmt  Definition
#> -------------------------------------
#>   1 Vy        3   8  24.16 c1
#>   2 Vx        3   8  24.16 c2
#> 
#> Reading a data frame to create internal data matrix.
#> Maximum number of cases: 100
#> Created a data matrix with 2 variables and 100 cases.
#> ----------------------------------------------------------------------------
#> qreg(...)=...
#> Quantal response models. Current memory: 391668 bytes.
#> 
#> Model: binary logit.
#> 
#> Variables (cross-section)
#> -------------------------
#> Y   : Vy  
#> X1  : Vx  
#> 
#> Checking available data (pmin=1)
#> Number of cases with valid data: 100
#> 
#> Categories of dependent variable.
#> Maximum number of categories: 100
#> 
#> Index               0         1   (Weighted)
#> Category            0         1  Observations
#> ---------------------------------------------
#> Wave 1   N      41.00     59.00        100.00
#>          Pct    41.00     59.00  
#> 
#> Maximum likelihood estimation.
#> Algorithm 5: Newton (I)
#> 
#> Number of model parameters: 2
#> Type of covariance matrix: 2
#> Maximum number of iterations: 20
#> Convergence criterion: 1
#> Tolerance for norm of final gradient: 1e-06
#> Mue of Armijo condition: 0.2
#> Minimum of step size value: 1e-10
#> Scaling factor: -1
#> 
#> 
#> Convergence reached in 5 iterations.
#> Number of function evaluations: 6 (6,6)
#> 
#> Maximum of log likelihood: -60.764
#> Norm of final gradient vector: 4.21004e-09
#> Last absolute change of function value: 8.73702e-11
#> Last relative change in parameters: 2.89833e-05
#> 
#> Idx Cat Term   Variable                      Coeff                    Error                  C/Error  Signif
#> ------------------------------------------------------------------------------------------------------------
#>   1   1  I     Intercept        0.4058623920111429       0.2195944416257131       1.8482361803260641  0.9354
#>   2   1  X     Vx               0.9012703568192907       0.2672256852091420       3.3726935946068162  0.9993
#> 
#> Log likelihood (starting values):     -69.3147180559945895 
#> Log likelihood (final estimates):     -60.7639890351186409 
#> 
#> ----------------------------------------------------------------------------
#> Current memory: 390032 bytes. Max memory used: 394519 bytes.
#> End of program. Sun Sep 20 03:44:26 2026
#> 
#> --- stderr ---
#> 
#>   Iter    Function Value     Norm of Gradient  Par Change   FCall
#>     1   6.9314718055995e+01  1.8886024681e+01          --       1 (1,1)
#>     2   6.0927896207017e+01  2.2170219810e+00  1.0000e+00       2 (2,2)
#>     3   6.0764474063353e+01  1.1593901524e-01  1.4768e-01       3 (3,3)
#>     4   6.0763989040428e+01  3.8218230965e-04  8.6593e-03       4 (4,4)
#>     5   6.0763989035119e+01  4.2100428432e-09  2.8983e-05       5 (5,5)
```
