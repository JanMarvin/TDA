# TDA's manual, and two notes about the session

`tda_help` looks a keyword up in `tda.hlp`, the manual that ships with
TDA – the same text the standalone program prints for `help=...`. With
no keyword it lists every topic.

## Usage

``` r
tda_help(topic = "")

tda_time()

tda_memory()
```

## Arguments

- topic:

  a keyword. TDA appends its wildcard, so a prefix is enough; an empty
  topic lists every keyword.

## Value

The text TDA printed, invisibly; it is also printed.

## Details

`tda_time` and `tda_memory` report what TDA reports for its `time` and
`mem` commands.

## See also

Other TDA infrastructure:
[`tda_ccnt()`](https://janmarvin.github.io/TDA/reference/tda_ccnt.md),
[`tda_output()`](https://janmarvin.github.io/TDA/reference/tda_estimates.md),
[`tda_read_table()`](https://janmarvin.github.io/TDA/reference/tda_read_table.md),
[`tda_run()`](https://janmarvin.github.io/TDA/reference/tda_run.md),
[`tda_write_data()`](https://janmarvin.github.io/TDA/reference/tda_write_data.md)

## Examples

``` r
tda_help("rate")
#> rate: The rate command estimates transition rate models. The syntax is:
#> 
#>     rate(
#>         tp=...,         time periods
#>         xa(j,k)=...,    variables for A-term of model
#>         xb(j,k)=...,    variables for B-term of model
#>         xc(j,k)=...,    variables for C-term of model
#>         xd(j,k)=...,    variables for D-term of model
#>         deg=...,        degree of polynomial rate, def. 0
#>         kgam=...,       coefficient for gamma model, def. 1
#>         mix=...,        selection of mixture model
#>         grp=...,        groups for stratified Cox model
#>         rrisk,          print table with relative risks
#>         ppar=...,       output file with estimated parameters
#>         pcov=...,       output file with covariance matrix
#>         pres=...,       generalized residuals
#>         prate=...,      calculation of estimated rates
#>         tfmt=...,       print format for results, def. 10.4
#>         mfmt=...,       print format for pcov, pres and prate, def. 12.4
#>         dsv=...,        input file with starting values
#>         con=...,        linear parameter constraints
#>                         in addition, one can use all parameters to
#>                         control maximum likelihood estimation.
#>         mplog=...,      write loglikelihood value into matrix
#>         mppar=...,      write parameters into matrix
#>         mpcov=...,      write covariance matrix into matrix
#>         mpgrad=...,     write gradients into matrix
#>     ) = model_number;
#> 
#> For available models, see: rate model numbers
#> 
#> See also: frml, mpcov, mpgrad, mplog, mppar
tda_time()
#> Current time: Sun Sep 20 03:01:06 2026
tda_memory()
#> Current memory: 390032 bytes.
#> Currently requested memory: 390032 (390032) bytes.
#> Current memory: 390032 bytes. Max memory used: 390032 bytes.
```
