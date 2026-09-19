# Methods for fitted TDA models

Standard methods for the objects the model functions return. See
[`tda_rate`](tda_rate.md) for what those objects are.

## Usage

``` r
# S3 method for class 'tda_mds'
print(x, ...)

# S3 method for class 'tda_freq2'
print(x, ...)

# S3 method for class 'tda_ragged'
print(x, n = 10L, ...)

# S3 method for class 'tda_graph'
print(x, ...)

# S3 method for class 'tda_iv'
print(x, ...)

# S3 method for class 'tda_ivstat'
print(x, ...)

# S3 method for class 'tda_ps'
print(x, ...)

# S3 method for class 'tda_ps'
plot(x, ...)

# S3 method for class 'summary.tda_fit'
print(x, ...)

# S3 method for class 'tda_ple'
print(x, n = 10L, ...)

# S3 method for class 'tda_ltb'
summary(object, n = 10L, ...)

# S3 method for class 'summary.tda_ltb'
print(x, ...)

# S3 method for class 'tda_fit'
print(x, ...)

# S3 method for class 'tda_ltb'
print(x, ...)

# S3 method for class 'tda_table'
print(x, ...)

# S3 method for class 'tda_table'
plot(x, ...)

# S3 method for class 'tda_expr'
print(x, ...)

# S3 method for class 'tda_sdvd'
plot(x, ...)

# S3 method for class 'tda_sd_ps3'
plot(x, newpage = TRUE, ...)

# S3 method for class 'tda_spatial'
print(x, ...)

# S3 method for class 'tda_result'
print(x, ...)

# S3 method for class 'tda_result'
summary(object, ...)

# S3 method for class 'tda_result'
logLik(object, ...)

# S3 method for class 'tda_result'
coef(object, ...)
```

## Arguments

- x, object:

  a fitted model or result.

- ...:

  passed on.

- n:

  number of rows to show.

## Value

As for the corresponding generic.

## Examples

``` r
set.seed(41)
d <- data.frame(t = round(rexp(50, 0.1), 1) + 0.5,
                 s = rbinom(50, 1, 0.85), x = rnorm(50))
f <- tda_rate(Surv(t, s) ~ x, d, model = "exponential")
f            # print.tda_fit
#> Call: tda_rate(formula = Surv(t, s) ~ x, data = d, model = "exponential")
#> 
#> Episodes: 50
#> Model:   exponential (TDA 2)
#> logLik (starting values): -141.7427
#> logLik:  -141.4876
#> Converged in 4 iterations
#> 
#>  Idx SN Org Des MT Variable   Coeff  Error  C/Error Signif
#>    1  1   0   1  A Constant -2.5480 0.1588 -16.0415 1.0000
#>    2  1   0   1  A        x  0.1087 0.1526   0.7122 0.5237
coef(f)
#>   Constant          x 
#> -2.5479654  0.1086693 
summary(f)   # summary.tda_fit / print.summary.tda_fit
#> Call: tda_rate(formula = Surv(t, s) ~ x, data = d, model = "exponential")
#> 
#> Episodes: 50
#> Model:   exponential
#> Converged in 4 iterations
#> 
#>          Estimate Std. Error  C/Error Pr(>|t|)    
#> Constant -2.54797    0.15884 -16.0415   <2e-16 ***
#> x         0.10867    0.15258   0.7122   0.4763    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> logLik (starting values): -141.7427
#> logLik: -141.4876   df: 2   AIC: 286.9752
#> Null model (exponential, constant rate) logLik: -141.7427
#> LR test against it: chi2 = 0.510277 on 1 df, p = 0.475
```
