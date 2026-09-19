# Extract a command's result section from TDA output

TDA's protocol is structured: a header, then each command's echo
followed by its result, with separator lines between commands.
`tda_payload` returns the result section of the last command (or the
`n`-th from the end), without the surrounding protocol, which is what an
example or a quick look usually wants.

## Usage

``` r
tda_payload(output, n = 1L)
```

## Arguments

- output:

  a character vector of protocol lines, or a `tda_run` result (its
  `$output` is used).

- n:

  which command's section, counted from the end; default 1, the last
  command.

## Value

the lines of that section, invisibly empty if the output has no
separator structure.

## Examples

``` r
r <- tda_run("int(ab=0,1, fmt=12.8) = x*x;")
cat(tda_payload(r), sep = "\n")   # the integral, 1/3
#> int(...)=...
#> Numerical integration. Current memory: 390032 bytes.
#> Function definition:
#> fn      = x*x
#> Function argument: x
#> Integration interval: 0 -- 1
#> Method 1 (QNG).
#> Relative error:  1.00000e-04
#> Approximation:   0.33333333 
#> Number of function calls: 21
```
