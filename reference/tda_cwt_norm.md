# Rescale case weights the way TDA's cwt(wnorm=) does

`tda_cwt_norm` reproduces TDA's `cwt(wnorm=s)=W;` (or the bare-flag
`cwt(wnorm)=W;`) as a plain vector transform, for `weights=` on
[`tda_qreg`](tda_qreg.md), [`tda_glm`](tda_glm.md),
[`tda_ple`](tda_ltb.md), [`tda_ltb`](tda_ltb.md) and
[`tda_rate`](tda_rate.md) – none of which expose `wnorm` as its
argument, since this is all it does: fitted coefficients come back
identical whether the weights are rescaled this way first or not, only
the standard errors change (rescaling is a constant multiplier on the
whole weighted log-likelihood, so it moves the curvature at the optimum,
not where the optimum is). `cwt(wnorm=s)=W;` on the raw weights and
plain `cwt=W;` on `tda_cwt_norm(W, s)` give the same fit, to the last
printed digit.

## Usage

``` r
tda_cwt_norm(w, s = length(w))
```

## Arguments

- w:

  a numeric vector of raw case weights.

- s:

  the target sum for the rescaled weights – TDA's `wnorm=s`. Left at its
  default, `length(w)` (the number of cases), it matches the bare-flag
  form, `cwt(wnorm)=W;`, TDA's "keep the same effective sample size"
  correction.

## Value

`w`, rescaled so it sums to `s`.

## Examples

``` r
w <- c(2, 2, 2, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10)  # sums to 74
tda_cwt_norm(w)         # cwt(wnorm)=W;    -- rescaled to sum 13 (length(w))
#>  [1] 0.245283 0.245283 0.245283 1.226415 1.226415 1.226415 1.226415 1.226415
#>  [9] 1.226415 1.226415 1.226415 1.226415 1.226415
tda_cwt_norm(w, 100)    # cwt(wnorm=100)=W;
#>  [1] 1.886792 1.886792 1.886792 9.433962 9.433962 9.433962 9.433962 9.433962
#>  [9] 9.433962 9.433962 9.433962 9.433962 9.433962
```
