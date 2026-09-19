# Maximum likelihood with a user-defined likelihood

`fml` maximises a log-likelihood you write yourself, which is what TDA
offers that a fixed catalogue of models cannot. The likelihood is
written in an expression language, evaluated inside TDA rather than as
an R function – but as of this version, that expression can be written
in ordinary R syntax and translated automatically, rather than TDA's.

## Usage

``` r
tda_fml(
  definitions,
  data,
  start = NULL,
  control = NULL,
  constraints = NULL,
  residuals = FALSE,
  residual_vars = NULL,
  protocol = FALSE,
  options = list(),
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- definitions:

  either an unevaluated [`{ }`](https://rdrr.io/r/base/Paren.html) block
  of R assignments (see Details) or a character vector of TDA's
  assignments, the last defining `fn`, the contribution of one case to
  the log-likelihood. Names that appear on the right but are neither
  variables in `data` nor defined earlier are the parameters to be
  estimated.

- data:

  a data frame; its columns are the variables the definitions can refer
  to.

- start:

  optional named vector of starting values, e.g. `c(b0 = 0, b1 = 1)`.

- control:

  convergence settings from [`tda_control`](tda_control.md).

- constraints:

  optional linear constraints on the parameters – `fml`'s `con=`, the
  same machinery [`tda_qreg`](tda_qreg.md)'s `constraints` uses.
  **Refers to parameters by their position, `b1`/`b2`/..., 1-based, in
  the order they are first introduced across `definitions` – never by
  whatever name was actually given them.** A parameter named `b0` in
  `definitions` is still `b1` here if it is the first one introduced;
  TDA's error for a literal `"b0 = ..."` constraint is “Error in
  parameter index”, since `b0` is not a valid 1-based position. Run once
  without `constraints` and check
  [`tda_estimates()`](tda_estimates.md)'s `Idx` column for which
  position is which parameter before writing one. Must be strictly fewer
  constraints than parameters (TDA's rule; its error reads “number of
  constraints should be less than number of parameters”).

- residuals:

  ask TDA to also compute, per case, its contribution to the
  log-likelihood at the converged parameters – `fml`'s `pres=`. **Not a
  classical observed-minus-fitted residual** – `fml` has no built-in
  notion of a fitted value to subtract from, unlike
  [`tda_lsreg`](tda_lsreg.md)'s `residuals`. For
  `fn = -0.5 * (y - xb)^2`, say, this is that expression's value per
  case, not `y - xb` itself.

- residual_vars:

  with `residuals = TRUE`, extra columns from `data` to write alongside
  `fn`'s value, one column per case, in the order given – `fml`'s `v=`.
  `$residuals` then comes back as a data frame (`fn` plus each named
  column) instead of a plain vector. Has no effect without
  `residuals = TRUE`.

- protocol:

  ask TDA to also write its iteration-by-iteration diagnostic log –
  `fml`'s `prot=`: the iteration table (function value, gradient norm,
  parameter change per step), the parameter vector, and the covariance
  matrix at convergence. Not one consistent shape across every algorithm
  and situation, so returned as `$protocol`, the file's text lines,
  rather than force-parsed into a table that would not fit every case.

- options:

  a named list of further TDA options, passed through.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

An object of class `tda_fit`. With `residuals = TRUE`, also carries
`$residuals` (see `residuals` above).

## Writing the likelihood in R syntax

Give `definitions` as an unevaluated
[`{ }`](https://rdrr.io/r/base/Paren.html) block of plain R assignments
instead of a character vector, and it is translated to TDA's own
expression syntax for you – never actually run as R code (`DES`, `rate`
and the rest only ever exist inside TDA, so it cannot be).
`ifelse(cond, yes, no)` becomes TDA's `if()` (TDA's `if` takes a 2- or
3-argument `if(cond, then)`/`if(cond, then, else)` shape matching
`ifelse` exactly, not R's control-flow `if`), and R's comparison
operators (`<`, `>`, `<=`, `>=`, `==`, `!=`) become TDA's comparison
functions (`lt`, `gt`, `le`, `ge`, `eq`, `ne`). Arithmetic
(`+ - * / ^`), `&`/`|`, and ordinary function calls (`exp`, `log`,
`negbin`, ...) already mean the same thing in both and pass through
unchanged.

**These six TDA functions are not the strict comparisons their R
equivalents are.** Every one of them – `lt`, `le`, `gt`, `ge`, `eq`,
`ne` alike – carries a fixed tolerance, `EPSI2 = 1000 * DBL_EPSILON`
(about `2.22e-13`): `le`, say, computes `x <= y + EPSI2`, not `x <= y`.
This can matter on real data: a case in `examples/exam/rrdat.1` has
`PRESN/PRES - 1` equal to `0.2` exactly in exact arithmetic, and the two
land on different sides of `0.2` once represented as IEEE 754 doubles –
`ge(x, 0.2)` inside TDA still calls it true, R's own `x >= 0.2` on the
identical computation does not. Applied consistently across all six, at
least, but nowhere in TDA's documentation: not in the manual, not in
[`tda_help()`](tda_help.md)'s operator table, which describes `le` as
plainly as “true if less than, or equal” with no footnote at all. A
hand-written R replication of a TDA calculation involving a comparison
near a round number (`==`, or a boundary like `0.2` here) can therefore
disagree with TDA's output at that boundary, for a reason with nothing
to do with which language computed it – this is a property of TDA's
expression evaluator, present whether the comparison was written by hand
in TDA's syntax or reached it through this translation.


    tda_fml({
        rate = exp(a0 + COHO2 * a1 + COHO3 * a2 + W * a3)
        l1 = ifelse(DES, log(rate), 0)
        fn = l1 - rate * DUR
    }, data = d, start = c(-4, 0, 0, 0))

The older, direct form – a character vector of TDA's assignment syntax –
still works exactly as before, for anything the translator does not
(yet) cover, or for pasting a `.cf` file's text in directly.


    tda_fml(c("xb = b0 + x * b1",
              "fn = -0.5 * (y - xb)^2"), data = d)

estimates `b0` and `b1` by least squares written as a likelihood. TDA's
functions are available, so `fn = negbin(alpha, gamma, y)` and the rest
work.

## See also

Other sequence analysis: [`tda_evalfi()`](tda_evalfi.md),
[`tda_frml()`](tda_frml.md), [`tda_seq_info()`](tda_seq_info.md),
[`tda_seqgc()`](tda_seqgc.md), [`tda_seqm()`](tda_seqm.md),
[`tda_seqmd()`](tda_seqmd.md), [`tda_seqpe()`](tda_seqpe.md),
[`tda_seqpm()`](tda_seqpm.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x = 1:20, y = 2 + 0.5 * (1:20) + rnorm(20))
tda_fml(c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"), d)
#> Call: tda_fml(definitions = c("xb = b0 + x * b1", "fn = -0.5 * (y - xb)^2"), 
#>     data = d)
#> 
#> Cases: 20
#> logLik (starting values): -651.8385
#> logLik:  -7.768427
#> Converged in 2 iterations
#> 
#>  Idx Parameter  Coeff  Error C/Error Signif
#>    1        b0 1.9639 0.4645  4.2277 1.0000
#>    2        b1 0.5216 0.0388 13.4504 1.0000

# the same model, written in R syntax instead
tda_fml({
    xb = b0 + x * b1
    fn = -0.5 * (y - xb)^2
}, d)
#> Call: tda_fml(definitions = {
#>     xb = b0 + x * b1
#>     fn = -0.5 * (y - xb)^2
#> }, data = d)
#> 
#> Cases: 20
#> logLik (starting values): -651.8385
#> logLik:  -7.768427
#> Converged in 2 iterations
#> 
#>  Idx Parameter  Coeff  Error C/Error Signif
#>    1        b0 1.9639 0.4645  4.2277 1.0000
#>    2        b1 0.5216 0.0388 13.4504 1.0000

# a compound Poisson (negative binomial) count regression -- TDA has no
# separate command for this (unlike its plain Poisson, which is a glm()
# family), but negbin(alpha, gamma, k) is a built-in log-likelihood
# function in its expression language, reachable through fml like any
# other: gamma is the mean, alpha = 1/sigma the dispersion
set.seed(7)
d2 <- data.frame(x = rnorm(200))
d2$ndi <- rnbinom(200, mu = exp(0.5 + 0.8 * d2$x), size = 2)
coef(tda_fml({
    xb = b0 + x * b1
    gamma = exp(xb)
    sigma = exp(sig)
    alpha = 1 / sigma
    fn = negbin(alpha, gamma, ndi)
}, d2, start = list(b0 = 0, b1 = 0, sig = 0)))
#>         b0         b1        sig 
#>  0.4208311  0.8130428 -0.4504466 
```
