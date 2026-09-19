# TDA's random number generator

`tda_rng` is a stateful generator reproducing TDA's `rd`/ `rdn` bit for
bit. `tda_runif`/`tda_rnorm` are convenience wrappers around a fresh
generator, for when only one of `rd`/ `rdn` is needed; `tda_rng` itself
is what a `nvar()` block using *both* together needs, since they draw
from one continuous, shared stream, not two independent ones.

## Usage

``` r
tda_rng(seed = 13421773)

tda_runif(n, seed = 13421773)

tda_rnorm(n, seed = 13421773)
```

## Arguments

- seed:

  the starting value. TDA's default is 13421773.

- n:

  how many numbers to draw, for `tda_runif`/`tda_rnorm`.

## Value

`tda_rng` returns a list with two functions, `$rd(n = 1, a = 0, b = 1)`
and `$rdn(n = 1)`, each drawing `n` values from the shared stream (a
single value by default, for calling one at a time in an interleaved
loop) and returning a length-`n` numeric vector (or a single number, for
the default `n = 1`). `$rd()`'s `a`/`b` are the range – TDA's
`rd(a, b)`; `$rd()` alone is `[0, 1)`, matching `rd` with no arguments.
`tda_runif`/`tda_rnorm` return a plain numeric vector of length `n`
directly.

## Details

`random1()` is the multiplicative generator behind `rd()`: the seed is
multiplied by 3125 modulo 2^26, in three steps of 25, 25 and 5.
`normal()` (behind `rdn()`; `rdn1` is a different pseudo-variable,
backed by a different function, `normal1()`, not this one) draws
`random1()` values in a rejection loop (Marsaglia's polar method, J.R.
Bell's ACM Algorithm 334) and returns *two* standard normal deviates per
accepted draw, the second cached and returned on the next call with no
further `random1()` draws consumed – so the number of underlying uniform
draws per `rdn()` call varies, and cannot be precomputed or batched
independently of `rd()` calls interleaved with it.

This is why `nvar()`'s evaluation order matters: it evaluates one case
fully, one variable at a time in definition order, before moving to the
next case – not all of one variable's values before the next. A block
defining `X = rd` then `E = rdn` draws `X`'s case 1, then `E`'s case 1
(however many `random1()` calls that costs), then `X`'s case 2, and so
on – `tda_runif(256)` followed separately by `tda_rnorm(256)` does not
reproduce this at all, even with the same seed, since it draws all of
`X` before touching `E`'s own share of the stream. Reproducing that
exact case-by-case interleaving is what `tda_rng`'s stateful
`$rd()`/`$rdn()` are for: call them in the same order, one case at a
time, and the shared state advances exactly as TDA's evaluation loop
would.

It matters because several of the shipped examples select a subsample
with `isel = le(rd(0,1), p)`, or generate simulated data with `rd`/`rdn`
together (`doc/npreg1.cf`'s worked example, say). Reproducing their
numbers needs the same stream, and R's generators will not do – this is
not a question of setting a seed. `tda_runif`/`tda_rnorm` reproduce
`random1()`/`normal()` bit for bit.

## Examples

``` r
tda_runif(5)
#> [1] 0.02910383 0.94947018 0.09430404 0.70012523 0.89135106

# X = rd, E = rdn in the same nvar() block -- interleaved, not two
# separate streams; compare doc/npreg1.cf's worked example.
gen <- tda_rng()
n <- 256
X <- E <- numeric(n)
for (i in seq_len(n)) {
    X[i] <- gen$rd()
    E[i] <- 0.32 * gen$rdn()
}
```
