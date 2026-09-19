# Where the manual's remaining commands live in R

Working through TDA's manual section by section, some commands have no
function of their own here. None of them is missing functionality; each
is covered in one of three ways.

## Details

**Subsumed by R itself.** The matrix definition commands exist because
TDA has no other way to get a matrix into a session; R does. `mdef`
(literal) is [`matrix()`](https://rdrr.io/r/base/matrix.html); `mdefc`
(constant) is `matrix(d, m, n)`; `mdefi` (rectangular identity) is
`diag(1, m, n)`; `mdeff` (from file) is `as.matrix(read.table(...))` or
[`tda_read_table`](https://janmarvin.github.io/TDA/reference/tda_read_table.md);
`mpr` (print a matrix expression) is printing the return value of any
matrix wrapper. The expression operators of section 5.2.5 – `nd`, `ndi`,
`td`, `chd`, `fd`, `mr`, `ndf`, `poisson`, `negbin`, `bivn`, `jul`,
`july`, `julm`, `juld` and relatives – are `pnorm`, `qnorm`, `pt`,
`pchisq`, `pf`, the Mill's ratio `dnorm(x)/(1 - pnorm(x))`, `dnorm`,
`dpois`/`dnbinom` (TDA returns the log density), `mvtnorm::pmvnorm`, and
R's date arithmetic. The engine's implementations of all of these are
compared against R's on every run of the differential fuzzer, so using
either side gives the same numbers; to evaluate one of TDA's operators
directly, pass the expression to
[`tda_evalf`](https://janmarvin.github.io/TDA/reference/tda_integrate.md).

**Delivered by another wrapper.** `mdefg` (adjacency matrix from the
current graph) is `tda_g_edges(g, format = "full_square")`, with `-1`
marking absent edges. `mbrr` is
[`tda_brr`](https://janmarvin.github.io/TDA/reference/tda_brr.md).
`mple` is
[`tda_mple`](https://janmarvin.github.io/TDA/reference/tda_mple.md).
`mlp`'s raw-tableau form is assembled internally by
[`tda_mlp`](https://janmarvin.github.io/TDA/reference/tda_mlp.md) from
the objective, constraints, and bounds.

**Reached through the session they belong to.** `mdefb` (a block of the
data matrix) requires a `dblock` context and is used through the
`commands=` argument of
[`tda_dblock`](https://janmarvin.github.io/TDA/reference/tda_dblock.md).
`mplog`, `mppar`, `mpcov` and `mpgrad` are parameters OF estimation
commands, not commands; their contents are what the fitted objects
already carry as `$logLik`, `$estimates` and `$vcov` (case-specific
gradients, `mpgrad`, can be requested through a command's `options=`
where supported).

Chapter 5 names that are not standalone commands:

- `mplog`, `mppar`, `mpcov`, `mpgrad`:

  parameters of the estimation commands that write the log likelihood,
  parameters, covariance matrix or gradient into a named TDA matrix. The
  R wrappers return the same quantities directly (`$logLik`,
  `$estimates`, `$vcov`), so these are only needed inside a hand-written
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md)
  matrix pipeline, where they pass through `options=` as usual.

- `mdef`, `mpr`, `mdefc`, `mdefi`, `mdeff`:

  matrix definition and printing, which R does natively: `matrix(...)`,
  [`print()`](https://rdrr.io/r/base/print.html), `matrix(d, m, n)`,
  `diag(1, m, n)` and `as.matrix(read.table(fname))`; every matrix
  wrapper accepts plain R matrices. `mdefg` (a graph's adjacency matrix)
  does compute something and has its wrapper,
  [`tda_mdefg`](https://janmarvin.github.io/TDA/reference/tda_mdefg.md).

- `mdefb`:

  extracts one dblock-defined block of the data matrix into a matrix –
  in R, subset the data frame; in a
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md)
  script it works after a `dblock` command as documented.

- `mbrr`:

  the manual's name for `brr`; see
  [`tda_brr`](https://janmarvin.github.io/TDA/reference/tda_brr.md).

- `sort` (command):

  see
  [`tda_esort`](https://janmarvin.github.io/TDA/reference/tda_esort.md)
  and the `sort=` option of `pdata` in
  [`tda_read_sys`](https://janmarvin.github.io/TDA/reference/tda_read_sys.md).

- `lag`, `suc`, `snum`, `rank`, `sort` (operator):

  operators of TDA's expression language, used inside `nvar` variable
  definitions – e.g. `L[6.0]=lag(Y,1)`, `R[6.0]=rank(Y)` – in a
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md)
  script. `snum`, `rank` and `sort` take the variable as an argument.
