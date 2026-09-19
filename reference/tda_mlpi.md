# 0-1 linear programming

Maximize a linear objective over 0-1 variables under linear constraints
– TDA's `mlpi` (CACM algorithm 449). The problem solved is \$\$\max\\
c'x + k \quad \mathrm{s.t.}\\ A x \ge b,\\ x \in \\0,1\\^n\$\$ with
integer coefficients (negative ones are fine). The command's source
states no objective; this is what it computes, established against
exhaustive search – see `examples/coverage/mlpiops.cf`.

## Usage

``` r
tda_mlpi(
  objective,
  constraints,
  bounds,
  constant = 0,
  max_solutions = 100,
  dir = tempfile("tda"),
  ...
)
```

## Arguments

- objective:

  integer coefficients `c`, one per variable (at least two variables).

- constraints:

  the constraint matrix `A`, one row per constraint (at least one),
  `length(objective)` columns.

- bounds:

  the right-hand sides `b`, one per constraint row.

- constant:

  `k`, added to the reported value (default 0).

- max_solutions:

  how many optimal solutions to return at most.

- dir:

  working directory.

- ...:

  passed to [`tda_run`](tda_run.md).

## Value

A list: `value`, the maximum; `solutions`, a matrix with one optimal 0-1
vector per row.

## See also

Other optimization: [`tda_boolean_min()`](tda_boolean_min.md),
[`tda_mlp()`](tda_mlp.md), [`tda_mls()`](tda_mls.md),
[`tda_mqap()`](tda_mqap.md), [`tda_mqp()`](tda_mqp.md)

## Examples

``` r
# maximize 2 x1 - 3 x2 + 4 x3 + x4  s.t.  sum(x) >= 2,
# -x1 + 2 x2 + x4 >= 1
tda_mlpi(c(2, -3, 4, 1),
        rbind(c(1, 1, 1, 1), c(-1, 2, 0, 1)),
        bounds = c(2, 1))
#> $value
#> [1] 5
#> 
#> $solutions
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    1    1
#> 
```
