# Quadratic assignment problem

Approximate solution of the quadratic assignment problem – TDA's `mqap`,
adapted from CACM algorithm 608 (D.H. West). Tries to find a permutation
`p` of `1:n` minimizing \$\$\sum_i \left( C\_{i,p(i)} + \sum_j F\_{ij}
D\_{p(i),p(j)} \right)\$\$ This is a heuristic, not an exact solver: the
value returned is an upper bound on the true minimum, and need not be
optimal for larger `n` (see `examples/coverage/mqapops.cf` for the
objective and small hand-verified instances). `flows` and `distances`
must have a zero main diagonal; a nonzero diagonal is silently reset to
zero with a warning from TDA itself.

## Usage

``` r
tda_mqap(flows, distances, costs = NULL, dir = tempfile("tda"), ...)
```

## Arguments

- flows:

  the flow matrix `F` (n x n).

- distances:

  the distance matrix `D` (n x n).

- costs:

  the fixed placement cost matrix `C` (n x n), default all zero.

- dir:

  working directory.

- ...:

  passed to
  [`tda_run`](https://janmarvin.github.io/TDA/reference/tda_run.md).

## Value

A list: `value`, the objective at the returned permutation;
`permutation`, integer vector `p` with `p[i]` the location assigned to
facility `i`.

## See also

Other optimization:
[`tda_boolean_min()`](https://janmarvin.github.io/TDA/reference/tda_boolean_min.md),
[`tda_mlp()`](https://janmarvin.github.io/TDA/reference/tda_mlp.md),
[`tda_mlpi()`](https://janmarvin.github.io/TDA/reference/tda_mlpi.md),
[`tda_mls()`](https://janmarvin.github.io/TDA/reference/tda_mls.md),
[`tda_mqp()`](https://janmarvin.github.io/TDA/reference/tda_mqp.md)

## Examples

``` r
F <- matrix(c(0,5,2,4, 5,0,3,0, 2,3,0,0, 4,0,0,0), 4, 4, byrow = TRUE)
D <- matrix(c(0,8,15,13, 8,0,9,10, 15,9,0,17, 13,10,17,0), 4, 4,
            byrow = TRUE)
tda_mqap(F, D)
#> $value
#> [1] 270
#> 
#> $permutation
#> [1] 2 1 4 3
#> 
```
