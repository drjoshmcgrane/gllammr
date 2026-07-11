# Parse observation-level and group-level survey weights

GLLAMM-style level-specific weights: `weights` may be a numeric vector
(level-1 / observation weights, as before) or a list with elements
`level1` (length n_obs) and/or `level2` (one weight per group, or per
observation but constant within group). Level-2 weights scale each
group's full likelihood contribution including its random-effects prior
(pseudo-likelihood for two-stage sampling designs).

## Usage

``` r
parse_level_weights(weights, n_obs, groups, n_groups)
```

## Arguments

- weights:

  NULL, numeric vector, or list(level1=, level2=)

- n_obs:

  Number of observations

- groups:

  0-indexed group index per observation

- n_groups:

  Number of groups

## Value

list(level1 = numeric(n_obs), level2 = numeric(n_groups))
