# Simulate from a fitted mixed-response model

Draws fresh shared random intercepts and simulates every outcome from
its fitted family (gaussian/binomial/poisson), returning a list of data
frames.

## Usage

``` r
# S3 method for class 'gllamm_mixed'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  Fitted `gllamm_mixed` object

- nsim:

  Number of replicates

- seed:

  Optional seed (stored as the `"seed"` attribute)

- ...:

  Unused
