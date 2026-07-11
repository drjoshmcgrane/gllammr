# Simulate from a fitted ordinal model

Draws fresh random effects from their estimated distribution for each
replicate (population-level simulation) and samples categories from the
implied response distribution.

## Usage

``` r
# S3 method for class 'gllamm_ordinal'
simulate(object, nsim = 1, seed = NULL, newdata = NULL, ...)
```

## Arguments

- object:

  A fitted `gllamm_ordinal` object

- nsim:

  Number of simulations (default: 1)

- seed:

  Optional random seed (stored in the `"seed"` attribute)

- newdata:

  Optional new data frame with the model covariates and grouping
  variables

- ...:

  Additional arguments (currently unused)

## Value

A data frame with `nsim` columns of simulated category responses
(integer codes 1..K), with a `"seed"` attribute.
