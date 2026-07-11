# Simulate from a fitted multinomial model

Draws fresh random effects for every term and samples categories from
the baseline-category logit probabilities.

## Usage

``` r
# S3 method for class 'gllamm_multinomial'
simulate(object, nsim = 1, seed = NULL, newdata = NULL, ...)
```

## Arguments

- object:

  Fitted `gllamm_multinomial` object

- nsim:

  Number of replicates

- seed:

  Optional seed (stored as the `"seed"` attribute)

- newdata:

  Optional new data frame

- ...:

  Unused
