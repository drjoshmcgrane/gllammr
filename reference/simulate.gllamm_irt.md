# Simulate response matrices from a fitted IRT model

Draws fresh person abilities from the fitted ability distribution and
samples item responses - dichotomous or polytomous - from the estimated
item parameters (parametric bootstrap / posterior predictive style).
Works for both EM and Laplace fits.

## Usage

``` r
# S3 method for class 'gllamm_irt'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  Fitted `gllamm_irt` object

- nsim:

  Number of replicates

- seed:

  Optional seed (stored as the `"seed"` attribute)

- ...:

  Unused

## Value

A list of `nsim` simulated response matrices, with a `"seed"` attribute.
