# Simulate indicator data from a fitted SEM

Draws from the fitted multivariate-normal implied distribution (means +
implied covariance) of the indicators.

## Usage

``` r
# S3 method for class 'gllamm_sem'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  Fitted `gllamm_sem` object (ML method)

- nsim:

  Number of replicates

- seed:

  Optional seed (stored as the `"seed"` attribute)

- ...:

  Unused
