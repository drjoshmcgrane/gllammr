# Simulate from a fitted NPML model

Draws group-level intercepts from the estimated discrete (mass-point)
distribution and responses from the fitted family.

## Usage

``` r
# S3 method for class 'gllamm_npml'
simulate(object, nsim = 1, seed = NULL, newdata = NULL, ...)
```

## Arguments

- object:

  Fitted `gllamm_npml` object

- nsim:

  Number of replicates

- seed:

  Optional seed (stored as the `"seed"` attribute)

- newdata:

  Optional new data frame

- ...:

  Unused
