# Simulate response matrices from a fitted latent class model

Draws class memberships from the estimated prevalences and item
responses from the class-conditional distributions (binary, categorical,
and gaussian indicators).

## Usage

``` r
# S3 method for class 'gllamm_lca'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  Fitted `gllamm_lca` object

- nsim:

  Number of replicates

- seed:

  Optional seed (stored as the `"seed"` attribute)

- ...:

  Unused
