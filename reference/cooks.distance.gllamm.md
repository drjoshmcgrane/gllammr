# Group-level Cook's distance for GLLAMM models

Influence of each cluster on the fixed effects, computed by refitting
the model with the cluster deleted: D_j = (beta - beta\_(-j))' V^(-1)
(beta - beta\_(-j)) / p, with V the estimated covariance of the fixed
effects. Case deletion at the cluster level is the standard influence
measure for mixed models; observation-level deletion would break the
random-effects structure.

## Usage

``` r
# S3 method for class 'gllamm'
cooks.distance(model, max_groups = 50, ...)
```

## Arguments

- model:

  A fitted `gllamm` object (from
  [`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md))

- max_groups:

  Refuse to run for more clusters than this (each cluster costs one
  model refit); raise the limit explicitly for large data

- ...:

  Additional arguments (currently unused)

## Value

Named vector of Cook's distances, one per cluster. Clusters whose
deletion refit fails get `NA`.
