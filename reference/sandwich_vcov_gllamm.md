# Cluster-robust (sandwich) covariance for GLLAMM fixed effects

Computes the cluster-robust covariance \\H^{-1} M H^{-1}\\ where \\H\\
is the observed information of the outer parameters and \\M = \sum_j s_j
s_j'\\ accumulates per-cluster score vectors. Because clusters are
independent, the Laplace marginal log-likelihood decomposes by cluster;
each cluster's score is evaluated by rebuilding its objective on the
cluster's rows and differentiating at the full-data estimates (with the
cluster's random effects re-profiled internally).

## Usage

``` r
sandwich_vcov_gllamm(object)
```

## Arguments

- object:

  A two-level GLMM fitted via
  [`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
  (gaussian, binomial, poisson, or Gamma family, single random-effects
  term)

## Value

Covariance matrix for the full outer parameter vector, with a `"fixed"`
attribute holding the fixed-effects block.

## See also

[`vcov.gllamm`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm-class.md)
with `type = "sandwich"`
