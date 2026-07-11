# EM estimation for latent class models

Classic EM for finite mixtures with conditionally independent indicators
(the poLCA algorithm). All M-steps are closed-form (weighted proportions
for binary/categorical indicators, weighted moments for gaussian ones),
so each iteration is one N x K posterior computation plus a handful of
cross-products - the E-step runs on BLAS matrix products. When
`order_edges` is supplied, binary probabilities and gaussian means are
additionally constrained to be nondecreasing along the given class
partial order via isotonic regression in the M-step.

## Usage

``` r
fit_lca_em(
  Y,
  nclass,
  item_type,
  n_cats,
  weights = NULL,
  n_starts = 3,
  max_iter = 1000,
  tol = 1e-08,
  order_edges = NULL,
  item_edges = NULL,
  structure = "free"
)
```
