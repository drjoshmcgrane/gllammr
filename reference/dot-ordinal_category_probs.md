# Category probabilities for every ordinal link

Mirrors the likelihood construction in gllamm_ordinal.hpp. For the PPO
link `eta` must be an n x (K-1) matrix of per-threshold linear
predictors; for all other links a length-n vector.

## Usage

``` r
.ordinal_category_probs(eta, thresholds, link, K)
```
