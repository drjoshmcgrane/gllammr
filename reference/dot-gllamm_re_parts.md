# Per-term random-effects pieces of a fitted GLMM

Normalizes the two stored shapes (single-term: per-group list of
coefficient vectors; multi-term: per-term list of group x coef matrices)
into aligned per-term lists of design matrices, group indices,
covariance matrices, and BLUP matrices.

## Usage

``` r
.gllamm_re_parts(object, data)
```
