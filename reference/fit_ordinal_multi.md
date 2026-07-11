# Ordinal model with multiple random-effects terms

Internal engine behind fit_ordinal for crossed/multiple random-effects
terms (links logit, probit, acl, crl_forward, crl_backward). The
random-effects layout mirrors fit_tmb_gllamm_multi: one combined sparse
Z, term-major u, per-term (possibly correlated) covariance.

## Usage

``` r
fit_ordinal_multi(
  formula,
  data,
  link,
  link_code,
  parsed,
  model_data,
  y_numeric,
  n_categories,
  category_labels,
  weights,
  control
)
```
