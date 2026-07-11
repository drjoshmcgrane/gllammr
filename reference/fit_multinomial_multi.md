# Multinomial model with multiple random-effects terms

Internal engine behind fit_multinomial for crossed/multiple random
effects. Layout mirrors fit_ordinal_multi / fit_tmb_gllamm_multi; as in
the single-term template, the random effects act as a common shifter on
every non-reference category.

## Usage

``` r
fit_multinomial_multi(
  formula,
  data,
  reference,
  parsed,
  model_data,
  y_numeric,
  n_categories,
  category_labels,
  control
)
```
