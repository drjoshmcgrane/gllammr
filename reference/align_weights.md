# Align observation weights with listwise-deleted model data

make_model_matrices() drops rows with missing values in any model
variable; user-supplied weights validated against the original data
length must be subset to the retained rows. Level-specific weight lists
are passed through untouched (they reference data columns and are
resolved downstream).

## Usage

``` r
align_weights(weights, model_data)
```
