# Test Item Covariate Effects

Test whether adding item covariates significantly improves model fit

## Usage

``` r
test_item_covariates(
  response_matrix,
  item_data,
  difficulty_formula = ~1,
  discrimination_formula = ~1,
  model = c("Rasch", "2PL", "GRM"),
  ...
)
```

## Arguments

- response_matrix:

  Matrix of item responses

- item_data:

  Data frame of item-level covariates

- difficulty_formula:

  Formula for difficulty regression

- discrimination_formula:

  Formula for discrimination regression

- model:

  IRT model type

- ...:

  Additional arguments passed to fit_eirt

## Value

A list with:

- full_model:

  EIRT model with covariates

- null_model:

  EIRT model without covariates (intercept only)

- comparison:

  Model comparison results

## Examples

``` r
if (FALSE) { # \dontrun{
item_data <- data.frame(
  word_freq = rnorm(20),
  length = rpois(20, 5)
)

# Test if word frequency matters
result <- test_item_covariates(
  responses,
  item_data,
  difficulty_formula = ~ word_freq + length,
  model = "Rasch"
)

print(result$comparison)
} # }
```
