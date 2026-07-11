# Predict Item Difficulties from Covariates

Compute predicted difficulties based on item covariate model

## Usage

``` r
predict_difficulty(object, newdata = NULL)
```

## Arguments

- object:

  A gllamm_eirt object

- newdata:

  Optional new item data for predictions

## Value

Vector of predicted difficulties (fixed effects only, no residuals)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_eirt(responses, item_data,
                difficulty_formula = ~ word_freq)

# Predicted difficulties for fitted data
pred_diff <- predict_difficulty(fit)

# Predicted difficulties for new items
new_items <- data.frame(word_freq = c(-1, 0, 1))
pred_new <- predict_difficulty(fit, newdata = new_items)
} # }
```
