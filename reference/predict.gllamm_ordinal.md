# Predict method for ordinal models

Obtain predictions from a fitted ordinal regression model

## Usage

``` r
# S3 method for class 'gllamm_ordinal'
predict(
  object,
  newdata = NULL,
  type = c("class", "probs", "cumprobs", "marginal"),
  n_sim = 1000,
  ...
)
```

## Arguments

- object:

  A fitted ordinal model (class gllamm with ordinal_family)

- newdata:

  Optional new data frame for predictions

- type:

  Type of prediction:

  class

  :   Predicted class (modal category)

  probs

  :   Conditional probabilities for each category

  cumprobs

  :   Conditional cumulative probabilities

  marginal

  :   Marginal probabilities (population-averaged)

- n_sim:

  Number of Monte Carlo samples for marginal predictions (default: 1000)

- ...:

  Additional arguments (currently unused)

## Value

Depends on `type`:

- class: Vector of predicted classes

- probs: Matrix of probabilities (n_obs × n_categories)

- cumprobs: Matrix of cumulative probabilities

- marginal: Matrix of marginal probabilities (n_obs × n_categories)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit ordinal model
fit <- gllamm(rating ~ temp + (1 | judge),
              data = wine,
              family = ordinal(link = "logit"))

# Predicted classes
pred_class <- predict(fit, type = "class")

# Conditional probabilities
pred_probs <- predict(fit, type = "probs")

# Marginal probabilities (population-averaged)
pred_marg <- predict(fit, type = "marginal")
} # }
```
