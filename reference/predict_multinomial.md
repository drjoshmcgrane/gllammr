# Predict method for multinomial models

Obtain predictions from a fitted multinomial logit model with random
effects

## Usage

``` r
# S3 method for class 'gllamm_multinomial'
predict(
  object,
  newdata = NULL,
  type = c("class", "probs", "marginal"),
  n_sim = 1000,
  ...
)
```

## Arguments

- object:

  A fitted multinomial model

- newdata:

  Optional new data frame for predictions

- type:

  Type of prediction:

  class

  :   Predicted class (modal category)

  probs

  :   Conditional probabilities for each category

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

- marginal: Matrix of marginal probabilities (n_obs × n_categories)
