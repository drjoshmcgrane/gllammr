# Predict method for IRT models

Obtain predictions from a fitted IRT model

## Usage

``` r
# S3 method for class 'gllamm_irt'
predict(
  object,
  newdata = NULL,
  type = c("probability", "ability", "marginal"),
  ability = NULL,
  n_sim = 1000,
  ...
)
```

## Arguments

- object:

  A fitted IRT model (class gllamm_irt)

- newdata:

  Optional specification of items/persons for predictions. Can be:

  - NULL: Predict for all original items

  - Numeric vector: Item indices to predict

  - Data frame with 'item' column: Specific items

- type:

  Type of prediction:

  probability

  :   Item response probabilities P(Y=1\|θ)

  ability

  :   Person ability estimates (θ)

  marginal

  :   Marginal item response probabilities E\[P(Y=1\|θ)\]

- ability:

  Optional vector of ability values for which to compute probabilities.
  If NULL and type="probability", uses estimated abilities from fitted
  model.

- n_sim:

  Number of Monte Carlo samples for marginal predictions (default: 1000)

- ...:

  Additional arguments (currently unused)

## Value

Depends on `type`:

- probability: Matrix of probabilities (n_persons × n_items) or vector
  if ability specified

- ability: Vector of person abilities

- marginal: Vector of marginal probabilities (one per item)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit 2PL model
responses <- matrix(rbinom(1000, 1, 0.6), 100, 10)
fit <- fit_irt(responses, model = "2PL")

# Person abilities
abilities <- predict(fit, type = "ability")

# Item response probabilities for each person
probs <- predict(fit, type = "probability")

# Marginal item response probabilities (population-level)
marg_probs <- predict(fit, type = "marginal")
} # }
```
