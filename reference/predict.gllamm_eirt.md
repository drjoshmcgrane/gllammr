# Predict method for EIRT models

Obtain predictions from a fitted Explanatory IRT model.

## Usage

``` r
# S3 method for class 'gllamm_eirt'
predict(
  object,
  newdata = NULL,
  type = c("probability", "ability", "difficulty", "discrimination", "marginal"),
  ability = NULL,
  n_sim = 1000,
  ...
)
```

## Arguments

- object:

  A fitted EIRT model (class gllamm_eirt)

- newdata:

  Optional data frame with item covariates for new items. Must include
  all variables used in the difficulty and discrimination formulas. For
  polytomous models, new-item predictions use the predicted item
  location with step deviations at their population value of zero
  (PCM/GPCM) or the threshold regression (LPCM); GRM threshold spacing
  is item-specific and cannot be predicted for new items.

- type:

  Type of prediction:

  probability

  :   Item response probabilities. Dichotomous: persons x items matrix
      of P(Y=1). Polytomous: list of persons x categories matrices, one
      per item.

  ability

  :   Person ability estimates (theta)

  difficulty

  :   Predicted item difficulties

  discrimination

  :   Predicted item discriminations

  marginal

  :   Marginal response probabilities, integrating ability over its
      population distribution. Dichotomous: vector of E\[P(Y=1)\].
      Polytomous: items x categories matrix of E\[P(Y=k)\].

- ability:

  Optional vector of ability values. If NULL, uses estimated abilities.

- n_sim:

  Number of Monte Carlo samples for marginal predictions (default: 1000)

- ...:

  Additional arguments (currently unused)

## Value

Depends on `type`
