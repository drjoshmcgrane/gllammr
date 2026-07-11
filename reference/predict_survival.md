# Predict method for survival models

Obtain predictions from a fitted survival model with random effects

## Usage

``` r
# S3 method for class 'gllamm_survival'
predict(
  object,
  newdata = NULL,
  type = c("lp", "risk", "survival", "hazard", "marginal_survival", "marginal_hazard"),
  times = NULL,
  n_sim = 1000,
  ...
)
```

## Arguments

- object:

  A fitted survival model

- newdata:

  Optional new data frame for predictions

- type:

  Type of prediction:

  lp

  :   Linear predictor (η)

  risk

  :   Relative risk exp(η)

  survival

  :   Survival probability at specified times

  hazard

  :   Hazard at specified times

  marginal_survival

  :   Marginal survival probability (population-averaged)

  marginal_hazard

  :   Marginal hazard (population-averaged)

- times:

  Times at which to evaluate survival/hazard (required for
  survival/hazard types)

- n_sim:

  Number of Monte Carlo samples for marginal predictions (default: 1000)

- ...:

  Additional arguments (currently unused)

## Value

Depends on `type`
