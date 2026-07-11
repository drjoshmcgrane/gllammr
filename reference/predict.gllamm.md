# Predict method for GLLAMM models

Obtain predictions from a fitted GLLAMM model

## Usage

``` r
# S3 method for class 'gllamm'
predict(
  object,
  newdata = NULL,
  type = c("response", "link", "random", "marginal"),
  re.form = NULL,
  n_sim = 1000,
  se.fit = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted `gllamm` object

- newdata:

  Optional new data frame for predictions. If omitted, fitted values
  from the original data are returned.

- type:

  Type of prediction:

  response

  :   Predictions on the response scale (default)

  link

  :   Predictions on the link scale (linear predictor)

  random

  :   Random effects only

  marginal

  :   Population-averaged predictions (integrating over random effects)

- re.form:

  Formula for random effects to include. Use `NA` or `~0` to exclude all
  random effects (population-level predictions). Ignored when
  `type = "marginal"`.

- n_sim:

  Number of Monte Carlo samples for marginal predictions (default:
  1000). Only used when `type = "marginal"`.

- se.fit:

  Logical; return standard errors for marginal predictions? (default:
  FALSE). Only used when `type = "marginal"`.

- ...:

  Additional arguments (currently unused)

## Value

A vector of predictions

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- gllamm(y ~ x + (1 | group), data = mydata)

# Fitted values (default - conditional on random effects)
pred1 <- predict(fit)

# Population-level predictions (fixed effects only, u=0)
pred2 <- predict(fit, re.form = NA)

# Marginal predictions (population-averaged, integrating over u)
pred3 <- predict(fit, type = "marginal")

# Marginal predictions with standard errors
pred4 <- predict(fit, type = "marginal", se.fit = TRUE)

# Marginal predictions for new data
pred5 <- predict(fit, newdata = newdata, type = "marginal")
} # }
```
