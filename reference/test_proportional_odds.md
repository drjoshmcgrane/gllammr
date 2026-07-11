# Test Proportional Odds Assumption

Perform a likelihood ratio test of the proportional odds assumption by
comparing a proportional odds model to a partial proportional odds model

## Usage

``` r
test_proportional_odds(object, data = NULL)
```

## Arguments

- object:

  A fitted ordinal regression model (gllamm_ordinal)

- data:

  The data frame used to fit the model (required when the fitted object
  does not store its data)

## Value

An object of class `po_test` with components:

- statistic:

  Likelihood ratio test statistic

- df:

  Degrees of freedom for the test

- p_value:

  P-value from chi-squared distribution

- conclusion:

  Text interpretation of the test result

- models:

  List containing the base and PPO models

## Details

The proportional odds assumption states that the effect of covariates is
the same across all thresholds. This function fits a partial
proportional odds (PPO) model where each threshold can have different
covariate effects and tests whether this provides a significantly better
fit.

The test statistic is: \$\$LRT = 2(logLik\_{PPO} - logLik\_{PO})\$\$

which follows a chi-squared distribution with degrees of freedom equal
to the difference in number of parameters.

If p \< 0.05, the proportional odds assumption is rejected, suggesting
that covariate effects vary across thresholds.

## Note

This function currently only works for models with logit or probit
links. It will not work with ACL, CRL, or already-PPO models.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit proportional odds model
fit_po <- fit_ordinal(rating ~ temp + (1 | judge),
                      data = wine, link = "logit")

# Test proportional odds assumption
po_test <- test_proportional_odds(fit_po)
print(po_test)
} # }
```
