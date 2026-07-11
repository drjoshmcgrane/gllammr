# Compare Explanatory IRT Models

Compare two or more EIRT models using likelihood ratio tests and
information criteria

## Usage

``` r
compare_eirt(..., test = c("LRT", "none"))
```

## Arguments

- ...:

  Two or more gllamm_eirt objects to compare

- test:

  Type of test: "LRT" for likelihood ratio test, "none" for just IC
  comparison

## Value

A data frame with model comparison statistics

## Details

This function compares EIRT models to test the importance of item
covariates. Models are compared using:

- Log-likelihood

- AIC and BIC

- Likelihood ratio test (if nested models)

For nested models (e.g., with and without a predictor), the LRT tests
whether adding the predictor significantly improves fit.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit model with predictor
fit1 <- fit_eirt(responses, item_data,
                 difficulty_formula = ~ word_freq,
                 model = "Rasch")

# Fit model without predictor
fit0 <- fit_eirt(responses, item_data,
                 difficulty_formula = ~ 1,
                 model = "Rasch")

# Compare models
compare_eirt(fit0, fit1)
} # }
```
