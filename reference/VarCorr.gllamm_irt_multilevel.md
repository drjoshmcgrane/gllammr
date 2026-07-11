# Extract Variance Components from Multi-Level IRT Models

Extract variance components from fitted multi-level IRT models

## Usage

``` r
# S3 method for class 'gllamm_irt_multilevel'
VarCorr(x, ...)
```

## Arguments

- x:

  A fitted multi-level IRT model (class gllamm_irt_multilevel)

- ...:

  Additional arguments (not used)

## Value

A data frame with variance components

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit multi-level model
fit <- fit_irt(responses, model = "2PL",
               person_data = data, random = ~ (1 | class))

# Extract variance components
VarCorr(fit)
} # }
```
