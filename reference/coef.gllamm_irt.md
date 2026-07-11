# Extract Coefficients from IRT Models

Extract item parameters from fitted IRT models

## Usage

``` r
# S3 method for class 'gllamm_irt'
coef(object, type = c("item", "person"), ...)
```

## Arguments

- object:

  A fitted IRT model

- type:

  Type of coefficients: "item" for item parameters, "person" for person
  abilities, "random" for random effects (multi-level only)

- ...:

  Additional arguments passed to specific methods

## Value

Item parameters (data frame), person abilities (vector), or random
effects (list)

## Examples

``` r
if (FALSE) { # \dontrun{
# Item parameters
coef(fit, type = "item")

# Person abilities
coef(fit, type = "person")

# Random effects (multi-level only)
coef(fit, type = "random")
} # }
```
