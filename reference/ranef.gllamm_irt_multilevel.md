# Extract Random Effects from Multi-Level IRT Models

Extract estimated random effects (group-level deviations) from fitted
models

## Usage

``` r
# S3 method for class 'gllamm_irt_multilevel'
ranef(object, level = NULL, ...)
```

## Arguments

- object:

  A fitted multi-level IRT model

- level:

  Which random effect level to extract. If NULL, returns all levels.

- ...:

  Additional arguments (not used)

## Value

A named vector or matrix of random effects

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract class effects
ranef(fit, level = "class")

# Extract all random effects
ranef(fit)
} # }
```
