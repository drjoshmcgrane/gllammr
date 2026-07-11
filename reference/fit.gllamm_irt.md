# Fit Statistics for IRT Models

Fit Statistics for IRT Models

## Usage

``` r
# S3 method for class 'gllamm_irt'
fit(object, compute_item_fit = TRUE, compute_person_fit = TRUE, ...)
```

## Arguments

- object:

  A gllamm_irt object

- compute_item_fit:

  Compute item fit statistics (S-X^2) (default: TRUE)

- compute_person_fit:

  Compute person fit statistics (outfit/infit) (default: TRUE)

- ...:

  Additional arguments

## Value

Object of class `fit_irt` with IRT-specific fit statistics
