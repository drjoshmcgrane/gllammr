# Predict from a fitted mixed-response model

Conditional predictions (fixed effects plus estimated group effects,
loading 1 on every outcome) for each outcome of the joint model.

## Usage

``` r
# S3 method for class 'gllamm_mixed'
predict(object, newdata = NULL, type = c("response", "link"), ...)
```

## Arguments

- object:

  Fitted `gllamm_mixed` object

- newdata:

  Optional data frame; groups unseen at fit time get a group effect of
  zero (population-level prediction)

- type:

  `"response"` (default) or `"link"`

- ...:

  Unused

## Value

Named list with one prediction vector per outcome
