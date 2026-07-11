# Predict from a fitted NPML model

Conditional predictions for NPML fits. Each training group gets its
posterior-mean intercept over the estimated mass points; groups unseen
at fit time (or absent from `newdata`) get the prior mean
`sum(masses * locations)`.

## Usage

``` r
# S3 method for class 'gllamm_npml'
predict(object, newdata = NULL, type = c("response", "link"), ...)
```

## Arguments

- object:

  Fitted `gllamm_npml` object

- newdata:

  Optional data frame (defaults to the training data)

- type:

  `"response"` (default) or `"link"`

- ...:

  Unused

## Value

Numeric vector of predictions
