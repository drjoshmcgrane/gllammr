# Predicted category probabilities for polytomous IRT fits

Returns the model-implied category probabilities at each person's
estimated ability: a list (one element per item) of persons x categories
matrices, or the persons x items matrix of expected scores with
`type = "expected"`.

## Usage

``` r
# S3 method for class 'gllamm_irt_poly'
predict(
  object,
  type = c("probs", "expected", "marginal", "probability", "ability"),
  ...
)
```

## Arguments

- object:

  Fitted polytomous `gllamm_irt` model

- type:

  "probs" (default) or "expected"

- ...:

  Unused
