# Total abilities for multi-level IRT fits

Total abilities for multi-level IRT fits

## Usage

``` r
# S3 method for class 'gllamm_irt_multilevel'
abilities(object, composite = FALSE, ...)
```

## Arguments

- object:

  A fitted multi-level IRT model

- composite:

  If TRUE, return person deviation plus group effects (total ability);
  if FALSE (default), the person deviation only

- ...:

  Additional arguments (not used)

## Value

A named vector of person abilities
