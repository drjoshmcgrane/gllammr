# GLLAMM model object class

S3 class for fitted GLLAMM models

## Usage

``` r
# S3 method for class 'gllamm'
print(x, ...)

# S3 method for class 'gllamm'
summary(object, ...)

# S3 method for class 'gllamm'
coef(object, ...)

# S3 method for class 'gllamm'
vcov(object, which = "fixed", type = c("model", "sandwich"), ...)

# S3 method for class 'gllamm'
logLik(object, ...)

# S3 method for class 'gllamm'
fitted(object, ...)

# S3 method for class 'gllamm'
residuals(object, type = c("response", "pearson", "deviance"), ...)
```

## Arguments

- x:

  A gllamm object

- ...:

  Additional arguments

- object:

  A gllamm object

- which:

  Which covariance block to return: "fixed" (default) or "all"

- type:

  Covariance type: "model" (default) or "sandwich" (cluster-robust)
