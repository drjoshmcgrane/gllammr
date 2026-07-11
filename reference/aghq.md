# Adaptive quadrature integration specification

Use as `gllamm(..., integration = aghq(k))` to integrate the random
intercept by adaptive Gauss-Hermite quadrature with `k` nodes instead of
the Laplace approximation. Currently supports two-level random-intercept
models with gaussian, binomial, or poisson families. Laplace (the
default) is equivalent to `aghq(1)`.

## Usage

``` r
aghq(k = 15)
```

## Arguments

- k:

  Number of quadrature nodes (default 15)

## Value

An object of class `gllamm_integration`
