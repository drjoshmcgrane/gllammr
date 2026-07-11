# Simulate from a GLLAMM model

Simulate response data from a fitted GLLAMM model. Random effects are
drawn fresh from their estimated distribution for every replicate
(population-level simulation), both for the original data and for
`newdata`.

## Usage

``` r
# S3 method for class 'gllamm'
simulate(object, nsim = 1, seed = NULL, newdata = NULL, ...)
```

## Arguments

- object:

  A fitted `gllamm` object

- nsim:

  Number of simulations (default: 1)

- seed:

  Optional random seed (stored in the `"seed"` attribute, following the
  [`simulate`](https://rdrr.io/r/stats/simulate.html) contract)

- newdata:

  Optional new data frame containing the covariates and grouping
  variables of the model formula

- ...:

  Additional arguments (currently unused)

## Value

A data frame with `nsim` columns, one simulated response vector per
column, with a `"seed"` attribute.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- gllamm(y ~ x + (1 | group), data = mydata)

# Single simulation
sim1 <- simulate(fit)

# Multiple simulations
sim10 <- simulate(fit, nsim = 10)
} # }
```
