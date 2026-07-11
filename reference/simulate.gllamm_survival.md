# Simulate from a fitted parametric frailty survival model

Draws fresh frailties per replicate and generates uncensored event times
by inverse transform from the fitted exponential or Weibull hazard
(censoring schemes are design choices, so simulated times are
uncensored).

## Usage

``` r
# S3 method for class 'gllamm_survival'
simulate(object, nsim = 1, seed = NULL, newdata = NULL, ...)
```

## Arguments

- object:

  Fitted `gllamm_survival` object

- nsim:

  Number of replicates

- seed:

  Optional seed (stored as the `"seed"` attribute)

- newdata:

  Optional new data frame

- ...:

  Unused
