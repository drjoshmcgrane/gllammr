# Simulate response matrices from a fitted explanatory IRT model

Draws fresh person abilities (plus group effects when present is not
supported - simulation is at the population level) and samples responses
from the estimated item parameters.

## Usage

``` r
# S3 method for class 'gllamm_eirt'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  Fitted `gllamm_eirt` object

- nsim:

  Number of replicates

- seed:

  Optional seed (stored as the `"seed"` attribute)

- ...:

  Unused
