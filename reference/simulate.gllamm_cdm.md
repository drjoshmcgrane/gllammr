# Simulate response matrices from a fitted cognitive diagnosis model

Draws attribute profiles from the estimated prevalences and item
responses from the per-profile item kernels.

## Usage

``` r
# S3 method for class 'gllamm_cdm'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  Fitted `gllamm_cdm` object

- nsim:

  Number of replicates

- seed:

  Optional seed (stored as the `"seed"` attribute)

- ...:

  Unused
