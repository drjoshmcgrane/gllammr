# Reconstruct per-item threshold parameters from a fitted polytomous EIRT model

Mirrors the threshold construction in the TMB template
(gllamm_eirt.hpp): GRM uses ordered sum-to-zero deviations around the
item location, PCM/GPCM use sum-to-zero step deviations around the item
location, LPCM uses the threshold regression plus step residuals.

## Usage

``` r
eirt_item_thresholds(object)
```

## Arguments

- object:

  Fitted gllamm_eirt model (polytomous)

## Value

List of per-item threshold vectors on the absolute scale
