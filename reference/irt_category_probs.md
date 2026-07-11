# Category response probabilities for polytomous IRT models

Single source of truth for the category probability math, mirroring the
likelihood in the TMB template (gllamm_irt_poly.hpp / gllamm_eirt.hpp).
Used by plotting, DIF displays, and marginal predictions so the formulas
cannot drift apart.

## Usage

``` r
irt_category_probs(model, theta, thresholds, discrimination = 1)
```

## Arguments

- model:

  One of "GRM", "PCM", "GPCM", "NRM"

- theta:

  Numeric vector of ability values

- thresholds:

  Numeric vector of item parameters: ordered thresholds (GRM), free step
  difficulties (PCM/GPCM), or category intercepts (NRM, reference
  category omitted)

- discrimination:

  Item discrimination (ignored for PCM)

## Value

Matrix \[length(theta) x K\] of category probabilities; rows sum to 1
