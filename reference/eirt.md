# Explanatory IRT Family

Create a family object for explanatory item response models through the
unified
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
interface. The response is the persons x items matrix passed as the
first argument of
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md);
`data` (optional) carries person-level variables and `random`
person-level random effects. See
[`fit_eirt`](https://drjoshmcgrane.github.io/gllammr/reference/fit_eirt.md).

## Usage

``` r
eirt(
  item_data,
  difficulty_formula = ~1,
  discrimination_formula = ~1,
  threshold_formula = NULL,
  step_formula = NULL,
  step_data = NULL,
  model = c("Rasch", "2PL", "GRM", "PCM", "GPCM"),
  item_residuals = TRUE
)
```

## Arguments

- item_data:

  Data frame of item covariates (one row per item)

- difficulty_formula:

  Item-covariate formula for difficulty

- discrimination_formula:

  Item-covariate formula for (log) discrimination (2PL/GRM/GPCM)

- threshold_formula:

  Optional threshold regression (LPCM framework)

- step_formula:

  Optional step-level covariate formula (PCM/GPCM)

- step_data:

  Data frame for step_formula (one row per item-step cell)

- model:

  "Rasch", "2PL", "GRM", "PCM", or "GPCM"

- item_residuals:

  Random item residuals around the regression (LLTM-plus-error; default
  TRUE)

## Value

An object of class `eirt_family`

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- gllamm(resp, family = eirt(item_data,
                                  difficulty_formula = ~ btype + mode))
} # }
```
