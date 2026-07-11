# IRT Family for Item Response Theory Models

Create a family object for fitting IRT models through the unified
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
interface. The response is a persons x items matrix passed as the first
argument of
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md).

## Usage

``` r
irt(
  model = c("Rasch", "2PL", "3PL", "GRM", "PCM", "GPCM", "NRM"),
  mc_items = NULL,
  se = TRUE
)
```

## Arguments

- model:

  IRT model type: "Rasch", "2PL", "3PL" (dichotomous) or "GRM", "PCM",
  "GPCM", "NRM" (polytomous)

- mc_items:

  For 3PL only: which items have guessing parameters (NULL = all;
  logical or integer index vector)

- se:

  Compute standard errors (default TRUE); see
  [`fit_irt`](https://drjoshmcgrane.github.io/gllammr/reference/fit_irt.md)

## Value

A family object of class `irt_family`

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- gllamm(response_matrix, family = irt("2PL"))
# Multi-level IRT: persons nested in classes
fit_ml <- gllamm(response_matrix, data = person_data,
                 family = irt("Rasch"), random = ~ (1 | class))
} # }
```
