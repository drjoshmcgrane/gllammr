# Cognitive Diagnosis Family for Q-Matrix Models

Create a family object for fitting cognitive diagnosis models through
the unified
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
interface. The response is a persons x items binary matrix passed as the
first argument of
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md);
the Q-matrix and model options are carried by the family object. See
[`fit_cdm`](https://drjoshmcgrane.github.io/gllammr/reference/fit_cdm.md)
for details of the models and arguments.

## Usage

``` r
cdm(Q, model = c("gdina", "dina", "dino"), hierarchy = NULL, monotone = TRUE)
```

## Arguments

- Q:

  Binary Q-matrix (items x attributes)

- model:

  "gdina" (default), "dina", or "dino"

- hierarchy:

  Optional attribute hierarchy (list of prerequisite pairs)

- monotone:

  Enforce monotonicity in the attributes (default TRUE)

## Value

An object of class `cdm_family`

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- gllamm(Y, family = cdm(Q, model = "dina"))
} # }
```
