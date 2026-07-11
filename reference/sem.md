# SEM Family for Structural Equation Models

Create a family object for structural equation models through the
unified
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
interface. The data frame is passed as the first argument of
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md);
the measurement and structural models live in the family object. See
[`fit_sem`](https://drjoshmcgrane.github.io/gllammr/reference/fit_sem.md).

## Usage

``` r
sem(measurement, structural = NULL, missing = c("listwise", "fiml"), se = TRUE)
```

## Arguments

- measurement:

  Named list of one-sided indicator formulas

- structural:

  Optional list of structural formulas (latent and/or observed
  predictors)

- missing:

  "listwise" (default) or "fiml"

- se:

  Compute standard errors (default TRUE)

## Value

An object of class `sem_family`

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- gllamm(d, family = sem(
  measurement = list(f1 = ~ x1 + x2 + x3, f2 = ~ y1 + y2 + y3),
  structural = list(f2 ~ f1)))
} # }
```
