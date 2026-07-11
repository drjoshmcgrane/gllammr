# Parametric Frailty Survival Family

Create a family object for parametric survival models with shared
(log-normal) frailties through the unified
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
interface; the formula uses `Surv(time, event)` on the left-hand side.
See
[`fit_survival`](https://drjoshmcgrane.github.io/gllammr/reference/fit_survival.md).

## Usage

``` r
survival_family(distribution = c("exponential", "weibull"))
```

## Arguments

- distribution:

  "exponential" (default) or "weibull"

## Value

An object of class `survival_family`

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- gllamm(Surv(time, status) ~ x + (1 | clinic), data = d,
              family = survival_family("weibull"))
} # }
```
