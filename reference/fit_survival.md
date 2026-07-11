# Fit Parametric Survival Models with Random Effects (Frailty)

Fits exponential or Weibull proportional-hazards models with log-normal
frailty (normally distributed random effects on the log-hazard scale),
supporting right censoring.

## Usage

``` r
fit_survival(
  formula,
  data,
  distribution = c("weibull", "exponential"),
  weights = NULL,
  start = NULL,
  control = list()
)
```

## Arguments

- formula:

  Formula of the form `Surv(time, event) ~ x + (1 | group)`. The
  left-hand side names the time and event (1 = event, 0 = censored)
  variables; the survival package is not required.

- data:

  Data frame

- distribution:

  "weibull" (default) or "exponential"

- weights:

  Optional case weights

- start:

  Optional starting values

- control:

  Optimization control list

## Value

An object of class `gllamm_survival`

## Details

The cumulative hazard is \\H(t \mid x, u) = (\lambda t)^{shape}\\ with
\\\lambda = \exp(x'\beta + z'u)\\; the exponential model fixes shape
= 1. Under this parameterization the accelerated-failure-time
coefficients of
[`survival::survreg`](https://rdrr.io/pkg/survival/man/survreg.html)
correspond to \\-\beta\\ and its scale to \\1/shape\\.

The exponential frailty model is likelihood-equivalent to a Poisson GLMM
on the event indicator with offset \\\log(t)\\, which is used in the
package validation suite.

## Parameterization

The exponential model has hazard \\\exp(\eta)\\, so coefficients are log
hazard ratios. The Weibull model is parameterized as \\S(t) =
\exp(-(\exp(\eta)\\t)^{shape})\\: \\\eta\\ scales time
(accelerated-failure-time form), and the log hazard ratio for a
covariate is \\shape \times \beta\\. This matches
[`survival::survreg`](https://rdrr.io/pkg/survival/man/survreg.html)
with \\\beta = -\beta\_{AFT}\\ and \\shape = 1/scale\\.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_survival(Surv(time, status) ~ age + (1 | center),
                    data = d, distribution = "weibull")
} # }
```
