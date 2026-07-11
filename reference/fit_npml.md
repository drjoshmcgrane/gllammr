# Fit Two-Level GLMMs by Nonparametric Maximum Likelihood (NPML)

Replaces the normal random-intercept distribution with a discrete
distribution on `k` estimated mass points (locations and masses). The
marginal likelihood is an exact finite mixture - no integral
approximation - making NPML robust to misspecification of the
random-effects distribution.

## Usage

``` r
fit_npml(
  formula,
  data,
  k = 2,
  family = stats::gaussian(),
  weights = NULL,
  n_starts = 3,
  start = NULL,
  control = list()
)
```

## Arguments

- formula:

  Model formula `y ~ x + (1 | group)`

- data:

  Data frame

- k:

  Number of mass points (default 2)

- family:

  gaussian(), binomial(), or poisson()

- weights:

  Optional observation weights

- n_starts:

  Random restarts (mixtures are prone to local optima)

- start:

  Optional starting values

- control:

  Optimization control list

## Value

An object of class `gllamm_npml`

## Details

The fixed-effects design drops its intercept: the mass-point locations
play the role of component intercepts (the usual NPML identification, as
in npmlreg). Masses are softmax-parameterized.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_npml(y ~ x + (1 | group), data = d, k = 3,
                family = stats::binomial())
} # }
```
