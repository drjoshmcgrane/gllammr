# Fit Binomial Regression Models with Random Effects

Fit logistic, probit, or complementary log-log models for binary or
binomial responses. This function can be called directly or through
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
with `family = binomial()`. The
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
interface is recommended for consistency with other model types.

## Usage

``` r
fit_binomial(
  formula,
  data,
  link = c("logit", "probit", "cloglog"),
  weights = NULL,
  start = NULL,
  control = list()
)
```

## Arguments

- formula:

  Formula with syntax: y ~ x + (terms \| group)

- data:

  Data frame

- link:

  Link function: "logit" (default), "probit", or "cloglog"

- weights:

  Optional vector of case weights (one per observation)

- start:

  Optional starting values

- control:

  Control parameters for optimization

## Value

An object of class `gllamm_binomial` with components:

- coefficients:

  List of fixed effects and random effect variances

- logLik:

  Log-likelihood at convergence

- AIC:

  Akaike Information Criterion

- BIC:

  Bayesian Information Criterion

- convergence:

  Convergence information

- link:

  Link function used

- n_obs:

  Number of observations

- fitted_values:

  Fitted probabilities

- tmb_obj:

  TMB object for further inference

- tmb_opt:

  Optimization result

- tmb_sdr:

  Standard errors via sdreport

## Details

For binary response Y in {0, 1} or binomial response Y/n, the model is:

**Logit link (default - logistic regression):** \$\$P(Y=1\|x) =
\frac{1}{1 + \exp(-x'\beta - Z'u)}\$\$

**Probit link:** \$\$P(Y=1\|x) = \Phi(x'\beta + Z'u)\$\$ where \\\Phi\\
is the standard normal CDF.

**Complementary log-log (cloglog) link:** \$\$P(Y=1\|x) = 1 -
\exp(-\exp(x'\beta + Z'u))\$\$

The cloglog link is asymmetric and particularly useful for:

- Rare events (when baseline probability is low)

- Survival analysis with discrete time

- When modeling hazards that are proportional

- Grouped survival data from an underlying Poisson process

Random effects u are assumed multivariate normal with mean 0.

## Note

The recommended interface is
`gllamm(formula, data, family = binomial(link))`. This function is also
available for direct use with the `link` argument.

## See also

[`gllamm`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md),
[`binomial`](https://drjoshmcgrane.github.io/gllammr/reference/binomial.md),
[`ordinal`](https://drjoshmcgrane.github.io/gllammr/reference/ordinal.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate binary data
set.seed(123)
n_groups <- 20
n_per_group <- 10
data <- data.frame(
  group = rep(1:n_groups, each = n_per_group),
  x = rnorm(n_groups * n_per_group),
  y = rbinom(n_groups * n_per_group, 1, 0.5)
)

# Recommended: Use gllamm() with binomial() family
fit1 <- gllamm(y ~ x + (1 | group),
               data = data,
               family = binomial(link = "logit"))
summary(fit1)

# Probit link
fit2 <- gllamm(y ~ x + (1 | group),
               data = data,
               family = binomial(link = "probit"))

# Complementary log-log for rare events
data$rare_event <- rbinom(nrow(data), 1, 0.05)
fit3 <- gllamm(rare_event ~ x + (1 | group),
               data = data,
               family = binomial(link = "cloglog"))
summary(fit3)

# Alternative: Call fit_binomial() directly
fit4 <- fit_binomial(y ~ x + (1 | group),
                     data = data,
                     link = "logit")
} # }
```
