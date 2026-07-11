# Fit Ordinal Regression Models with Random Effects

Fit proportional odds or cumulative probit models for ordinal responses.
This function can be called directly or through
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
with `family = ordinal()`. The
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
interface is recommended for consistency with other model types.

## Usage

``` r
fit_ordinal(
  formula,
  data,
  link = c("logit", "probit", "acl", "crl_forward", "crl_backward", "ppo"),
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

  Link function: "logit" (proportional odds) or "probit"

- weights:

  Optional vector of case weights (one per observation)

- start:

  Optional starting values

- control:

  Control parameters

## Value

An object of class `gllamm_ordinal`

## Details

For ordinal response Y with K categories (1, 2, ..., K), the model is:

Proportional odds (logit link): \$\$P(Y \le k \| x) = \frac{1}{1 +
\exp(-(\tau_k - x'\beta))}\$\$

Cumulative probit: \$\$P(Y \le k \| x) = \Phi(\tau_k - x'\beta)\$\$

where \\\tau_1 \< \tau_2 \< ... \< \tau\_{K-1}\\ are threshold
parameters.

## Note

The recommended interface is
`gllamm(formula, data, family = ordinal(link))`. This function is also
available for direct use with the `link` argument.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate ordinal data
data$satisfaction <- factor(sample(1:5, 100, replace = TRUE),
                            ordered = TRUE,
                            levels = 1:5)

# Recommended: Use gllamm() with ordinal() family
fit1 <- gllamm(satisfaction ~ age + (1 | clinic),
               data = data,
               family = ordinal(link = "logit"))
summary(fit1)

# Alternative: Call fit_ordinal() directly
fit2 <- fit_ordinal(satisfaction ~ age + (1 | clinic),
                    data = data,
                    link = "logit")
summary(fit2)
} # }
```
