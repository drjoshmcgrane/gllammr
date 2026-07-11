# Ordinal Family for Proportional and Non-Proportional Odds Models

Create a family object for ordinal regression models with various link
functions

## Usage

``` r
ordinal(
  link = c("logit", "probit", "acl", "crl_forward", "crl_backward", "ppo")
)
```

## Arguments

- link:

  Link function to use:

  "logit"

  :   Proportional odds (cumulative logit) - default

  "probit"

  :   Cumulative probit

  "acl"

  :   Adjacent category logit

  "crl_forward"

  :   Forward continuation ratio logit

  "crl_backward"

  :   Backward continuation ratio logit

  "ppo"

  :   Partial proportional odds (non-proportional)

## Value

A family object of class `ordinal_family` with components:

- family:

  Character: "ordinal"

- link:

  Character: name of link function

- link_code:

  Integer: numeric code for TMB (1-6)

## Details

The ordinal family supports several models for ordered categorical
responses:

**Proportional Odds (logit):** \$\$P(Y \le k \| x) = \frac{1}{1 +
\exp(-(\tau_k - x'\beta))}\$\$

**Cumulative Probit:** \$\$P(Y \le k \| x) = \Phi(\tau_k - x'\beta)\$\$

**Adjacent Category Logit (ACL):** Models the log-odds of adjacent
categories: \$\$\log\frac{P(Y=k)}{P(Y=k-1)} = \alpha_k + x'\beta\$\$

**Continuation Ratio Logit (CRL):** Forward version models sequential
decisions: \$\$\log\frac{P(Y=k \| Y \ge k)}{P(Y\>k \| Y \ge k)} =
\tau_k - x'\beta\$\$

Backward version reverses the conditioning.

**Partial Proportional Odds (PPO):** Relaxes the proportional odds
assumption by allowing different covariate effects per threshold:
\$\$P(Y \le k \| x) = F(\tau_k - x'\beta_k)\$\$

## Examples

``` r
if (FALSE) { # \dontrun{
# Proportional odds model (default)
family1 <- ordinal()
family2 <- ordinal(link = "logit")

# Adjacent category logit
family3 <- ordinal(link = "acl")

# Partial proportional odds
family4 <- ordinal(link = "ppo")

# Use with gllamm() - recommended interface
fit <- gllamm(rating ~ temp + (1 | judge),
              data = wine,
              family = ordinal(link = "logit"))

# Or use fit_ordinal() directly
fit2 <- fit_ordinal(rating ~ temp + (1 | judge),
                    data = wine,
                    link = "acl")
} # }
```
