# Binomial Family for Binary and Binomial Outcomes

Create a family object for binomial regression models with various link
functions

## Usage

``` r
binomial(link = c("logit", "probit", "cloglog"))
```

## Arguments

- link:

  Link function to use:

  "logit"

  :   Logistic regression (default) - symmetric S-curve

  "probit"

  :   Probit regression - based on normal distribution

  "cloglog"

  :   Complementary log-log - asymmetric, suitable for rare events and
      survival data

## Value

A family object of class `binomial_family` with components:

- family:

  Character: "binomial"

- link:

  Character: name of link function

- link_code:

  Integer: numeric code for TMB (1-3)

## Details

The binomial family is used for binary (0/1) or binomial (count/total)
responses.

**Logit Link (default):** \$\$P(Y=1\|x) = \frac{1}{1 +
\exp(-x'\beta)}\$\$ This is the standard logistic regression. The logit
link is symmetric and appropriate when the probability of success and
failure are equally likely to change with covariates.

**Probit Link:** \$\$P(Y=1\|x) = \Phi(x'\beta)\$\$ where \\\Phi\\ is the
standard normal CDF. This link is also symmetric and yields similar
results to logit but with slightly different tail behavior.

**Complementary Log-Log (cloglog) Link:** \$\$P(Y=1\|x) = 1 -
\exp(-\exp(x'\beta))\$\$ This link is *asymmetric* and is particularly
useful for:

- Rare events (when P(Y=1) is small)

- Survival analysis with discrete time intervals

- Gompertz or extreme value distributions

- When the hazard is proportional (as in survival models)

The cloglog link arises naturally when modeling grouped survival data or
when events follow a Poisson process over time intervals.

## See also

[`gllamm`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md),
[`ordinal`](https://drjoshmcgrane.github.io/gllammr/reference/ordinal.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Logistic regression (default)
family1 <- binomial()
family2 <- binomial(link = "logit")

# Probit regression
family3 <- binomial(link = "probit")

# Complementary log-log for rare events
family4 <- binomial(link = "cloglog")

# Use with gllamm() - recommended interface
fit <- gllamm(outcome ~ age + treatment + (1 | clinic),
              data = mydata,
              family = binomial(link = "logit"))

# Rare event with cloglog
fit_rare <- gllamm(rare_disease ~ exposure + (1 | region),
                   data = epi_data,
                   family = binomial(link = "cloglog"))
} # }
```
