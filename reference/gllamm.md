# Fit Generalized Linear Latent and Mixed Models

Main function for fitting GLLAMM models. Supports multilevel generalized
linear models, factor models, IRT models, latent class models, and more.

## Usage

``` r
gllamm(
  formula,
  data,
  family = gaussian(),
  weights = NULL,
  random = NULL,
  integration = NULL,
  start = NULL,
  control = list(),
  ...
)
```

## Arguments

- formula:

  A two-sided formula with syntax: `y ~ x + (terms | group)`. The
  right-hand side contains fixed effects and random effects
  specifications. Random effects are specified using `(term | group)`
  for correlated random effects or `(term || group)` for uncorrelated
  random effects. Nested random effects can be specified as
  `(1 | level1/level2)`.

- data:

  A data frame containing the variables in the formula.

- family:

  A family object selecting the model class. Every model in the package
  is reachable here:

  - [`gaussian()`](https://rdrr.io/r/stats/family.html) (default),
    [`binomial()`](https://drjoshmcgrane.github.io/gllammr/reference/binomial.md),
    [`poisson()`](https://rdrr.io/r/stats/family.html),
    [`Gamma()`](https://rdrr.io/r/stats/family.html) - GLMMs with any
    number of crossed/nested random-effects terms and random slopes

  - `ordinal(link)` - cumulative, adjacent-category, continuation-ratio,
    and partial-proportional-odds models

  - `multinomial(reference)` - baseline-category logit

  - `irt(model)` / `eirt(item_data, ...)` - (explanatory) item response
    models; first argument is the response matrix

  - `lca(nclass, ordering)` / `cdm(Q, model)` - latent class and
    cognitive diagnosis models; first argument is the response matrix

  - `sem(measurement, structural)` - structural equation models; first
    argument is the data frame

  - `mixed_response(...)` - joint mixed-type outcomes; first argument is
    the shared random-effects formula

  - `survival_family(distribution)` - parametric frailty survival with
    `Surv(time, event)` on the left-hand side

  - `ranking(case)` - rank-ordered (exploded) logit

  The latent distribution is normal with Laplace integration by default;
  `integration = aghq(k)` requests adaptive quadrature and
  `integration = npml(k)` a nonparametric (mass-point) latent
  distribution.

- weights:

  Optional case weights: a numeric vector of observation (level-1)
  weights, or a list with elements `level1` and/or `level2` for survey
  designs with weights at both levels. Under the default Laplace fit,
  level-2 weights must be integer frequency weights; they are
  implemented by exact replication of whole groups, so results are
  identical to fitting the duplicated data. Non-integer level-2 weights
  require `integration = aghq(k)`, which weights each group's log
  marginal likelihood directly and supports arbitrary weights.

- random:

  For matrix-response families
  ([`irt()`](https://drjoshmcgrane.github.io/gllammr/reference/irt.md),
  [`lca()`](https://drjoshmcgrane.github.io/gllammr/reference/lca.md)):
  an optional person-level random-effects formula such as
  `~ (1 | class)`.

- integration:

  Optional integration specification; `aghq(k)` selects adaptive
  Gauss-Hermite quadrature with `k` nodes for two-level random-intercept
  models (default is the Laplace approximation).

- start:

  Optional named list of starting values for parameters.

- control:

  A list of control parameters for the optimization algorithm:

  eval.max

  :   Maximum number of function evaluations (default: 1000)

  iter.max

  :   Maximum number of iterations (default: 500)

  trace

  :   Integer controlling printed output (default: 0)

- ...:

  Additional arguments (reserved for future use).

## Value

An object of class `gllamm` with components:

- coefficients:

  List with `fixed` (fixed effects coefficients) and `random_var`
  (random effects variance components)

- vcov:

  List with `fixed` (variance-covariance matrix of fixed effects) and
  `all` (full variance-covariance matrix)

- random_effects:

  List of random effects predictions by group

- fitted.values:

  Fitted values on the response scale

- residuals:

  Response residuals

- y:

  Response vector

- X:

  Fixed effects design matrix

- logLik:

  Log-likelihood at convergence

- AIC:

  Akaike Information Criterion

- BIC:

  Bayesian Information Criterion

- n_obs:

  Number of observations

- n_params:

  Number of parameters

- n_groups:

  Number of groups for each random effects term

- convergence:

  List with convergence information

- call:

  The matched call

- formula:

  The model formula

- family:

  The GLM family

- data:

  The data frame (if requested)

- random_terms:

  List of random effects specifications

## Details

gllammr uses Template Model Builder (TMB) for efficient computation via
automatic differentiation. The random effects are integrated out using
Laplace approximation, providing fast and accurate inference.

The formula syntax follows lme4 conventions:

- `(1 | group)` - Random intercept for group

- `(x | group)` - Random intercept and slope for x

- `(x || group)` - Uncorrelated random intercept and slope

- `(1 | level1/level2)` - Nested random effects

- `(1 | group1) + (1 | group2)` - Crossed random effects

## Missing data

Formula-based models (GLMMs, ordinal, multinomial, survival, NPML, mixed
responses) listwise-delete rows with missing values in any model
variable, with a warning; supplied weights are aligned automatically.
Matrix-response latent variable models
([`irt()`](https://drjoshmcgrane.github.io/gllammr/reference/irt.md),
[`lca()`](https://drjoshmcgrane.github.io/gllammr/reference/lca.md),
[`cdm()`](https://drjoshmcgrane.github.io/gllammr/reference/cdm.md), and
the corresponding `fit_*` functions) use all observed responses -
item-level missingness is handled under the MAR assumption by the
marginal likelihood itself. `fit_sem` offers full-information maximum
likelihood via `missing = "fiml"`. `fit_rank` treats missing ranks as
deliberately unranked alternatives (partial rankings), not as missing
data.

## References

Rabe-Hesketh, S., Skrondal, A., & Pickles, A. (2004). GLLAMM Manual.
U.C. Berkeley Division of Biostatistics Working Paper Series.

Skrondal, A., & Rabe-Hesketh, S. (2004). Generalized Latent Variable
Modeling: Multilevel, Longitudinal, and Structural Equation Models.
Chapman & Hall/CRC.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic random intercept model
data(sleepstudy, package = "lme4")
fit1 <- gllamm(Reaction ~ Days + (1 | Subject),
               data = sleepstudy)
summary(fit1)

# Random intercept and slope
fit2 <- gllamm(Reaction ~ Days + (Days | Subject),
               data = sleepstudy)
summary(fit2)

# Extract components
fixef(fit2)        # Fixed effects
ranef(fit2)        # Random effects
VarCorr(fit2)      # Variance components
fitted(fit2)       # Fitted values
residuals(fit2)    # Residuals

# Ordinal regression (proportional odds)
data$satisfaction <- ordered(sample(1:5, nrow(data), replace = TRUE))
fit3 <- gllamm(satisfaction ~ x + (1 | group),
               data = data,
               family = ordinal(link = "logit"))

# Ordinal with adjacent category logit
fit4 <- gllamm(satisfaction ~ x + (1 | group),
               data = data,
               family = ordinal(link = "acl"))
} # }
```
