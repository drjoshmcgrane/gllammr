# GLLAMMR Quick Reference

## Installation

```r
# From GitHub (development version)
remotes::install_github("yourusername/GLLAMMR")

# Load package
library(GLLAMMR)
```

## Basic Usage

### Simple Random Intercept

```r
fit <- gllamm(y ~ x + (1 | group), data = mydata)
summary(fit)
```

### Random Intercept and Slope

```r
fit <- gllamm(y ~ x + (x | group), data = mydata)
```

### Uncorrelated Random Effects

```r
fit <- gllamm(y ~ x + (x || group), data = mydata)
```

### Nested Random Effects

```r
fit <- gllamm(y ~ x + (1 | level1/level2), data = mydata)
```

### Crossed Random Effects

```r
fit <- gllamm(y ~ x + (1 | group1) + (1 | group2), data = mydata)
```

## Formula Syntax

| Syntax | Description |
|--------|-------------|
| `y ~ x` | Fixed effect only |
| `(1 \| g)` | Random intercept |
| `(x \| g)` | Random intercept + slope (correlated) |
| `(x \|\| g)` | Random intercept + slope (uncorrelated) |
| `(1 \| a/b)` | Nested: b within a |
| `(1 \| a) + (1 \| b)` | Crossed random effects |

## GLM Families

```r
# Gaussian (default)
gllamm(y ~ x + (1|g), data, family = gaussian())

# Binomial (coming in Phase 2)
gllamm(y ~ x + (1|g), data, family = binomial(link = "logit"))

# Poisson (coming in Phase 2)
gllamm(y ~ x + (1|g), data, family = poisson(link = "log"))
```

## Extracting Results

```r
# Fixed effects
fixef(fit)

# Random effects (BLUPs)
ranef(fit)

# Variance components
VarCorr(fit)

# All coefficients
coef(fit)

# Variance-covariance matrix
vcov(fit)

# Log-likelihood
logLik(fit)

# Information criteria
AIC(fit)
BIC(fit)

# Fitted values
fitted(fit)

# Residuals
residuals(fit, type = "response")  # or "pearson", "deviance"
```

## Prediction

```r
# Predictions with random effects (default)
predict(fit)

# Population-level (no random effects)
predict(fit, re.form = NA)

# New data (coming soon)
predict(fit, newdata = newdata)
```

## Simulation

```r
# Single simulation
sim <- simulate(fit, nsim = 1)

# Multiple simulations
sims <- simulate(fit, nsim = 100, seed = 123)
```

## Model Comparison

```r
fit1 <- gllamm(y ~ x + (1|g), data)
fit2 <- gllamm(y ~ x + z + (1|g), data)

# Compare AICs
AIC(fit1, fit2)

# Likelihood ratio test (if nested)
anova(fit1, fit2)  # Coming soon
```

## Diagnostics

```r
# Plot diagnostics (coming in Phase 8)
plot(fit)

# Check convergence
fit$convergence$converged
fit$convergence$message

# Residual plots
plot(fitted(fit), residuals(fit))
qqnorm(residuals(fit))
```

## Control Parameters

```r
fit <- gllamm(
  y ~ x + (1|g),
  data = mydata,
  control = list(
    eval.max = 1000,  # Max function evaluations
    iter.max = 500,   # Max iterations
    trace = 1         # Print optimization progress
  )
)
```

## Common Issues

### TMB Not Compiled

```r
# Compile TMB template
TMB::compile(system.file("src/gllamm_gaussian.cpp", package = "GLLAMMR"))
```

### Convergence Issues

```r
# Try different starting values
fit <- gllamm(y ~ x + (1|g), data, start = my_starts)

# Increase iterations
fit <- gllamm(y ~ x + (1|g), data,
              control = list(iter.max = 1000))
```

### Missing Variables

```r
# Check data
str(mydata)
summary(mydata)

# Ensure all variables in formula exist in data
```

## Getting Help

```r
# Function help
?gllamm
?predict.gllamm
?simulate.gllamm

# Package overview
help(package = "GLLAMMR")

# Vignettes
browseVignettes("GLLAMMR")
vignette("getting-started", package = "GLLAMMR")
```

## Example Datasets

```r
# lme4 sleepstudy
data(sleepstudy, package = "lme4")
fit <- gllamm(Reaction ~ Days + (Days|Subject), data = sleepstudy)
```

## Citation

```r
citation("GLLAMMR")
```

## Current Limitations (Phase 1)

- Only Gaussian family with identity link
- Only single random intercept term
- No random slopes yet (coming Week 2)
- No other GLM families (coming Phase 2)

## Upcoming Features

**Phase 2 (Weeks 5-8)**:
- Random slopes
- Binomial and Poisson families
- 3+ levels
- Adaptive quadrature

**Phase 4 (Weeks 11-14)**:
- IRT models (Rasch, 2PL, 3PL)
- Factor analysis

**Phase 5 (Weeks 15-17)**:
- Latent class models

See `ROADMAP.md` for complete development plan.
