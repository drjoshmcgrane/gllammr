# GLLAMMR User Guide

**Version 0.5.0** - Comprehensive Latent Variable Modeling in R

---

## What is GLLAMMR?

GLLAMMR (Generalized Linear Latent and Mixed Models in R) is a unified framework for fitting a wide range of statistical models involving latent variables and hierarchical structures. It provides:

- **Generalized Linear Mixed Models (GLMMs)** for clustered/nested data
- **Item Response Theory (IRT)** for educational and psychological testing
- **Latent Class Analysis (LCA)** for identifying hidden subgroups

All models use the TMB (Template Model Builder) backend for fast, efficient computation.

---

## Installation

```r
# Install dependencies
install.packages(c("TMB", "Matrix"))

# Load package (after compilation)
library(GLLAMMR)
```

### Compiling TMB Templates

Before first use, compile the C++ templates:

```r
# From package directory
TMB::compile("src/gllamm_gaussian.cpp")
TMB::compile("src/gllamm_gaussian_slopes.cpp")
TMB::compile("src/gllamm_binomial.cpp")
TMB::compile("src/gllamm_poisson.cpp")
TMB::compile("src/gllamm_irt.cpp")
TMB::compile("src/gllamm_latent_class.cpp")
```

Or use the helper:

```r
compile_gllamm_tmb()
```

---

## Quick Start Guide

### 1. Generalized Linear Mixed Models

#### Gaussian (Linear) Mixed Model

```r
# Random intercept model
fit1 <- gllamm(score ~ time + treatment + (1 | student),
               data = mydata,
               family = gaussian())

summary(fit1)

# Extract components
fixef(fit1)     # Fixed effects
ranef(fit1)     # Random effects (BLUPs)
VarCorr(fit1)   # Variance components
```

#### Random Slopes

```r
# Random intercept and slope for time
fit2 <- gllamm(score ~ time + treatment + (time | student),
               data = mydata)

summary(fit2)
```

#### Uncorrelated Random Effects

```r
# Random intercept and slope (uncorrelated)
fit3 <- gllamm(score ~ time + treatment + (time || student),
               data = mydata)
```

#### Binomial (Logistic) Regression

```r
# Binary outcome
fit4 <- gllamm(passed ~ hours_studied + (1 | school),
               data = mydata,
               family = binomial(link = "logit"))

summary(fit4)

# Fitted probabilities
fitted_probs <- fitted(fit4)
```

#### Poisson (Count) Regression

```r
# Count outcome
fit5 <- gllamm(num_arrests ~ poverty + (1 | neighborhood),
               data = mydata,
               family = poisson(link = "log"))

summary(fit5)
```

---

### 2. Item Response Theory Models

#### Rasch Model (1-Parameter Logistic)

The Rasch model assumes all items have equal discrimination:

$$P(Y_{ij} = 1) = \frac{1}{1 + \exp(-(θ_i - b_j))}$$

```r
# Create response matrix (persons × items)
# Each cell is 0 (incorrect) or 1 (correct)
responses <- matrix(c(
  1, 0, 1, 1, 0,
  0, 0, 0, 1, 0,
  1, 1, 1, 1, 1,
  # ... more rows
), ncol = 5, byrow = TRUE)

colnames(responses) <- paste0("Item", 1:5)

# Fit Rasch model
fit_rasch <- fit_irt(responses, model = "Rasch")

print(fit_rasch)
summary(fit_rasch)

# Item difficulties
difficulties <- fit_rasch$item_parameters$difficulty

# Person abilities
abilities <- fit_rasch$person_abilities

# Plot item difficulties
plot(difficulties, main = "Item Difficulty Parameters")
abline(h = 0, lty = 2, col = "red")
```

#### 2-Parameter Logistic (2PL) Model

Adds discrimination parameters:

```r
fit_2pl <- fit_irt(responses, model = "2PL")

# Item parameters
params_2pl <- fit_2pl$item_parameters
print(params_2pl)

# Higher discrimination = better differentiation
# between high and low ability students
plot(params_2pl$discrimination, params_2pl$difficulty,
     xlab = "Discrimination", ylab = "Difficulty",
     main = "Item Parameter Space")
```

#### 3-Parameter Logistic (3PL) Model

Adds pseudo-guessing parameter:

```r
fit_3pl <- fit_irt(responses, model = "3PL")

# Guessing parameters (lower asymptote)
guessing <- fit_3pl$item_parameters$guessing
print(guessing)
```

#### Model Comparison

```r
# Compare models
cat("Rasch AIC:", fit_rasch$AIC, "\n")
cat("2PL AIC:", fit_2pl$AIC, "\n")
cat("3PL AIC:", fit_3pl$AIC, "\n")

# Lower AIC is better
```

#### Practical Example: Educational Test

```r
# Simulate a 20-item math test
set.seed(123)
n_students <- 500
n_items <- 20

# Items range from easy (-2) to hard (+2)
true_difficulty <- seq(-2, 2, length.out = n_items)

# Student abilities ~ N(0, 1)
true_ability <- rnorm(n_students, 0, 1)

# Generate responses
test_data <- matrix(NA, n_students, n_items)
for (i in 1:n_students) {
  for (j in 1:n_items) {
    prob <- plogis(true_ability[i] - true_difficulty[j])
    test_data[i, j] <- rbinom(1, 1, prob)
  }
}

# Fit Rasch model
fit_math <- fit_irt(test_data, model = "Rasch")

# Identify struggling students (bottom 10%)
struggling_cutoff <- quantile(fit_math$person_abilities, 0.10)
struggling_students <- which(fit_math$person_abilities < struggling_cutoff)

cat("Number of struggling students:", length(struggling_students), "\n")

# Identify problematic items (too hard/easy)
extreme_items <- which(abs(fit_math$item_parameters$difficulty) > 2)
cat("Extreme items:", extreme_items, "\n")
```

---

### 3. Latent Class Analysis

#### Basic 2-Class Model

```r
# Create binary indicator data
# Rows = individuals, Columns = items
# Each cell is 0 or 1

data_lca <- matrix(c(
  1, 1, 0, 1, 0,
  0, 0, 1, 0, 1,
  1, 1, 1, 1, 0,
  # ... more rows
), ncol = 5, byrow = TRUE)

colnames(data_lca) <- paste0("Symptom", 1:5)

# Fit 2-class model
fit_2class <- fit_lca(data_lca, nclass = 2)

print(fit_2class)

# Class membership probabilities
cat("Class 1:", fit_2class$class_probs[1], "\n")
cat("Class 2:", fit_2class$class_probs[2], "\n")

# Item response probabilities by class
print(fit_2class$item_probs)
```

#### Model Selection

```r
# Fit models with different numbers of classes
fit_1 <- fit_lca(data_lca, nclass = 1)
fit_2 <- fit_lca(data_lca, nclass = 2)
fit_3 <- fit_lca(data_lca, nclass = 3)
fit_4 <- fit_lca(data_lca, nclass = 4)

# Compare BIC (lower is better)
bics <- c(fit_1$BIC, fit_2$BIC, fit_3$BIC, fit_4$BIC)
plot(1:4, bics, type = "b",
     xlab = "Number of Classes",
     ylab = "BIC",
     main = "Model Selection")

# Select best model
best_k <- which.min(bics)
cat("Best number of classes:", best_k, "\n")
```

#### Examining Class Membership

```r
# Get most likely class for each person
modal_classes <- fit_2class$modal_class
table(modal_classes)

# Examine posterior probabilities
head(fit_2class$posterior)

# How certain are classifications?
max_prob <- apply(fit_2class$posterior, 1, max)
hist(max_prob, main = "Classification Certainty",
     xlab = "Maximum Posterior Probability")

# Uncertain cases (max prob < 0.7)
uncertain <- which(max_prob < 0.7)
cat("Number of uncertain classifications:", length(uncertain), "\n")
```

#### Interpreting Classes

```r
# Examine item profiles
profiles <- fit_2class$item_probs

# Visualize
barplot(t(profiles), beside = TRUE,
        legend.text = paste("Class", 1:2),
        main = "Item Response Profiles",
        xlab = "Item", ylab = "P(Endorse)")

# Example interpretation:
# Class 1: High on all symptoms -> "Severe" group
# Class 2: Low on all symptoms -> "Mild" group
```

#### Practical Example: Depression Subtypes

```r
set.seed(456)
n_patients <- 500
symptoms <- c("Sadness", "Fatigue", "Insomnia",
              "Appetite", "Guilt", "Concentration")

# Three subtypes:
# 1. Severe (high on all)
# 2. Moderate (medium)
# 3. Mild (low on all)

severe_probs <- rep(0.85, 6)
moderate_probs <- rep(0.50, 6)
mild_probs <- rep(0.15, 6)

# Simulate
true_subtype <- sample(1:3, n_patients, replace = TRUE,
                       prob = c(0.2, 0.5, 0.3))

depression_data <- matrix(NA, n_patients, 6)
for (i in 1:n_patients) {
  probs <- switch(true_subtype[i],
                  severe_probs, moderate_probs, mild_probs)
  depression_data[i, ] <- rbinom(6, 1, probs)
}

colnames(depression_data) <- symptoms

# Fit 3-class model
fit_depression <- fit_lca(depression_data, nclass = 3)

print(fit_depression)

# Examine profiles
print(round(fit_depression$item_probs, 2))

# Class sizes
table(fit_depression$modal_class)
```

---

## Advanced Usage

### Model Diagnostics

#### GLMM Residuals

```r
fit <- gllamm(y ~ x + (1 | group), data = mydata)

# Different types of residuals
resid_response <- residuals(fit, type = "response")
resid_pearson <- residuals(fit, type = "pearson")
resid_deviance <- residuals(fit, type = "deviance")

# Diagnostic plots
par(mfrow = c(2, 2))
plot(fitted(fit), resid_response, main = "Residuals vs Fitted")
abline(h = 0, lty = 2, col = "red")

qqnorm(resid_response)
qqline(resid_response, col = "red")
```

#### IRT Item Fit

```r
fit_irt <- fit_irt(responses, model = "2PL")

# Items with extreme difficulties
extreme <- which(abs(fit_irt$item_parameters$difficulty) > 2.5)

# Items with low discrimination
low_discrim <- which(fit_irt$item_parameters$discrimination < 0.5)

cat("Problematic items:", union(extreme, low_discrim), "\n")
```

#### LCA Entropy

```r
# Measure classification quality
compute_entropy <- function(posterior) {
  posterior[posterior == 0] <- 1e-10
  entropy_i <- -rowSums(posterior * log(posterior))
  max_entropy <- log(ncol(posterior))
  rel_entropy <- 1 - mean(entropy_i) / max_entropy
  return(rel_entropy)
}

fit_lca <- fit_lca(data, nclass = 3)
entropy <- compute_entropy(fit_lca$posterior)
cat("Entropy:", round(entropy, 3), "\n")

# Entropy > 0.8 indicates good separation
```

### Prediction

```r
# GLMM predictions
fit <- gllamm(y ~ x + (1 | group), data = train_data)

# With random effects (training data)
pred_train <- predict(fit)

# Without random effects (population-average)
pred_pop <- predict(fit, re.form = NA)

# For new data (coming soon)
# pred_new <- predict(fit, newdata = test_data)
```

### Simulation

```r
# Simulate from fitted model
fit <- gllamm(y ~ x + (1 | group), data = mydata)

# Single simulation
sim1 <- simulate(fit, nsim = 1, seed = 123)

# Multiple simulations for uncertainty
sims <- simulate(fit, nsim = 100, seed = 123)
```

---

## Choosing the Right Model

### Decision Tree

**Do you have clustered/nested data?**
- Yes → Use **GLMM** (`gllamm()`)
  - Continuous outcome → `family = gaussian()`
  - Binary outcome → `family = binomial()`
  - Count outcome → `family = poisson()`

**Do you have test/questionnaire data?**
- Yes → Use **IRT** (`fit_irt()`)
  - Equal discrimination → Rasch
  - Varying discrimination → 2PL
  - Multiple choice (guessing) → 3PL

**Do you want to identify subgroups?**
- Yes → Use **LCA** (`fit_lca()`)
  - Fit models with 2, 3, 4,... classes
  - Select based on BIC
  - Interpret class profiles

---

## Model Comparison

### Information Criteria

```r
# Lower is better
cat("AIC:", AIC(fit), "\n")
cat("BIC:", BIC(fit), "\n")

# Compare two models
fit1 <- gllamm(y ~ x + (1 | group), data = mydata)
fit2 <- gllamm(y ~ x + z + (1 | group), data = mydata)

# BIC difference > 10 = strong evidence
bic_diff <- BIC(fit2) - BIC(fit1)
cat("BIC difference:", bic_diff, "\n")
```

### Likelihood Ratio Test (for nested models)

```r
# Coming soon
# anova(fit1, fit2)
```

---

## Tips & Best Practices

### GLMMs

1. **Start simple**: Begin with random intercepts, add slopes if needed
2. **Check convergence**: `fit$convergence$converged` should be TRUE
3. **Scale predictors**: Standardize continuous predictors for stability
4. **Sample size**: Need at least 5-10 obs per group for reliable RE

### IRT

1. **Sample size**:
   - Rasch: 200+ persons, 10+ items
   - 2PL: 500+ persons, 15+ items
   - 3PL: 1000+ persons, 20+ items
2. **Item selection**: Remove items with extreme difficulties
3. **Dimensionality**: Check assumption of unidimensionality
4. **Person fit**: Examine unusual response patterns

### LCA

1. **Model selection**: Use BIC, interpretability, and theory
2. **Random starts**: Use multiple starts (default: 3)
3. **Local independence**: Check assumption within classes
4. **Sample size**: 200+ observations, 5+ indicators
5. **Class size**: Avoid very small classes (< 5% of sample)

---

## Common Errors & Solutions

### Error: "TMB template not compiled"

**Solution**:
```r
TMB::compile("src/gllamm_gaussian.cpp")
# Or the appropriate template for your model
```

### Error: "Optimization failed"

**Solutions**:
- Check for missing data
- Try different starting values
- Increase iterations: `control = list(iter.max = 1000)`
- Simplify model (remove random slopes)

### Warning: "Failed to compute standard errors"

**Solutions**:
- Model may be too complex for data
- Try simpler random effects structure
- Check for perfect separation (binomial)
- Increase sample size

### LCA: "All optimization attempts failed"

**Solutions**:
- Reduce number of classes
- Increase random starts: `control = list(n_starts = 10)`
- Check data quality (too sparse?)
- Try different starting values

---

## Getting Help

### Documentation

```r
# Function help
?gllamm
?fit_irt
?fit_lca

# Package overview
help(package = "GLLAMMR")

# Vignettes
browseVignettes("GLLAMMR")
vignette("getting-started")
vignette("irt-models")
vignette("latent-class")
```

### Reporting Issues

- GitHub: https://github.com/yourusername/GLLAMMR/issues
- Include reproducible example
- Provide session info: `sessionInfo()`

---

## References

### GLLAMM Framework
- Rabe-Hesketh et al. (2004). GLLAMM Manual.
- Skrondal & Rabe-Hesketh (2004). Generalized Latent Variable Modeling.

### Item Response Theory
- Rasch (1960). Probabilistic Models.
- Lord (1980). Applications of IRT.
- Embretson & Reise (2000). IRT for Psychologists.

### Latent Class Analysis
- Lazarsfeld & Henry (1968). Latent Structure Analysis.
- Collins & Lanza (2010). Latent Class and Latent Transition Analysis.
- Hagenaars & McCutcheon (2002). Applied LCA.

---

## Citation

```r
citation("GLLAMMR")
```

If you use GLLAMMR in publications, please cite:

> GLLAMMR: Generalized Linear Latent and Mixed Models in R.
> Version 0.5.0 (2026). https://github.com/yourusername/GLLAMMR

---

**Version**: 0.5.0
**Last Updated**: 2026-02-06
**Maintainer**: Josh
**Contributors**: Claude Sonnet 4.5
