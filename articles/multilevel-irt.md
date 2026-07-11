# Multi-Level IRT Models in gllammr

## Introduction

Note: the code in this vignette is illustrative and is not evaluated
when the vignette is built.

This vignette demonstrates how to fit multi-level (hierarchical) Item
Response Theory (IRT) models using gllammr. Multi-level IRT models are
essential when test-takers are clustered within groups (e.g., students
within classrooms, patients within hospitals, or repeated measurements
within individuals).

### Why Multi-Level IRT?

Traditional IRT models assume independence across persons. However, when
persons are nested within groups, this assumption is often violated:

- **Students within classrooms**: Students in the same class share the
  same teacher and learning environment
- **Patients within hospitals**: Patients at the same hospital share
  institutional factors
- **Longitudinal assessments**: Multiple observations from the same
  person are correlated

Ignoring clustering can lead to:

- Underestimated standard errors
- Biased parameter estimates
- Invalid statistical inferences

### Model Specification

A multi-level 2PL model decomposes person ability into multiple
components:

``` math
\theta_i = \theta_{0i} + \sum_{g=1}^G u_{g[i]}
```

Where:

- $`\theta_i`$ is the total ability for person $`i`$
- $`\theta_{0i}`$ is the person-specific deviation (individual effect)
- $`u_{g[i]}`$ are random effects for groups $`g`$ (class effects,
  school effects, etc.)

The item response function becomes:

``` math
P(Y_{ij} = 1 | \theta_i) = \text{logit}^{-1}(a_j(\theta_i - b_j))
```

where the total ability $`\theta_i`$ includes both person and group
components.

## Basic Example: Students Nested in Classes

Let’s simulate data from 200 students nested in 20 classes:

``` r

library(gllammr)
set.seed(42)

# Setup
n_students <- 200
n_classes <- 20
n_items <- 15

# Create person-level data
person_data <- data.frame(
  student_id = 1:n_students,
  class_id = rep(1:n_classes, each = n_students / n_classes)
)

# Generate abilities with class effects
sigma_person <- 1.0   # Person-level SD
sigma_class <- 0.5    # Class-level SD

theta_0 <- rnorm(n_students, 0, sigma_person)  # Individual deviations
u_class <- rnorm(n_classes, 0, sigma_class)    # Class effects
theta_total <- theta_0 + u_class[person_data$class_id]

# Generate item parameters
difficulty <- rnorm(n_items, 0, 1)
discrimination <- rlnorm(n_items, 0, 0.3)

# Generate responses
responses <- matrix(NA, n_students, n_items)
for (i in 1:n_students) {
  for (j in 1:n_items) {
    p <- plogis(discrimination[j] * (theta_total[i] - difficulty[j]))
    responses[i, j] <- rbinom(1, 1, p)
  }
}

head(person_data)
```

### Fitting the Standard Model (Ignoring Clustering)

First, let’s fit a standard 2PL model that ignores the class structure:

``` r

fit_standard <- fit_irt(responses, model = "2PL")
print(fit_standard)
```

The standard model estimates a single variance component for person
abilities.

### Fitting the Multi-Level Model

Now, let’s properly account for the class structure using the
`person_data` and `random` arguments:

``` r

fit_multilevel <- fit_irt(
  response_matrix = responses,
  model = "2PL",
  person_data = person_data,
  random = ~ (1 | class_id)
)

print(fit_multilevel)
```

The output shows:

1.  **Variance Components**: Separate variances for class effects,
    person effects, and residual
2.  **Intraclass Correlations (ICCs)**: Proportion of variance at each
    level
3.  **Group Information**: Number of classes detected

### Extracting Variance Components

Use
[`VarCorr()`](https://drjoshmcgrane.github.io/gllammr/reference/VarCorr.md)
to get detailed variance decomposition:

``` r

vc <- VarCorr(fit_multilevel)
print(vc)
```

This shows:

- **class_id**: Variance due to differences between classes
  ($`\sigma^2_{class}`$)
- **Person**: Variance due to differences between students within
  classes ($`\sigma^2_{person}`$)
- **Residual**: Residual variance from logistic distribution
  ($`\pi^2/3 \approx 3.29`$)

### Computing Intraclass Correlations

The ICC tells us how much of the total variance is attributable to each
level:

``` r

# All ICCs
icc_all <- icc(fit_multilevel)
print(icc_all)

# Specific level
icc_class <- icc(fit_multilevel, level = "class_id")
cat("\nICC for class level:", round(icc_class, 3), "\n")
```

**Interpretation**: If ICC(class) = 0.10, then 10% of the total variance
in abilities is between classes, while 90% is within classes
(individuals + residual).

### Extracting Random Effects

Get the estimated class effects:

``` r

# All random effects
class_effects <- ranef(fit_multilevel, level = "class_id")
head(class_effects, 10)

# Plot class effects
plot(class_effects, main = "Estimated Class Effects",
     xlab = "Class", ylab = "Effect (u_class)")
abline(h = 0, lty = 2, col = "red")
```

Positive values indicate classes with higher-than-average ability,
negative values indicate lower-than-average ability.

### Extracting Person Abilities

For multi-level models, you can extract:

1.  **Person-level deviations** ($`\theta_{0i}`$): Individual effects
    within their class
2.  **Composite abilities** ($`\theta_i = \theta_{0i} + u_{class}`$):
    Total ability including class effect

``` r

# Person-level deviations only
theta_0 <- abilities(fit_multilevel, composite = FALSE)

# Composite abilities (person + class)
theta_composite <- abilities(fit_multilevel, composite = TRUE)

# Compare
plot(theta_0, theta_composite,
     xlab = "Person Deviation (theta_0)",
     ylab = "Composite Ability (theta_0 + u_class)",
     main = "Person vs Composite Abilities")
abline(0, 1, lty = 2, col = "red")
```

Students from the same class will have the same vertical offset from the
diagonal.

### Model Comparison

Compare log-likelihoods to see if the multi-level model fits better:

``` r

cat("Standard Model:\n")
cat("  Log-likelihood:", round(fit_standard$logLik, 2), "\n")
cat("  AIC:", round(fit_standard$AIC, 2), "\n\n")

cat("Multi-level Model:\n")
cat("  Log-likelihood:", round(fit_multilevel$logLik, 2), "\n")
cat("  AIC:", round(fit_multilevel$AIC, 2), "\n\n")

# Likelihood ratio test
lr_stat <- 2 * (fit_multilevel$logLik - fit_standard$logLik)
df_diff <- 1  # One additional variance parameter
p_value <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)

cat("Likelihood Ratio Test:\n")
cat("  LR statistic:", round(lr_stat, 2), "\n")
cat("  df:", df_diff, "\n")
cat("  p-value:", format.pval(p_value), "\n")
```

A significant p-value indicates that accounting for clustering
significantly improves model fit.

## Advanced: Multiple Levels

### Nested Structure: Schools \> Classes \> Students

When students are nested in classes, which are nested in schools:

``` r

# Create nested structure
person_data_nested <- data.frame(
  student_id = 1:200,
  school_id = rep(1:5, each = 40),
  class_id = rep(1:20, each = 10)
)

# Nested notation (school/class)
fit_nested1 <- fit_irt(
  responses,
  model = "2PL",
  person_data = person_data_nested,
  random = ~ (1 | school_id/class_id)
)

# Equivalent explicit notation
fit_nested2 <- fit_irt(
  responses,
  model = "2PL",
  person_data = person_data_nested,
  random = ~ (1 | school_id) + (1 | class_id)
)

# Both specifications are equivalent
print(fit_nested1)
```

This model decomposes ability into three components:

``` math
\theta_i = \theta_{0i} + u_{school[i]} + u_{class[i]}
```

### Crossed Effects: Students × Time

For longitudinal data where students are measured at multiple time
points:

``` r

# Simulate longitudinal data
n_students <- 100
n_timepoints <- 3
n_items <- 10

long_data <- expand.grid(
  student_id = 1:n_students,
  time = 1:n_timepoints
)
long_data$obs_id <- 1:nrow(long_data)

# Generate responses (rows = observations, cols = items)
responses_long <- matrix(rbinom(nrow(long_data) * n_items, 1, 0.6),
                         nrow(long_data), n_items)

# Crossed random effects
fit_crossed <- fit_irt(
  responses_long,
  model = "2PL",
  person_data = long_data,
  random = ~ (1 | student_id) + (1 | time)
)

print(fit_crossed)
```

This separates student effects (consistent individual differences) from
time effects (shared learning/fatigue across all students at each
occasion).

### Partial Nesting: Some Students Not in Classes

Some students might not belong to any class (e.g., independent study):

``` r

# Set some students to NA (not in any class)
person_data_partial <- person_data
person_data_partial$class_id[191:200] <- NA

# Fit with partial nesting
fit_partial <- fit_irt(
  responses,
  model = "2PL",
  person_data = person_data_partial,
  random = ~ (1 | class_id)
)

print(fit_partial)
# Note: Number of groups will be less than 20 (excludes NA class)
```

Students with `NA` class get no class effect (only $`\theta_{0i}`$),
while others get $`\theta_{0i} + u_{class}`$.

## Polytomous Multi-Level IRT

Multi-level models work with all IRT model types, including polytomous
models:

``` r

# Simulate 4-category responses
responses_poly <- matrix(sample(1:4, n_students * n_items, replace = TRUE),
                         n_students, n_items)

# Fit multi-level GRM
fit_grm_ml <- fit_irt(
  responses_poly,
  model = "GRM",
  person_data = person_data,
  random = ~ (1 | class_id)
)

print(fit_grm_ml)
```

All S3 methods
([`VarCorr()`](https://drjoshmcgrane.github.io/gllammr/reference/VarCorr.md),
[`icc()`](https://drjoshmcgrane.github.io/gllammr/reference/icc.md),
[`ranef()`](https://drjoshmcgrane.github.io/gllammr/reference/ranef.md),
[`abilities()`](https://drjoshmcgrane.github.io/gllammr/reference/abilities.md))
work the same way for polytomous models.

## Practical Recommendations

### When to Use Multi-Level IRT

Use multi-level IRT when:

- Test-takers are clearly clustered (schools, hospitals, sites)
- ICC \> 0.05 (more than 5% variance between groups)
- Ignoring clustering leads to artificially narrow confidence intervals

### Choosing Random Effects Structure

1.  **Start simple**: Begin with single-level random intercepts
2.  **Build up**: Add levels only if theoretically justified
3.  **Check ICCs**: If ICC \< 0.01, clustering may be negligible
4.  **Use LRT**: Compare nested models via likelihood ratio tests

### Interpretation Caveats

- **Composite abilities**: Use for individual predictions (accounts for
  group membership)
- **Person deviations**: Use for within-group comparisons (removes group
  effect)
- **Class effects**: Interpretable only if classes are random sample
  from population

### Convergence Tips

If the model fails to converge:

1.  Check for sparse groups (too few students per class)
2.  Simplify the model (fewer random effects levels)
3.  Use informative starting values
4.  Consider whether the data truly have sufficient information to
    estimate all variance components

## Summary

Multi-level IRT in gllammr:

- **Syntax**: Use `person_data` + `random = ~ (1 | group)` with any IRT
  model
- **Structures**: Supports nested, crossed, and partially nested designs
- **Methods**:
  [`VarCorr()`](https://drjoshmcgrane.github.io/gllammr/reference/VarCorr.md),
  [`icc()`](https://drjoshmcgrane.github.io/gllammr/reference/icc.md),
  [`ranef()`](https://drjoshmcgrane.github.io/gllammr/reference/ranef.md),
  [`abilities()`](https://drjoshmcgrane.github.io/gllammr/reference/abilities.md)
  for comprehensive output
- **Models**: Works with Rasch, 2PL, 3PL, GRM, PCM, GPCM, NRM

For more details, see:

- [`?fit_irt`](https://drjoshmcgrane.github.io/gllammr/reference/fit_irt.md)
  for full parameter documentation
- [`?VarCorr.gllamm_irt_multilevel`](https://drjoshmcgrane.github.io/gllammr/reference/VarCorr.gllamm_irt_multilevel.md)
  for variance component extraction
- [`?abilities`](https://drjoshmcgrane.github.io/gllammr/reference/abilities.md)
  for person ability extraction options

## Session Info

``` r

sessionInfo()
```
