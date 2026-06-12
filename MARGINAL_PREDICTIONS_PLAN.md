# Marginal Predictions Implementation Plan

## Overview

Marginal predictions (population-averaged predictions) integrate over the random effects distribution, as opposed to conditional predictions that condition on specific random effect values.

**Key Distinction:**
- **Conditional**: E[Y | X, u] - predictions given specific random effects u
- **Marginal**: E[Y | X] = ∫ E[Y | X, u] p(u) du - predictions averaged over random effects

This is critical for nonlinear models where marginal ≠ conditional predictions.

---

## Mathematical Background

### For GLMMs with link function g

**Conditional prediction:**
```
η_i = X_i'β + Z_i'u
μ_i|u = g^{-1}(η_i)
```

**Marginal prediction (no closed form for most links):**
```
μ_i = E[g^{-1}(X_i'β + Z_i'u)]
    = ∫ g^{-1}(X_i'β + Z_i'u) p(u) du
```

**Approximation Methods:**
1. **Gaussian quadrature** - accurate for 1D random effects
2. **Monte Carlo integration** - flexible for multi-dimensional u
3. **Laplace approximation** - fast but less accurate
4. **Adaptive Gaussian-Hermite quadrature** - best for low-dimensional u

We'll use **Monte Carlo integration** as it's flexible and works for all models.

---

## Implementation Strategy

### R Interface Design

Add `type` parameter to all predict() methods:

```r
predict.gllamm <- function(object, newdata = NULL,
                           type = c("link", "response", "marginal"),
                           n_sim = 1000, ...)
```

- **type = "link"**: Linear predictor (η = X'β + Z'u)
- **type = "response"**: Conditional mean (μ|u = g^{-1}(η))
- **type = "marginal"**: Population-averaged (μ = E[g^{-1}(η)])

### Monte Carlo Integration Algorithm

For each observation in newdata:

```r
# 1. Draw samples from random effects distribution
u_samples <- mvrnorm(n_sim, mu = 0, Sigma = Sigma_u)

# 2. Compute conditional prediction for each sample
mu_conditional <- sapply(u_samples, function(u) {
  eta <- X'beta + Z'u
  g_inverse(eta)  # Apply inverse link
})

# 3. Average across samples
mu_marginal <- mean(mu_conditional)
```

---

## Models Requiring Implementation

### Priority 1: GLMM Models

**Families:** Gaussian, Binomial, Poisson

**Why needed:** Nonlinear link functions (except Gaussian-identity)

**Implementation:**
- File: `R/predict.R` (new) or modify existing predict.gllamm()
- Method: Monte Carlo integration over u ~ N(0, Σ_u)
- For Gaussian with identity link: marginal = conditional (no integration needed)

**Example:**
```r
# Binomial GLMM with logit link
fit <- gllamm(y ~ x + (1 | group), data = data, family = binomial())

# Conditional prediction (at u = 0, i.e., "average" group)
pred_cond <- predict(fit, newdata, type = "response")  # E[Y|X, u=0]

# Marginal prediction (averaged over groups)
pred_marg <- predict(fit, newdata, type = "marginal")  # E[Y|X]
```

### Priority 2: Ordinal Models

**Models:** All ordinal link functions

**Why needed:** Nonlinear cumulative probabilities

**Implementation:**
- Extend predict.gllamm() to handle ordinal_family
- Return probability matrix (n_obs × n_categories) for marginal predictions
- Monte Carlo over u ~ N(0, Σ_u)

**Example:**
```r
fit <- gllamm(rating ~ temp + (1 | judge), family = ordinal(link = "logit"))

# Marginal probabilities for each category
pred_marg <- predict(fit, newdata, type = "marginal")
# Returns matrix: columns are P(Y=1), P(Y=2), ..., P(Y=K)
```

### Priority 3: IRT/EIRT Models

**Models:** All IRT and EIRT variants

**Why needed:** Population-level item response curves

**Implementation:**
- Add predict.gllamm_irt() and predict.gllamm_eirt()
- Monte Carlo over θ ~ N(0, σ²_θ)
- Return marginal item response probabilities

**Use cases:**
- Population-level difficulty estimates
- Expected test scores for random person
- Item characteristic curves averaged over ability distribution

**Example:**
```r
fit <- fit_irt(responses, model = "2PL")

# Conditional: P(correct | θ = 0)
pred_cond <- predict(fit, newdata, type = "response")

# Marginal: E[P(correct | θ)] over θ ~ N(0, σ²)
pred_marg <- predict(fit, newdata, type = "marginal")
```

### Priority 4: Multinomial Models

**Models:** Baseline category logit

**Why needed:** Nonlinear softmax transformation

**Implementation:**
- Monte Carlo over u ~ N(0, Σ_u)
- Return probability matrix for all categories

### Priority 5: Survival Models

**Models:** Exponential, Weibull with random effects

**Why needed:** Marginal hazard and survival curves

**Implementation:**
- Monte Carlo over u ~ N(0, Σ_u)
- Return marginal hazard, survival probability, or median survival time

---

## Implementation Details

### Step 1: Create Utility Function for MC Integration

File: **R/marginal_utils.R** (NEW)

```r
#' Monte Carlo integration for marginal predictions
#' @keywords internal
mc_integrate_marginal <- function(X, Z, beta, Sigma_u,
                                  inv_link_fn, n_sim = 1000) {
  n_obs <- nrow(X)
  n_random <- ncol(Z)

  # Draw samples from random effects distribution
  u_samples <- MASS::mvrnorm(n_sim, mu = rep(0, n_random), Sigma = Sigma_u)

  # Compute predictions for each sample
  predictions <- matrix(NA, n_obs, n_sim)
  for (s in 1:n_sim) {
    u <- u_samples[s, ]
    eta <- X %*% beta + Z %*% u
    predictions[, s] <- inv_link_fn(eta)
  }

  # Average across samples
  marginal_pred <- rowMeans(predictions)

  # Optional: return SE via simulation
  marginal_se <- apply(predictions, 1, sd) / sqrt(n_sim)

  list(
    fit = marginal_pred,
    se = marginal_se
  )
}
```

### Step 2: Extend predict.gllamm()

File: **R/predict.R** (NEW or modify existing)

```r
#' Predict method for gllamm objects
#' @export
predict.gllamm <- function(object, newdata = NULL,
                           type = c("link", "response", "marginal"),
                           n_sim = 1000, se.fit = FALSE, ...) {
  type <- match.arg(type)

  # Handle newdata
  if (is.null(newdata)) {
    X <- object$X
    Z <- object$Z  # May need to extract from object
  } else {
    # Parse newdata to get X and Z
    # ... (use existing formula infrastructure)
  }

  beta <- object$coefficients$fixed

  # Type: link
  if (type == "link") {
    eta <- X %*% beta
    # For conditional at u=0, just return eta
    return(as.vector(eta))
  }

  # Type: response (conditional at u=0)
  if (type == "response") {
    eta <- X %*% beta
    inv_link <- get_inverse_link(object$family)
    return(inv_link(eta))
  }

  # Type: marginal
  if (type == "marginal") {
    # Special case: Gaussian with identity link
    if (object$family$family == "gaussian" && object$family$link == "identity") {
      # Marginal = conditional for this case
      eta <- X %*% beta
      return(as.vector(eta))
    }

    # General case: Monte Carlo integration
    Sigma_u <- extract_random_vcov(object)
    inv_link <- get_inverse_link(object$family)

    result <- mc_integrate_marginal(
      X = X,
      Z = Z,  # Need to construct Z from newdata
      beta = beta,
      Sigma_u = Sigma_u,
      inv_link_fn = inv_link,
      n_sim = n_sim
    )

    if (se.fit) {
      return(list(fit = result$fit, se.fit = result$se))
    } else {
      return(result$fit)
    }
  }
}
```

### Step 3: Helper Functions

```r
#' Get inverse link function
#' @keywords internal
get_inverse_link <- function(family) {
  if (inherits(family, "binomial_family")) {
    switch(family$link,
      logit = plogis,
      probit = pnorm,
      cloglog = function(x) 1 - exp(-exp(x))
    )
  } else if (family$family == "poisson") {
    exp
  } else if (family$family == "gaussian") {
    identity
  } else {
    stop("Unknown family")
  }
}

#' Extract random effects variance-covariance matrix
#' @keywords internal
extract_random_vcov <- function(object) {
  # Extract Sigma_u from fitted object
  # This depends on how it's stored in the object
  # May need to reconstruct from log_sigma_u and theta parameters

  # For simple random intercept:
  if (length(object$coefficients$random_var) == 1) {
    sigma_u <- sqrt(object$coefficients$random_var[[1]])
    return(matrix(sigma_u^2, 1, 1))
  }

  # For multiple random effects:
  # Need to extract correlation matrix and SDs
  # ... (implementation depends on object structure)
}
```

### Step 4: Implement for Ordinal

```r
#' @export
predict.gllamm_ordinal <- function(object, newdata = NULL,
                                   type = c("class", "probs", "marginal"),
                                   n_sim = 1000, ...) {
  type <- match.arg(type)

  if (type == "marginal") {
    # Return marginal probability matrix (n_obs × n_categories)
    # Monte Carlo over u ~ N(0, Sigma_u)

    # For each observation and each MC sample:
    # 1. Draw u
    # 2. Compute conditional probabilities for all categories
    # 3. Average across samples

    # ... implementation ...
  }
}
```

### Step 5: Implement for IRT

```r
#' @export
predict.gllamm_irt <- function(object, newdata = NULL,
                               type = c("ability", "probability", "marginal"),
                               n_sim = 1000, ...) {
  type <- match.arg(type)

  if (type == "marginal") {
    # Marginal item response probability
    # E[P(correct | θ)] where θ ~ N(0, σ²_θ)

    # For each item:
    # 1. Draw θ samples
    # 2. Compute P(correct | θ, item_params)
    # 3. Average

    # ... implementation ...
  }
}
```

---

## Testing Strategy

### Unit Tests

**tests/testthat/test-marginal-predictions.R**

```r
test_that("Gaussian-identity: marginal equals conditional", {
  # For Gaussian with identity link, marginal = conditional
  fit <- gllamm(y ~ x + (1 | group), data = data, family = gaussian())

  pred_cond <- predict(fit, type = "response")
  pred_marg <- predict(fit, type = "marginal")

  expect_equal(pred_cond, pred_marg, tolerance = 1e-10)
})

test_that("Binomial-logit: marginal < conditional at u=0", {
  # Due to Jensen's inequality, E[g^{-1}(η)] < g^{-1}(E[η]) for convex g
  fit <- gllamm(y ~ x + (1 | group), data = data, family = binomial())

  pred_cond <- predict(fit, type = "response")  # At u=0
  pred_marg <- predict(fit, type = "marginal")

  # Marginal should be more conservative (closer to 0.5)
  expect_true(all(abs(pred_marg - 0.5) <= abs(pred_cond - 0.5)))
})

test_that("MC integration converges with more samples", {
  fit <- gllamm(y ~ x + (1 | group), data = data, family = binomial())

  pred_1k <- predict(fit, type = "marginal", n_sim = 1000)
  pred_10k <- predict(fit, type = "marginal", n_sim = 10000)

  # Should be very similar with more samples
  expect_equal(pred_1k, pred_10k, tolerance = 0.05)
})
```

### Integration Tests

Test with real data:
1. Fit binomial GLMM to clustered data
2. Compare marginal vs conditional predictions
3. Verify marginal predictions are more conservative (regression to mean)

---

## Performance Considerations

### Computational Cost

**n_sim = 1000** (default):
- For n = 100 observations: ~100,000 link function evaluations
- Fast for most link functions (plogis, exp)
- Bottleneck: Random number generation (MASS::mvrnorm)

**Optimization:**
- Vectorize across observations and samples
- Use efficient RNG (consider Rcpp if needed)
- Cache random samples if predicting on same data multiple times

### Accuracy vs Speed Trade-off

| n_sim | Speed | Accuracy | Recommendation |
|-------|-------|----------|----------------|
| 100 | Very fast | Low | Quick exploration |
| 1000 | Fast | Good | Default |
| 5000 | Moderate | Very good | Publication |
| 10000 | Slow | Excellent | High precision |

---

## Documentation Needs

### Function Documentation

```r
#' @param type Type of prediction:
#'   \describe{
#'     \item{"link"}{Linear predictor (η = X'β + Z'u|_{u=0})}
#'     \item{"response"}{Conditional mean (μ|u = g^{-1}(η)) at u=0}
#'     \item{"marginal"}{Population-averaged (μ = E[g^{-1}(η)]) integrated over u}
#'   }
#' @param n_sim Number of Monte Carlo samples for marginal predictions (default 1000)
#' @param se.fit Logical; return standard errors? (based on MC simulation variance)
```

### Vignette

Create `vignette("marginal-predictions")`:
- Explain conditional vs marginal
- When to use each type
- Interpretation of marginal predictions
- Examples with different families
- Computational considerations

---

## Implementation Phases

### Phase 1: Core Infrastructure (Priority 1)
**Duration:** 3-4 hours

1. Create `R/marginal_utils.R` with MC integration function
2. Create helper functions (get_inverse_link, extract_random_vcov)
3. Write unit tests for MC integration
4. Document utility functions

### Phase 2: GLMM predict() (Priority 1)
**Duration:** 3-4 hours

1. Create or modify `R/predict.R` with predict.gllamm()
2. Implement type = "link", "response", "marginal"
3. Handle newdata parsing
4. Write tests for all three types
5. Test with Gaussian, Binomial, Poisson

### Phase 3: Ordinal predict() (Priority 2)
**Duration:** 2-3 hours

1. Create predict.gllamm_ordinal() or extend predict.gllamm()
2. Handle ordinal_family dispatch
3. Return probability matrices for marginal
4. Write tests
5. Document

### Phase 4: IRT predict() (Priority 3)
**Duration:** 2-3 hours

1. Create predict.gllamm_irt()
2. Implement marginal item response probabilities
3. Write tests
4. Document

### Phase 5: Additional Models (Priority 4)
**Duration:** 3-4 hours

1. Multinomial predict()
2. Survival predict()
3. EIRT predict() (if different from IRT)

### Phase 6: Documentation & Polish (Priority 5)
**Duration:** 2-3 hours

1. Complete Rd documentation for all predict methods
2. Create marginal predictions vignette
3. Add examples to README
4. Final testing and validation

---

## Total Estimated Effort

**15-20 hours** for complete implementation across all models

**Next Step:** Implement Phase 1 (Core Infrastructure) and Phase 2 (GLMM predict())
