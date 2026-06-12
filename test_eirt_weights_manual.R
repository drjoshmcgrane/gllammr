#!/usr/bin/env Rscript

# Manual test for EIRT weights support
library(TMB)

# Load compiled TMB model
dyn.load(dynlib("src/gllamm_eirt"))

# Source R functions
source("R/eirt.R")

cat("===== Testing EIRT Weights Support =====\n\n")

# Test 1: Rasch EIRT with equal weights
cat("Test 1: Rasch EIRT - equal weights should match unweighted\n")
set.seed(123)
n_persons <- 40
n_items <- 10

# Create item data with predictor
item_data <- data.frame(
  item_id = 1:n_items,
  difficulty_pred = rnorm(n_items, 0, 0.5)
)

# Simulate abilities and item parameters
theta <- rnorm(n_persons, 0, 1)
gamma_true <- c(0, 0.7)  # Intercept and slope for difficulty formula

# Generate item difficulties from covariates
W_diff <- model.matrix(~ difficulty_pred, data = item_data)
difficulty <- as.vector(W_diff %*% gamma_true + rnorm(n_items, 0, 0.2))

# Simulate responses
responses <- matrix(0, n_persons, n_items)
for (i in 1:n_persons) {
  for (j in 1:n_items) {
    prob <- plogis(theta[i] - difficulty[j])
    responses[i, j] <- rbinom(1, 1, prob)
  }
}

# Fit without weights
fit_nowt <- fit_eirt(responses, item_data,
                     difficulty_formula = ~ difficulty_pred,
                     model = "Rasch")

# Fit with equal weights
fit_eqwt <- fit_eirt(responses, item_data,
                     difficulty_formula = ~ difficulty_pred,
                     weights = rep(1, n_persons),
                     model = "Rasch")

cat("  Log-likelihoods match:", all.equal(fit_nowt$logLik, fit_eqwt$logLik, tolerance = 1e-5), "\n")
cat("  Regression coefficients match:",
    all.equal(fit_nowt$regression_coefficients$difficulty,
              fit_eqwt$regression_coefficients$difficulty, tolerance = 1e-5), "\n")
cat("  Test 1:", ifelse(isTRUE(all.equal(fit_nowt$logLik, fit_eqwt$logLik, tolerance = 1e-5)),
                        "PASSED", "FAILED"), "\n\n")

# Test 2: GRM EIRT with variable weights
cat("Test 2: GRM EIRT - variable weights should work\n")
set.seed(456)
n_persons <- 50
n_items <- 8
n_categories <- 3

item_data <- data.frame(
  item_id = 1:n_items,
  abstractness = rnorm(n_items, 0, 1)
)

theta <- rnorm(n_persons, 0, 1)
gamma_true <- c(0, 0.5)
discrimination <- exp(rnorm(n_items, 0, 0.2))

W_diff <- model.matrix(~ abstractness, data = item_data)
difficulty <- as.vector(W_diff %*% gamma_true + rnorm(n_items, 0, 0.3))

# Simulate GRM responses
responses <- matrix(0, n_persons, n_items)
for (i in 1:n_persons) {
  for (j in 1:n_items) {
    # Create ordered thresholds around item difficulty
    tau <- difficulty[j] + c(-0.5, 0.5)
    a <- discrimination[j]

    p_exceed <- c(1, plogis(a * (theta[i] - tau)), 0)
    probs <- p_exceed[-length(p_exceed)] - p_exceed[-1]

    responses[i, j] <- sample(1:n_categories, 1, prob = probs)
  }
}

# Create variable weights
weights <- runif(n_persons, 0.5, 2)

# Fit with weights
fit_weighted <- fit_eirt(responses, item_data,
                        difficulty_formula = ~ abstractness,
                        weights = weights,
                        model = "GRM")

cat("  Converged:", fit_weighted$convergence$converged, "\n")
cat("  Log-likelihood finite:", is.finite(fit_weighted$logLik), "\n")
cat("  Regression coefficients exist:",
    !is.null(fit_weighted$regression_coefficients$difficulty), "\n")
cat("  Test 2:", ifelse(fit_weighted$convergence$converged && is.finite(fit_weighted$logLik),
                        "PASSED", "FAILED"), "\n\n")

# Test 3: LPCM with weights and threshold formula
cat("Test 3: LPCM with weights, difficulty formula, and threshold formula\n")
set.seed(789)
n_persons <- 50
n_items <- 8
n_categories <- 3

item_data <- data.frame(
  item_id = 1:n_items,
  difficulty_pred = rnorm(n_items, 0, 0.5),
  threshold_pred = rnorm(n_items, 0, 0.5)
)

theta <- rnorm(n_persons, 0, 1)

# Simulate LPCM data (simplified)
responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                   n_persons, n_items)

weights <- c(rep(1, 25), rep(2, 25))  # Some persons weighted more

fit_lpcm <- fit_eirt(responses, item_data,
                    difficulty_formula = ~ difficulty_pred,
                    threshold_formula = ~ threshold_pred,
                    weights = weights,
                    model = "LPCM")

cat("  Converged:", fit_lpcm$convergence$converged, "\n")
cat("  Has threshold coefficients:",
    !is.null(fit_lpcm$regression_coefficients$threshold), "\n")
cat("  Test 3:", ifelse(fit_lpcm$convergence$converged, "PASSED", "FAILED"), "\n\n")

cat("===== All EIRT Weights Tests Complete =====\n")
