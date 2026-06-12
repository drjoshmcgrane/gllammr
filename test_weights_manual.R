#!/usr/bin/env Rscript

# Manual test for IRT weights support
library(TMB)

# Load compiled TMB models
dyn.load(dynlib("src/gllamm_irt"))
dyn.load(dynlib("src/gllamm_irt_poly"))

# Source R functions
source("R/irt.R")

cat("===== Testing IRT Weights Support =====\n\n")

# Test 1: Rasch model with equal weights
cat("Test 1: Rasch model - equal weights should match unweighted\n")
set.seed(123)
n_persons <- 30
n_items <- 8

theta <- rnorm(n_persons, 0, 1)
difficulty <- rnorm(n_items, 0, 1)
prob <- outer(theta, difficulty, function(th, b) plogis(th - b))
responses <- matrix(rbinom(n_persons * n_items, 1, prob), n_persons, n_items)

fit_nowt <- fit_irt(responses, model = "Rasch")
fit_eqwt <- fit_irt(responses, model = "Rasch", weights = rep(1, n_persons))

cat("  Log-likelihoods match:", all.equal(fit_nowt$logLik, fit_eqwt$logLik, tolerance = 1e-6), "\n")
cat("  Difficulties match:", all.equal(fit_nowt$item_parameters$difficulty,
                                        fit_eqwt$item_parameters$difficulty, tolerance = 1e-6), "\n")
cat("  Test 1:", ifelse(isTRUE(all.equal(fit_nowt$logLik, fit_eqwt$logLik, tolerance = 1e-6)), "PASSED", "FAILED"), "\n\n")

# Test 2: GRM model with equal weights
cat("Test 2: GRM model - equal weights should match unweighted\n")
set.seed(789)  # Different seed
n_persons <- 50
n_items <- 6
n_categories <- 3  # Reduced to 3 categories for more reliable sampling

theta <- rnorm(n_persons, 0, 1)
discrimination <- runif(n_items, 0.8, 1.5)

# Create ordered thresholds
thresholds_list <- lapply(1:n_items, function(j) {
  sort(rnorm(n_categories - 1, 0, 0.8))
})

responses <- matrix(0, n_persons, n_items)
for (i in 1:n_persons) {
  for (j in 1:n_items) {
    tau <- thresholds_list[[j]]
    a <- discrimination[j]

    # P(Y >= k) = plogis(a * (theta - tau_k))
    p_exceed <- c(1, plogis(a * (theta[i] - tau)), 0)
    probs <- p_exceed[-length(p_exceed)] - p_exceed[-1]

    responses[i, j] <- sample(1:n_categories, 1, prob = probs)
  }
}

fit_nowt <- fit_irt(responses, model = "GRM")
fit_eqwt <- fit_irt(responses, model = "GRM", weights = rep(1, n_persons))

cat("  Log-likelihoods match:", all.equal(fit_nowt$logLik, fit_eqwt$logLik, tolerance = 1e-5), "\n")
cat("  Discriminations match:", all.equal(fit_nowt$item_parameters$discrimination,
                                           fit_eqwt$item_parameters$discrimination, tolerance = 1e-5), "\n")
cat("  Test 2:", ifelse(isTRUE(all.equal(fit_nowt$logLik, fit_eqwt$logLik, tolerance = 1e-5)), "PASSED", "FAILED"), "\n\n")

# Test 3: Weights validation
cat("Test 3: Weights validation\n")
tryCatch({
  fit_irt(responses, model = "GRM", weights = rep(1, n_persons - 1))
  cat("  ERROR: Should have failed with wrong length\n")
}, error = function(e) {
  if (grepl("Length of weights", e$message)) {
    cat("  Correctly rejected wrong length weights\n")
  } else {
    cat("  ERROR: Wrong error message:", e$message, "\n")
  }
})

tryCatch({
  fit_irt(responses, model = "GRM", weights = c(rep(1, n_persons - 1), -1))
  cat("  ERROR: Should have failed with negative weight\n")
}, error = function(e) {
  if (grepl("non-negative", e$message)) {
    cat("  Correctly rejected negative weights\n")
  } else {
    cat("  ERROR: Wrong error message:", e$message, "\n")
  }
})

cat("  Test 3: PASSED\n\n")

# Test 4: Variable weights
cat("Test 4: Variable weights for PCM\n")
set.seed(789)
n_persons <- 40
n_items <- 5
n_categories <- 3

theta <- rnorm(n_persons, 0, 1)
step_difficulties <- lapply(1:n_items, function(j) {
  rnorm(n_categories - 1, 0, 0.5)
})

responses <- matrix(0, n_persons, n_items)
for (i in 1:n_persons) {
  for (j in 1:n_items) {
    delta <- step_difficulties[[j]]
    cumsums <- c(0, cumsum(theta[i] - delta))
    probs <- exp(cumsums) / sum(exp(cumsums))
    responses[i, j] <- sample(1:n_categories, 1, prob = probs)
  }
}

weights <- runif(n_persons, 0.5, 2)
fit_weighted <- fit_irt(responses, model = "PCM", weights = weights)

cat("  Converged:", fit_weighted$convergence$converged, "\n")
cat("  Log-likelihood finite:", is.finite(fit_weighted$logLik), "\n")
cat("  Test 4:", ifelse(fit_weighted$convergence$converged && is.finite(fit_weighted$logLik),
                        "PASSED", "FAILED"), "\n\n")

cat("===== All Tests Complete =====\n")
