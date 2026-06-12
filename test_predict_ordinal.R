#!/usr/bin/env Rscript
# Manual test script for ordinal marginal predictions

library(MASS)

# Try to load package, otherwise source files
if (!require(GLLAMMR, quietly = TRUE)) {
  source("R/formula.R")
  source("R/families.R")
  source("R/tmb_interface_v2.R")
  source("R/ordinal.R")
  source("R/predict_ordinal.R")
  source("R/marginal_utils.R")

  # Load TMB templates if needed
  if (file.exists("src/gllamm_ordinal.so")) {
    dyn.load("src/gllamm_ordinal.so")
  }
}

set.seed(123)

cat("=== Testing Ordinal Marginal Predictions ===\n\n")

# Test 1: Proportional Odds Model with Random Intercept
cat("Test 1: Proportional Odds (Logit) with Random Intercept\n")
cat("--------------------------------------------------------\n")

n_groups <- 20
n_per_group <- 15
n_obs <- n_groups * n_per_group

test_data <- data.frame(
  y = sample(1:4, n_obs, replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2)),
  x1 = rnorm(n_obs),
  x2 = rnorm(n_obs),
  group = rep(1:n_groups, each = n_per_group)
)

tryCatch({
  fit_ord <- fit_ordinal(
    y ~ x1 + x2 + (1 | group),
    data = test_data,
    link = "logit"
  )

  cat("✓ Model fitted successfully\n")
  cat("  Thresholds:", round(fit_ord$coefficients$thresholds, 3), "\n")
  cat("  Fixed effects:", round(fit_ord$coefficients$fixed, 3), "\n")
  cat("  Random SD:", round(sqrt(fit_ord$coefficients$random_var[[1]]), 3), "\n\n")

  # Test conditional predictions
  cat("Testing conditional predictions (type='probs')...\n")
  pred_cond <- predict(fit_ord, type = "probs")
  cat("✓ Conditional predictions computed\n")
  cat("  Dimensions:", paste(dim(pred_cond), collapse=" x "), "\n")
  cat("  First obs probs:", round(pred_cond[1,], 3), "\n")
  cat("  Sum to 1?", all(abs(rowSums(pred_cond) - 1) < 1e-6), "\n\n")

  # Test marginal predictions
  cat("Testing marginal predictions (type='marginal')...\n")
  pred_marg <- predict(fit_ord, type = "marginal", n_sim = 500)
  cat("✓ Marginal predictions computed\n")
  cat("  Dimensions:", paste(dim(pred_marg), collapse=" x "), "\n")
  cat("  First obs probs:", round(pred_marg[1,], 3), "\n")
  cat("  Sum to 1?", all(abs(rowSums(pred_marg) - 1) < 1e-6), "\n\n")

  # Compare conditional vs marginal
  cat("Comparing conditional vs marginal predictions:\n")
  cat("  Mean absolute difference:", round(mean(abs(pred_cond - pred_marg)), 4), "\n")
  cat("  Expected: marginal should differ from conditional (Jensen's inequality)\n\n")

  # Test class predictions
  cat("Testing class predictions...\n")
  pred_class <- predict(fit_ord, type = "class")
  cat("✓ Class predictions computed\n")
  cat("  Class distribution:", table(pred_class), "\n\n")

  cat("✓ Test 1 PASSED\n\n")

}, error = function(e) {
  cat("✗ Test 1 FAILED:", conditionMessage(e), "\n\n")
})


# Test 2: Probit Link
cat("Test 2: Proportional Odds (Probit) with Random Intercept\n")
cat("----------------------------------------------------------\n")

tryCatch({
  fit_probit <- fit_ordinal(
    y ~ x1 + (1 | group),
    data = test_data,
    link = "probit"
  )

  cat("✓ Probit model fitted successfully\n")

  # Marginal predictions
  pred_marg_probit <- predict(fit_probit, type = "marginal", n_sim = 500)
  cat("✓ Marginal predictions computed\n")
  cat("  First obs probs:", round(pred_marg_probit[1,], 3), "\n")
  cat("  All sum to 1?", all(abs(rowSums(pred_marg_probit) - 1) < 1e-6), "\n\n")

  cat("✓ Test 2 PASSED\n\n")

}, error = function(e) {
  cat("✗ Test 2 FAILED:", conditionMessage(e), "\n\n")
})


# Test 3: New Data Prediction
cat("Test 3: Predictions on New Data\n")
cat("---------------------------------\n")

new_data <- data.frame(
  x1 = c(-2, 0, 2),
  x2 = c(-1, 0, 1),
  group = c(1, 1, 1)  # Use existing group
)

tryCatch({
  # Conditional
  pred_new_cond <- predict(fit_ord, newdata = new_data, type = "probs")
  cat("✓ Conditional predictions on new data\n")
  cat("  Dimensions:", paste(dim(pred_new_cond), collapse=" x "), "\n")

  # Marginal
  pred_new_marg <- predict(fit_ord, newdata = new_data, type = "marginal", n_sim = 500)
  cat("✓ Marginal predictions on new data\n")
  cat("  Dimensions:", paste(dim(pred_new_marg), collapse=" x "), "\n")

  # Show predictions for different covariate values
  cat("\nPredictions across covariate range:\n")
  for (i in 1:nrow(new_data)) {
    cat(sprintf("  x1=%.1f, x2=%.1f: Marginal P(Y=4) = %.3f\n",
                new_data$x1[i], new_data$x2[i], pred_new_marg[i, 4]))
  }

  cat("\n✓ Test 3 PASSED\n\n")

}, error = function(e) {
  cat("✗ Test 3 FAILED:", conditionMessage(e), "\n\n")
})


# Test 4: Cumulative Probabilities
cat("Test 4: Cumulative Probability Predictions\n")
cat("--------------------------------------------\n")

tryCatch({
  pred_cumprobs <- predict(fit_ord, type = "cumprobs")
  cat("✓ Cumulative probabilities computed\n")
  cat("  Dimensions:", paste(dim(pred_cumprobs), collapse=" x "), "\n")
  cat("  First obs cum probs:", round(pred_cumprobs[1,], 3), "\n")
  cat("  Should be increasing:", all(apply(pred_cumprobs, 1, function(x) all(diff(x) >= -1e-6))), "\n\n")

  cat("✓ Test 4 PASSED\n\n")

}, error = function(e) {
  cat("✗ Test 4 FAILED:", conditionMessage(e), "\n\n")
})


cat("=== Ordinal Marginal Predictions Testing Complete ===\n")
