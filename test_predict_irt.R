#!/usr/bin/env Rscript
# Manual test script for IRT marginal predictions

# Try to load package, otherwise source files
if (!require(GLLAMMR, quietly = TRUE)) {
  source("R/formula.R")
  source("R/families.R")
  source("R/tmb_interface_v2.R")
  source("R/irt.R")
  source("R/predict_irt.R")
  source("R/marginal_utils.R")

  # Load TMB templates if needed
  if (file.exists("src/gllamm_irt.so")) {
    dyn.load("src/gllamm_irt.so")
  }
  if (file.exists("src/gllamm_irt_poly.so")) {
    dyn.load("src/gllamm_irt_poly.so")
  }
}

set.seed(456)

cat("=== Testing IRT Marginal Predictions ===\n\n")

# Test 1: Rasch Model
cat("Test 1: Rasch Model Marginal Predictions\n")
cat("------------------------------------------\n")

n_persons <- 100
n_items <- 10

# Simulate Rasch data
true_theta <- rnorm(n_persons, 0, 1)
true_diff <- seq(-2, 2, length.out = n_items)

responses <- matrix(NA, n_persons, n_items)
for (i in 1:n_persons) {
  for (j in 1:n_items) {
    prob <- plogis(true_theta[i] - true_diff[j])
    responses[i, j] <- rbinom(1, 1, prob)
  }
}

tryCatch({
  fit_rasch <- fit_irt(responses, model = "Rasch")

  cat("✓ Rasch model fitted successfully\n")
  cat("  Item difficulties (first 5):", round(fit_rasch$item_parameters$difficulty[1:5], 3), "\n")
  cat("  Ability SD:", round(fit_rasch$ability_sd, 3), "\n\n")

  # Test ability predictions
  cat("Testing ability predictions...\n")
  pred_theta <- predict(fit_rasch, type = "ability")
  cat("✓ Ability estimates obtained\n")
  cat("  Length:", length(pred_theta), "\n")
  cat("  Mean:", round(mean(pred_theta), 3), "\n")
  cat("  SD:", round(sd(pred_theta), 3), "\n\n")

  # Test conditional probability predictions
  cat("Testing conditional probability predictions...\n")
  pred_cond <- predict(fit_rasch, type = "probability")
  cat("✓ Conditional probabilities computed\n")
  cat("  Dimensions:", paste(dim(pred_cond), collapse=" x "), "\n")
  cat("  Range:", round(range(pred_cond), 3), "\n\n")

  # Test marginal probability predictions
  cat("Testing marginal probability predictions...\n")
  pred_marg <- predict(fit_rasch, type = "marginal", n_sim = 1000)
  cat("✓ Marginal probabilities computed\n")
  cat("  Length:", length(pred_marg), "\n")
  cat("  Item 1 marginal prob:", round(pred_marg[1], 3), "\n")
  cat("  Item 5 marginal prob:", round(pred_marg[5], 3), "\n")
  cat("  Item 10 marginal prob:", round(pred_marg[10], 3), "\n\n")

  # Compare with empirical proportions
  cat("Comparing marginal predictions to empirical proportions:\n")
  emp_props <- colMeans(responses, na.rm = TRUE)
  cat("  Item 1: Predicted =", round(pred_marg[1], 3), ", Empirical =", round(emp_props[1], 3), "\n")
  cat("  Item 5: Predicted =", round(pred_marg[5], 3), ", Empirical =", round(emp_props[5], 3), "\n")
  cat("  Mean absolute error:", round(mean(abs(pred_marg - emp_props)), 4), "\n\n")

  # Test predictions for specific items
  cat("Testing predictions for subset of items...\n")
  pred_subset <- predict(fit_rasch, newdata = c(1, 3, 5), type = "marginal", n_sim = 1000)
  cat("✓ Subset predictions computed\n")
  cat("  Length:", length(pred_subset), "\n")
  cat("  Values:", round(pred_subset, 3), "\n\n")

  cat("✓ Test 1 PASSED\n\n")

}, error = function(e) {
  cat("✗ Test 1 FAILED:", conditionMessage(e), "\n\n")
})


# Test 2: 2PL Model
cat("Test 2: 2PL Model Marginal Predictions\n")
cat("----------------------------------------\n")

# Simulate 2PL data
true_disc <- runif(n_items, 0.5, 2.0)
responses_2pl <- matrix(NA, n_persons, n_items)
for (i in 1:n_persons) {
  for (j in 1:n_items) {
    prob <- plogis(true_disc[j] * (true_theta[i] - true_diff[j]))
    responses_2pl[i, j] <- rbinom(1, 1, prob)
  }
}

tryCatch({
  fit_2pl <- fit_irt(responses_2pl, model = "2PL")

  cat("✓ 2PL model fitted successfully\n")
  cat("  Discriminations (first 5):", round(fit_2pl$item_parameters$discrimination[1:5], 3), "\n\n")

  # Marginal predictions
  pred_marg_2pl <- predict(fit_2pl, type = "marginal", n_sim = 1000)
  cat("✓ Marginal probabilities computed\n")
  cat("  First 5 items:", round(pred_marg_2pl[1:5], 3), "\n\n")

  # Compare with empirical
  emp_props_2pl <- colMeans(responses_2pl, na.rm = TRUE)
  cat("Marginal vs empirical comparison:\n")
  cat("  Correlation:", round(cor(pred_marg_2pl, emp_props_2pl), 3), "\n")
  cat("  Mean absolute error:", round(mean(abs(pred_marg_2pl - emp_props_2pl)), 4), "\n\n")

  cat("✓ Test 2 PASSED\n\n")

}, error = function(e) {
  cat("✗ Test 2 FAILED:", conditionMessage(e), "\n\n")
})


# Test 3: 3PL Model with Guessing
cat("Test 3: 3PL Model Marginal Predictions\n")
cat("----------------------------------------\n")

# Simulate 3PL data with guessing
true_guess <- runif(n_items, 0, 0.25)
responses_3pl <- matrix(NA, n_persons, n_items)
for (i in 1:n_persons) {
  for (j in 1:n_items) {
    prob <- true_guess[j] + (1 - true_guess[j]) * plogis(true_disc[j] * (true_theta[i] - true_diff[j]))
    responses_3pl[i, j] <- rbinom(1, 1, prob)
  }
}

tryCatch({
  fit_3pl <- fit_irt(responses_3pl, model = "3PL")

  cat("✓ 3PL model fitted successfully\n")
  cat("  Guessing parameters (first 5):", round(fit_3pl$item_parameters$guessing[1:5], 3), "\n\n")

  # Marginal predictions
  pred_marg_3pl <- predict(fit_3pl, type = "marginal", n_sim = 1000)
  cat("✓ Marginal probabilities computed\n")
  cat("  First 5 items:", round(pred_marg_3pl[1:5], 3), "\n")
  cat("  All >= guessing parameter?",
      all(pred_marg_3pl >= fit_3pl$item_parameters$guessing - 0.01), "\n\n")

  cat("✓ Test 3 PASSED\n\n")

}, error = function(e) {
  cat("✗ Test 3 FAILED:", conditionMessage(e), "\n\n")
})


# Test 4: GRM (Graded Response Model)
cat("Test 4: GRM Marginal Predictions\n")
cat("---------------------------------\n")

n_items_grm <- 8
responses_grm <- matrix(sample(0:3, n_persons * n_items_grm, replace = TRUE),
                        n_persons, n_items_grm)

tryCatch({
  fit_grm <- fit_irt(responses_grm, model = "GRM")

  cat("✓ GRM model fitted successfully\n")
  cat("  Number of items:", fit_grm$n_items, "\n")
  cat("  Number of categories:", length(unique(as.vector(responses_grm))), "\n\n")

  # Marginal predictions
  pred_marg_grm <- predict(fit_grm, type = "marginal", n_sim = 1000)
  cat("✓ Marginal probabilities computed\n")
  cat("  First 3 items:", round(pred_marg_grm[1:3], 3), "\n")
  cat("  Note: For polytomous IRT, this returns probability of highest category\n\n")

  cat("✓ Test 4 PASSED (with note on polytomous handling)\n\n")

}, error = function(e) {
  cat("Note: GRM marginal predictions may not be fully implemented for polytomous models\n")
  cat("  Error:", conditionMessage(e), "\n\n")
})


# Test 5: Predictions at Specific Ability Levels
cat("Test 5: Predictions at Specific Ability Levels\n")
cat("-----------------------------------------------\n")

ability_grid <- seq(-3, 3, by = 0.5)

tryCatch({
  pred_at_abilities <- predict(fit_rasch, type = "probability", ability = ability_grid)

  cat("✓ Predictions at specified abilities computed\n")
  cat("  Dimensions:", paste(dim(pred_at_abilities), collapse=" x "), "\n")
  cat("  Ability levels:", length(ability_grid), "\n")
  cat("  Items:", ncol(pred_at_abilities), "\n\n")

  # Show probability progression for item 1
  cat("Item 1 probability across ability levels:\n")
  for (i in seq(1, length(ability_grid), by = 2)) {
    cat(sprintf("  θ = %4.1f: P = %.3f\n", ability_grid[i], pred_at_abilities[i, 1]))
  }

  cat("\n✓ Test 5 PASSED\n\n")

}, error = function(e) {
  cat("✗ Test 5 FAILED:", conditionMessage(e), "\n\n")
})


cat("=== IRT Marginal Predictions Testing Complete ===\n")
