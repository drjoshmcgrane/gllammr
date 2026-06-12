# Simple test of threshold-level EIRT
# Source files directly without package installation

# Load required packages
library(TMB)
library(Matrix)

# Source R files
source("R/eirt.R")

# Compile TMB template if needed
if (!file.exists("src/gllamm_eirt.so")) {
  TMB::compile("src/gllamm_eirt.cpp")
}
dyn.load(TMB::dynlib("src/gllamm_eirt"))

set.seed(123)

cat("=================================================================\n")
cat("Testing Threshold-Level Predictors in Polytomous EIRT\n")
cat("=================================================================\n\n")

# Simulate simple polytomous data
n_persons <- 100
n_items <- 10
n_categories <- 4

# Create item data with abstractness covariate
item_data <- data.frame(
  abstractness = rnorm(n_items, 0, 1)
)

cat("Item data created:\n")
print(head(item_data))
cat("\n")

# Simple simulation: abstractness affects both difficulty and threshold spacing
true_gamma <- 0.5  # Item-level effect
true_tau <- -0.3   # Threshold-level effect (negative = compressed thresholds)

# Generate person abilities
theta <- rnorm(n_persons, 0, 1)

# Generate responses
responses <- matrix(NA, n_persons, n_items)

for (i in 1:n_persons) {
  for (j in 1:n_items) {
    # Item difficulty based on abstractness
    difficulty <- true_gamma * item_data$abstractness[j]

    # Thresholds must be ordered: b1 < b2 < b3
    # Start with base thresholds
    base_thresholds <- c(-1, 0, 1)

    # Threshold spacing affected by abstractness
    # Positive abstractness compresses thresholds (tau is negative)
    tau_effect <- true_tau * item_data$abstractness[j]

    # Build ordered thresholds
    thresholds <- difficulty + base_thresholds
    # Modify spacing for thresholds 2 and 3
    thresholds[2] <- thresholds[2] + tau_effect * 0.3
    thresholds[3] <- thresholds[3] + tau_effect * 0.5

    # Ensure ordering (shouldn't be necessary but safe)
    if (thresholds[2] <= thresholds[1]) thresholds[2] <- thresholds[1] + 0.1
    if (thresholds[3] <= thresholds[2]) thresholds[3] <- thresholds[2] + 0.1

    # GRM probabilities
    a <- 1.2  # discrimination
    cum_probs <- c(0, plogis(a * (theta[i] - thresholds)), 1)
    cat_probs <- diff(cum_probs)

    # Check for valid probabilities
    if (any(cat_probs < 0) || abs(sum(cat_probs) - 1) > 1e-6) {
      # Fallback: uniform probabilities
      cat_probs <- rep(1/n_categories, n_categories)
    }

    # Sample response
    responses[i, j] <- sample(1:n_categories, 1, prob = cat_probs)
  }
}

cat("Responses generated: ", nrow(responses), "persons ×", ncol(responses), "items\n")
cat("Categories:", sort(unique(as.vector(responses))), "\n\n")

# Test 1: Fit model WITHOUT threshold predictors
cat("-------------------------------------------\n")
cat("Test 1: Model WITHOUT threshold predictors\n")
cat("-------------------------------------------\n")

fit1 <- tryCatch({
  fit_eirt(
    response_matrix = responses,
    item_data = item_data,
    difficulty_formula = ~ abstractness,
    discrimination_formula = ~ 1,
    threshold_formula = NULL,  # No threshold predictors
    model = "GRM"
  )
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  return(NULL)
})

if (!is.null(fit1)) {
  cat("✓ Model 1 fitted successfully\n")
  cat("  LogLik:", round(fit1$logLik, 2), "\n")
  cat("  AIC:", round(fit1$AIC, 2), "\n")
  cat("  Difficulty coefficients:\n")
  print(round(fit1$regression_coefficients$difficulty, 3))
  cat("  Threshold coefficients:",
      ifelse(is.null(fit1$regression_coefficients$threshold), "NULL", "NOT NULL"), "\n\n")
} else {
  cat("✗ Model 1 failed\n\n")
}

# Test 2: Fit model WITH threshold predictors
cat("-------------------------------------------\n")
cat("Test 2: Model WITH threshold predictors\n")
cat("-------------------------------------------\n")

fit2 <- tryCatch({
  fit_eirt(
    response_matrix = responses,
    item_data = item_data,
    difficulty_formula = ~ abstractness,
    discrimination_formula = ~ 1,
    threshold_formula = ~ abstractness,  # Threshold predictors!
    model = "GRM"
  )
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  return(NULL)
})

if (!is.null(fit2)) {
  cat("✓ Model 2 fitted successfully\n")
  cat("  LogLik:", round(fit2$logLik, 2), "\n")
  cat("  AIC:", round(fit2$AIC, 2), "\n")
  cat("  Difficulty coefficients:\n")
  print(round(fit2$regression_coefficients$difficulty, 3))
  cat("  Threshold coefficients:\n")
  if (!is.null(fit2$regression_coefficients$threshold)) {
    print(round(fit2$regression_coefficients$threshold, 3))
  } else {
    cat("    NULL (ERROR!)\n")
  }
  cat("\n")
} else {
  cat("✗ Model 2 failed\n\n")
}

# Compare models
cat("=================================================================\n")
cat("COMPARISON\n")
cat("=================================================================\n\n")

if (!is.null(fit1) && !is.null(fit2)) {
  cat("Model 1 (no threshold predictors):\n")
  cat("  LogLik:", round(fit1$logLik, 2), "\n")
  cat("  AIC:", round(fit1$AIC, 2), "\n\n")

  cat("Model 2 (with threshold predictors):\n")
  cat("  LogLik:", round(fit2$logLik, 2), "\n")
  cat("  AIC:", round(fit2$AIC, 2), "\n\n")

  delta_loglik <- fit2$logLik - fit1$logLik
  lrt_stat <- 2 * delta_loglik
  cat("LRT statistic:", round(lrt_stat, 2), "\n")

  if (fit2$AIC < fit1$AIC) {
    cat("✓ Model 2 (with threshold predictors) has better AIC\n")
  } else {
    cat("  Model 1 (without threshold predictors) has better AIC\n")
  }

  cat("\nParameter recovery:\n")
  cat("  True difficulty effect:", true_gamma, "\n")
  cat("  Estimated:", round(fit2$regression_coefficients$difficulty["abstractness"], 3), "\n")
  cat("  True threshold effect:", true_tau, "\n")
  if (!is.null(fit2$regression_coefficients$threshold)) {
    cat("  Estimated:", round(fit2$regression_coefficients$threshold["abstractness"], 3), "\n")
  }
}

cat("\n=================================================================\n")
cat("TEST COMPLETE\n")
cat("=================================================================\n")
