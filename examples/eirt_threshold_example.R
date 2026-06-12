# ============================================================================
# EIRT with Threshold-Level Predictors: Example
# Demonstrates threshold covariates in polytomous EIRT models
# ============================================================================

library(GLLAMMR)
set.seed(98765)

# ----------------------------------------------------------------------------
# 1. SIMULATE POLYTOMOUS DATA WITH THRESHOLD COVARIATES
# ----------------------------------------------------------------------------

n_persons <- 150
n_items <- 15
n_categories <- 4  # 4-point Likert scale (1, 2, 3, 4)

# Create item-level covariates
item_data <- data.frame(
  abstractness = rnorm(n_items, mean = 0, sd = 1),
  complexity = rnorm(n_items, mean = 0, sd = 1)
)

# True parameters
# Item-level difficulty (affects all thresholds)
true_gamma_0 <- 0.3
true_gamma_abstractness <- 0.5  # Abstract items harder

# Threshold-level effect
true_tau_abstractness <- -0.4   # Abstractness compresses thresholds

# Item discrimination
true_discrimination <- rep(1.2, n_items)

# Generate true difficulties (item-level)
true_difficulty <- true_gamma_0 + true_gamma_abstractness * item_data$abstractness

# Generate person abilities
theta <- rnorm(n_persons, mean = 0, sd = 1)

# Generate thresholds for each item
# Threshold spacing depends on abstractness
threshold_matrix <- matrix(NA, n_items, n_categories - 1)

for (j in 1:n_items) {
  # Base thresholds
  base_thresholds <- c(-1, 0, 1)

  # Modify spacing based on abstractness
  # More abstract items have compressed thresholds
  threshold_shift <- true_tau_abstractness * item_data$abstractness[j]

  threshold_matrix[j, 1] <- true_difficulty[j] + base_thresholds[1]
  threshold_matrix[j, 2] <- true_difficulty[j] + base_thresholds[2] + threshold_shift
  threshold_matrix[j, 3] <- true_difficulty[j] + base_thresholds[3] + 2 * threshold_shift
}

# Generate responses using GRM
responses <- matrix(NA, n_persons, n_items)

for (i in 1:n_persons) {
  for (j in 1:n_items) {
    # Cumulative probabilities
    cum_probs <- c(
      0,
      plogis(true_discrimination[j] * (theta[i] - threshold_matrix[j, 1])),
      plogis(true_discrimination[j] * (theta[i] - threshold_matrix[j, 2])),
      plogis(true_discrimination[j] * (theta[i] - threshold_matrix[j, 3])),
      1
    )

    # Category probabilities
    cat_probs <- diff(cum_probs)

    # Sample response
    responses[i, j] <- sample(1:n_categories, 1, prob = cat_probs)
  }
}

cat("Data simulated:\n")
cat("  Persons:", n_persons, "\n")
cat("  Items:", n_items, "\n")
cat("  Categories:", n_categories, "\n")
cat("  True parameters:\n")
cat("    Difficulty intercept:", true_gamma_0, "\n")
cat("    Difficulty ~ abstractness:", true_gamma_abstractness, "\n")
cat("    Threshold ~ abstractness:", true_tau_abstractness, "\n\n")


# ----------------------------------------------------------------------------
# 2. FIT EIRT MODELS: WITH AND WITHOUT THRESHOLD PREDICTORS
# ----------------------------------------------------------------------------

cat("=" , rep("=", 70), "\n", sep = "")
cat("FITTING MODELS\n")
cat("=" , rep("=", 70), "\n\n", sep = "")

# Model 1: Item-level predictors only (no threshold covariates)
cat("Model 1: Item-level difficulty predictors only...\n")
fit1 <- fit_eirt(
  response_matrix = responses,
  item_data = item_data,
  difficulty_formula = ~ abstractness,
  discrimination_formula = ~ 1,
  threshold_formula = NULL,  # No threshold predictors
  model = "GRM"
)
cat("  LogLik:", round(fit1$logLik, 2), "\n")
cat("  AIC:", round(fit1$AIC, 2), "\n\n")

# Model 2: Item-level + threshold-level predictors
cat("Model 2: Item-level AND threshold-level predictors...\n")
fit2 <- fit_eirt(
  response_matrix = responses,
  item_data = item_data,
  difficulty_formula = ~ abstractness,
  discrimination_formula = ~ 1,
  threshold_formula = ~ abstractness,  # Thresholds depend on abstractness
  model = "GRM"
)
cat("  LogLik:", round(fit2$logLik, 2), "\n")
cat("  AIC:", round(fit2$AIC, 2), "\n\n")


# ----------------------------------------------------------------------------
# 3. COMPARE MODELS
# ----------------------------------------------------------------------------

cat("\n", rep("=", 70), "\n", sep = "")
cat("MODEL COMPARISON\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Does abstractness affect threshold spacing?\n\n")

# Manual comparison
delta_loglik <- fit2$logLik - fit1$logLik
delta_df <- (length(fit2$tmb_obj$par) - length(fit1$tmb_obj$par))
lrt_stat <- 2 * delta_loglik
p_value <- pchisq(lrt_stat, df = delta_df, lower.tail = FALSE)

cat("Likelihood Ratio Test:\n")
cat("  LRT statistic:", round(lrt_stat, 2), "\n")
cat("  df:", delta_df, "\n")
cat("  p-value:", format.pval(p_value), "\n\n")

cat("Model fit comparison:\n")
cat("  Model 1 (no threshold predictors):\n")
cat("    LogLik:", round(fit1$logLik, 2), "\n")
cat("    AIC:", round(fit1$AIC, 2), "\n\n")

cat("  Model 2 (with threshold predictors):\n")
cat("    LogLik:", round(fit2$logLik, 2), "\n")
cat("    AIC:", round(fit2$AIC, 2), "\n\n")

if (p_value < 0.05) {
  cat("Conclusion: Threshold predictors significantly improve fit (p < 0.05)\n\n")
} else {
  cat("Conclusion: No significant improvement from threshold predictors\n\n")
}


# ----------------------------------------------------------------------------
# 4. EXAMINE PARAMETER RECOVERY
# ----------------------------------------------------------------------------

cat("\n", rep("=", 70), "\n", sep = "")
cat("PARAMETER RECOVERY (Model 2)\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Item-level difficulty regression:\n")
cat("  Intercept:\n")
cat("    True:", round(true_gamma_0, 3), "\n")
cat("    Estimated:", round(fit2$regression_coefficients$difficulty["(Intercept)"], 3), "\n\n")

cat("  Abstractness effect:\n")
cat("    True:", round(true_gamma_abstractness, 3), "\n")
cat("    Estimated:", round(fit2$regression_coefficients$difficulty["abstractness"], 3), "\n\n")

if (!is.null(fit2$regression_coefficients$threshold)) {
  cat("Threshold regression:\n")
  cat("  Abstractness effect on threshold spacing:\n")
  cat("    True:", round(true_tau_abstractness, 3), "\n")
  cat("    Estimated:", round(fit2$regression_coefficients$threshold["abstractness"], 3), "\n\n")
} else {
  cat("Threshold regression: NOT ESTIMATED\n\n")
}


# ----------------------------------------------------------------------------
# 5. PRINT FULL MODEL OUTPUT
# ----------------------------------------------------------------------------

cat("\n", rep("=", 70), "\n", sep = "")
cat("FULL MODEL OUTPUT\n")
cat(rep("=", 70), "\n\n", sep = "")

print(fit2)


# ----------------------------------------------------------------------------
# 6. SUMMARY AND INTERPRETATION
# ----------------------------------------------------------------------------

cat("\n\n", rep("=", 70), "\n", sep = "")
cat("SUMMARY AND INTERPRETATION\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Key Features Demonstrated:\n\n")

cat("1. ITEM-LEVEL PREDICTORS:\n")
cat("   - Difficulty regression: ~ abstractness\n")
cat("   - This shifts ALL thresholds for an item up or down\n")
cat("   - Abstract items are generally harder\n\n")

cat("2. THRESHOLD-LEVEL PREDICTORS:\n")
cat("   - Threshold regression: ~ abstractness\n")
cat("   - This affects SPACING between thresholds\n")
cat("   - Abstract items may have compressed or expanded threshold spacing\n\n")

cat("3. PRACTICAL APPLICATION:\n")
cat("   - Item-level: Overall item difficulty/easiness\n")
cat("   - Threshold-level: Response scale behavior (compressed vs expanded)\n")
cat("   - Example: Abstract items might make extreme categories harder to use\n\n")

cat("4. MODEL SELECTION:\n")
if (fit2$AIC < fit1$AIC) {
  cat("   - Model with threshold predictors has better AIC\n")
  cat("   - Threshold covariates improve model fit\n")
} else {
  cat("   - Model without threshold predictors has better AIC\n")
  cat("   - Threshold covariates may not be necessary for this data\n")
}

cat("\n\nExample complete!\n")
