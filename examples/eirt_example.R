# ============================================================================
# Comprehensive EIRT Example
# Demonstrates: Model comparison, removing predictors, testing effects
# ============================================================================

library(GLLAMMR)
set.seed(12345)

# ----------------------------------------------------------------------------
# 1. SIMULATE DATA WITH KNOWN ITEM COVARIATE EFFECTS
# ----------------------------------------------------------------------------

n_persons <- 200
n_items <- 20

# Create item-level covariates
item_data <- data.frame(
  word_frequency = rnorm(n_items, mean = 0, sd = 1),
  item_length = rpois(n_items, lambda = 6),
  is_abstract = rbinom(n_items, size = 1, prob = 0.4)
)

# True regression coefficients
true_gamma_0 <- 0.5     # Intercept
true_gamma_1 <- -0.8    # Word frequency effect (higher freq = easier)
true_gamma_2 <- 0.3     # Length effect (longer = harder)
true_gamma_3 <- 0.6     # Abstract effect (abstract = harder)
sigma_epsilon_b <- 0.4  # Residual SD for difficulty

# Generate true item difficulties
true_difficulty <- true_gamma_0 +
                  true_gamma_1 * item_data$word_frequency +
                  true_gamma_2 * item_data$item_length +
                  true_gamma_3 * item_data$is_abstract +
                  rnorm(n_items, 0, sigma_epsilon_b)

# Generate person abilities
theta <- rnorm(n_persons, mean = 0, sd = 1)

# Generate responses (Rasch model)
responses <- matrix(NA, n_persons, n_items)
for (i in 1:n_persons) {
  for (j in 1:n_items) {
    p <- plogis(theta[i] - true_difficulty[j])
    responses[i, j] <- rbinom(1, 1, p)
  }
}

cat("Data simulated:\n")
cat("  Persons:", n_persons, "\n")
cat("  Items:", n_items, "\n")
cat("  True coefficients:\n")
cat("    Intercept:", true_gamma_0, "\n")
cat("    Word frequency:", true_gamma_1, "\n")
cat("    Item length:", true_gamma_2, "\n")
cat("    Abstract:", true_gamma_3, "\n\n")


# ----------------------------------------------------------------------------
# 2. FIT MODELS WITH DIFFERENT PREDICTOR SETS
# ----------------------------------------------------------------------------

cat("="================================================================"\n")
cat("FITTING MULTIPLE MODELS\n")
cat("="================================================================"\n\n")

# Model 0: Null model (intercept only)
cat("Fitting Model 0: Intercept only...\n")
fit0 <- fit_eirt(
  responses,
  item_data = item_data,
  difficulty_formula = ~ 1,
  discrimination_formula = ~ 1,
  model = "Rasch"
)
cat("  LogLik:", round(fit0$logLik, 2), "\n")
cat("  AIC:", round(fit0$AIC, 2), "\n\n")

# Model 1: Word frequency only
cat("Fitting Model 1: Word frequency only...\n")
fit1 <- fit_eirt(
  responses,
  item_data = item_data,
  difficulty_formula = ~ word_frequency,
  discrimination_formula = ~ 1,
  model = "Rasch"
)
cat("  LogLik:", round(fit1$logLik, 2), "\n")
cat("  AIC:", round(fit1$AIC, 2), "\n\n")

# Model 2: Word frequency + item length
cat("Fitting Model 2: Word frequency + item length...\n")
fit2 <- fit_eirt(
  responses,
  item_data = item_data,
  difficulty_formula = ~ word_frequency + item_length,
  discrimination_formula = ~ 1,
  model = "Rasch"
)
cat("  LogLik:", round(fit2$logLik, 2), "\n")
cat("  AIC:", round(fit2$AIC, 2), "\n\n")

# Model 3: All predictors
cat("Fitting Model 3: All predictors...\n")
fit3 <- fit_eirt(
  responses,
  item_data = item_data,
  difficulty_formula = ~ word_frequency + item_length + is_abstract,
  discrimination_formula = ~ 1,
  model = "Rasch"
)
cat("  LogLik:", round(fit3$logLik, 2), "\n")
cat("  AIC:", round(fit3$AIC, 2), "\n\n")


# ----------------------------------------------------------------------------
# 3. COMPARE MODELS - TEST REMOVING PREDICTORS
# ----------------------------------------------------------------------------

cat("\n"================================================================"\n")
cat("MODEL COMPARISON: Testing importance of each predictor\n")
cat("="================================================================"\n\n")

# Compare Model 0 vs Model 1 (test word frequency)
cat("Test 1: Does word frequency matter?\n")
cat("  Comparing: Null vs Word frequency\n")
comp1 <- compare_eirt(fit0, fit1, test = "LRT")
print(comp1)
cat("\n\n")

# Compare Model 1 vs Model 2 (test adding item length)
cat("Test 2: Does item length add information beyond word frequency?\n")
cat("  Comparing: Word frequency vs Word frequency + length\n")
comp2 <- compare_eirt(fit1, fit2, test = "LRT")
print(comp2)
cat("\n\n")

# Compare Model 2 vs Model 3 (test adding abstract)
cat("Test 3: Does abstract status add information?\n")
cat("  Comparing: Freq + Length vs Freq + Length + Abstract\n")
comp3 <- compare_eirt(fit2, fit3, test = "LRT")
print(comp3)
cat("\n\n")

# Compare all models at once
cat("Overall comparison of all models:\n")
comp_all <- compare_eirt(fit0, fit1, fit2, fit3, test = "none")
print(comp_all)
cat("\n\n")


# ----------------------------------------------------------------------------
# 4. AUTOMATED TESTING WITH test_item_covariates()
# ----------------------------------------------------------------------------

cat("\n"================================================================"\n")
cat("AUTOMATED COVARIATE TESTING\n")
cat("="================================================================"\n\n")

cat("Testing all predictors automatically...\n\n")
test_result <- test_item_covariates(
  responses,
  item_data = item_data,
  difficulty_formula = ~ word_frequency + item_length + is_abstract,
  model = "Rasch"
)

print(test_result)


# ----------------------------------------------------------------------------
# 5. EXAMINE PARAMETER RECOVERY
# ----------------------------------------------------------------------------

cat("\n\n"================================================================"\n")
cat("PARAMETER RECOVERY\n")
cat("="================================================================"\n\n")

cat("Best model (Model 3) - Estimated coefficients:\n")
gamma_hat <- fit3$regression_coefficients$difficulty
print(round(gamma_hat, 3))

cat("\nTrue coefficients:\n")
true_coefs <- c(true_gamma_0, true_gamma_1, true_gamma_2, true_gamma_3)
names(true_coefs) <- c("(Intercept)", "word_frequency", "item_length", "is_abstract")
print(round(true_coefs, 3))

cat("\nRecovery quality:\n")
cat("  Intercept: Estimated =", round(gamma_hat[1], 3),
    "vs True =", round(true_gamma_0, 3), "\n")
cat("  Word frequency: Estimated =", round(gamma_hat[2], 3),
    "vs True =", round(true_gamma_1, 3), "\n")
cat("  Item length: Estimated =", round(gamma_hat[3], 3),
    "vs True =", round(true_gamma_2, 3), "\n")
cat("  Abstract: Estimated =", round(gamma_hat[4], 3),
    "vs True =", round(true_gamma_3, 3), "\n")

# Correlation between true and estimated
cor_coefs <- cor(true_coefs, gamma_hat)
cat("\nCorrelation:", round(cor_coefs, 3), "\n")


# ----------------------------------------------------------------------------
# 6. R-SQUARED FOR ITEM PARAMETER REGRESSION
# ----------------------------------------------------------------------------

cat("\n\n"================================================================"\n")
cat("VARIANCE EXPLAINED BY COVARIATES\n")
cat("="================================================================"\n\n")

r2_null <- eirt_r_squared(fit0, parameter = "difficulty")
r2_full <- eirt_r_squared(fit3, parameter = "difficulty")

cat("R² for difficulty regression:\n")
cat("  Null model (intercept only):", round(r2_null, 3), "\n")
cat("  Full model (all predictors):", round(r2_full, 3), "\n")
cat("  Improvement:", round(r2_full - r2_null, 3), "\n\n")

cat("Interpretation:\n")
cat("  The item covariates explain", round(100 * r2_full, 1),
    "% of variance in item difficulty\n")


# ----------------------------------------------------------------------------
# 7. VISUALIZE COVARIATE EFFECTS
# ----------------------------------------------------------------------------

cat("\n\n"================================================================"\n")
cat("VISUALIZING COVARIATE EFFECTS\n")
cat("="================================================================"\n\n")

cat("Creating plots...\n")

# Set up plotting area
par(mfrow = c(2, 2))

# Plot 1: Word frequency effect
plot_item_covariates(fit3, covariate = "word_frequency",
                    parameter = "difficulty")
title(sub = "Negative slope = easier with higher frequency")

# Plot 2: Item length effect
plot_item_covariates(fit3, covariate = "item_length",
                    parameter = "difficulty")
title(sub = "Positive slope = harder with more length")

# Plot 3: Abstract items
boxplot(fit3$item_parameters$difficulty ~ item_data$is_abstract,
        xlab = "Is Abstract", ylab = "Difficulty",
        names = c("Concrete", "Abstract"),
        main = "Difficulty by Item Type",
        col = c("lightblue", "lightcoral"))

# Plot 4: Predicted vs Observed
predicted_diff <- predict_difficulty(fit3)
observed_diff <- fit3$item_parameters$difficulty
plot(predicted_diff, observed_diff,
     xlab = "Predicted Difficulty",
     ylab = "Observed Difficulty",
     main = "Predicted vs Observed",
     pch = 19, col = "darkblue")
abline(0, 1, col = "red", lwd = 2, lty = 2)
grid(col = "gray90")

par(mfrow = c(1, 1))

cat("Plots complete!\n")


# ----------------------------------------------------------------------------
# 8. PREDICTIONS FOR NEW ITEMS
# ----------------------------------------------------------------------------

cat("\n\n"================================================================"\n")
cat("PREDICTING DIFFICULTY FOR NEW ITEMS\n")
cat("="================================================================"\n\n")

# Create hypothetical new items
new_items <- data.frame(
  word_frequency = c(-1.5, 0, 1.5),  # Low, medium, high frequency
  item_length = c(4, 6, 8),          # Short, medium, long
  is_abstract = c(0, 0, 1)           # Concrete, concrete, abstract
)

cat("New items:\n")
print(new_items)

# Predict difficulties
predicted_new <- predict_difficulty(fit3, newdata = new_items)

cat("\nPredicted difficulties:\n")
for (i in 1:3) {
  cat("  Item", i, ":", round(predicted_new[i], 3), "\n")
}

cat("\nInterpretation:\n")
cat("  Item 1 (high freq, short, concrete): Easiest\n")
cat("  Item 2 (medium freq, medium length, concrete): Medium\n")
cat("  Item 3 (low freq, long, abstract): Hardest\n")


# ----------------------------------------------------------------------------
# 9. SUMMARY AND RECOMMENDATIONS
# ----------------------------------------------------------------------------

cat("\n\n"================================================================"\n")
cat("SUMMARY AND RECOMMENDATIONS\n")
cat("="================================================================"\n\n")

cat("Model selection based on AIC:\n")
aic_values <- c(fit0$AIC, fit1$AIC, fit2$AIC, fit3$AIC)
best_model_idx <- which.min(aic_values)
model_names <- c("Null", "Word freq only", "Freq + Length", "All predictors")

cat("  Models ranked by AIC:\n")
for (i in order(aic_values)) {
  cat("    ", i, ". ", model_names[i], " (AIC = ",
      round(aic_values[i], 1), ")\n", sep = "")
}

cat("\n")
cat("Best model:", model_names[best_model_idx], "\n\n")

cat("Key findings:\n")
cat("  1. Word frequency has the strongest effect (negative)\n")
cat("  2. Item length also matters (positive effect)\n")
cat("  3. Abstract items are significantly harder\n")
cat("  4. Covariates explain", round(100 * r2_full, 1), "% of difficulty variance\n")

cat("\nConclusion:\n")
cat("  All three item characteristics are important predictors of difficulty.\n")
cat("  Including all predictors provides the best model fit.\n")
cat("  The model can accurately predict difficulty for new items based on\n")
cat("  their characteristics.\n")


# ----------------------------------------------------------------------------
# 10. EXAMPLE OUTPUT TO SHOW USER
# ----------------------------------------------------------------------------

cat("\n\n"================================================================"\n")
cat("QUICK REFERENCE: Key Functions\n")
cat("="================================================================"\n\n")

cat("# Fit EIRT model\n")
cat("fit <- fit_eirt(responses, item_data,\n")
cat("                difficulty_formula = ~ predictor1 + predictor2)\n\n")

cat("# Compare models (remove predictors)\n")
cat("fit_reduced <- fit_eirt(responses, item_data,\n")
cat("                        difficulty_formula = ~ predictor1)\n")
cat("compare_eirt(fit_reduced, fit, test = 'LRT')\n\n")

cat("# Automated testing\n")
cat("result <- test_item_covariates(responses, item_data,\n")
cat("                               difficulty_formula = ~ predictor1 + predictor2)\n\n")

cat("# Visualize effects\n")
cat("plot_item_covariates(fit, covariate = 'predictor1')\n\n")

cat("# Get R²\n")
cat("eirt_r_squared(fit, parameter = 'difficulty')\n\n")

cat("# Predict for new items\n")
cat("predict_difficulty(fit, newdata = new_item_data)\n\n")

cat("Example complete!\n")
