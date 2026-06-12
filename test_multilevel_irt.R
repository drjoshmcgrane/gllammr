# Test Multi-Level IRT Implementation

# Load package in development mode
devtools::load_all()

set.seed(123)

# Simulate nested data: 500 students in 20 classes
n_persons <- 500
n_items <- 15
n_classes <- 20

# Create person data with class membership
person_data <- data.frame(
  person_id = 1:n_persons,
  class_id = rep(1:n_classes, each = n_persons / n_classes)
)

# Simulate parameters
true_theta_0 <- rnorm(n_persons, 0, 1)      # Person deviations
true_u_class <- rnorm(n_classes, 0, 0.5)    # Class effects
true_difficulty <- rnorm(n_items, 0, 1)     # Item difficulties

# Compose total ability: theta = theta_0 + u_class[class]
true_theta <- true_theta_0 + true_u_class[person_data$class_id]

# Generate Rasch responses
responses <- matrix(NA, n_persons, n_items)
for (i in 1:n_persons) {
  for (j in 1:n_items) {
    p <- plogis(true_theta[i] - true_difficulty[j])
    responses[i, j] <- rbinom(1, 1, p)
  }
}

cat("=== Testing Multi-Level IRT ===\n\n")

# Test 1: Standard Rasch (no random effects) - should work as before
cat("Test 1: Standard Rasch model\n")
fit_standard <- fit_irt(responses, model = "Rasch")
cat("  Model class:", class(fit_standard), "\n")
cat("  Log-likelihood:", round(fit_standard$logLik, 2), "\n")
cat("  Person SD:", round(fit_standard$ability_sd, 3), "\n\n")

# Test 2: Multi-level Rasch with class random effects
cat("Test 2: Multi-level Rasch with class effects\n")
fit_multilevel <- fit_irt(
  responses,
  model = "Rasch",
  person_data = person_data,
  random = ~ (1 | class_id)
)
cat("  Model class:", class(fit_multilevel), "\n")
cat("  Log-likelihood:", round(fit_multilevel$logLik, 2), "\n")
cat("  Person SD:", round(fit_multilevel$ability_sd, 3), "\n")
cat("  Class SD:", round(fit_multilevel$random_effects$sigma_random, 3), "\n")
cat("  ICC (Class):", round(fit_multilevel$random_effects$icc["class_id"], 3), "\n")
cat("  ICC (Person):", round(fit_multilevel$random_effects$icc["Person"], 3), "\n\n")

# Test 3: Multi-level 2PL
cat("Test 3: Multi-level 2PL with class effects\n")
fit_2pl_multilevel <- fit_irt(
  responses,
  model = "2PL",
  person_data = person_data,
  random = ~ (1 | class_id)
)
cat("  Model class:", class(fit_2pl_multilevel), "\n")
cat("  Log-likelihood:", round(fit_2pl_multilevel$logLik, 2), "\n")
cat("  Class SD:", round(fit_2pl_multilevel$random_effects$sigma_random, 3), "\n\n")

# Test 4: Print method
cat("Test 4: Print method for multi-level model\n")
print(fit_multilevel)

# Test 5: Comparison - multi-level should fit better
cat("\nTest 5: Model comparison\n")
cat("  Standard Rasch AIC:", round(fit_standard$AIC, 2), "\n")
cat("  Multi-level Rasch AIC:", round(fit_multilevel$AIC, 2), "\n")
cat("  AIC improvement:", round(fit_standard$AIC - fit_multilevel$AIC, 2), "\n")

# Test 6: Parameter recovery
cat("\nTest 6: Parameter recovery\n")
cat("  True class SD: 0.5\n")
cat("  Estimated class SD:", round(fit_multilevel$random_effects$sigma_random, 3), "\n")
cat("  Correlation of estimated vs true total theta:",
    round(cor(fit_multilevel$random_effects$composite_theta, true_theta), 3), "\n")

cat("\n=== All tests completed successfully! ===\n")
