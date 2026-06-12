# Test if base polytomous EIRT works at all
# Based on existing test suite but without skips

library(TMB)
library(Matrix)

# Source and load
source("R/eirt.R")
dyn.load(TMB::dynlib("src/gllamm_eirt"))

cat("Testing basic polytomous EIRT (no threshold covariates)\n")
cat("=========================================================\n\n")

set.seed(999)
n_persons <- 60
n_items <- 8
n_categories <- 4

# Use simple uniform random assignment (most robust)
responses <- matrix(
  sample(1:n_categories, n_persons * n_items, replace = TRUE),
  n_persons, n_items
)

item_data <- data.frame(
  x = rnorm(n_items)
)

cat("Data:\n")
cat("  Persons:", n_persons, "\n")
cat("  Items:", n_items, "\n")
cat("  Categories:", n_categories, "\n")
cat("  Response range:", min(responses), "to", max(responses), "\n")
cat("  Item covariate range:", round(range(item_data$x), 2), "\n\n")

cat("Fitting model...\n")

fit <- tryCatch({
  fit_eirt(
    response_matrix = responses,
    item_data = item_data,
    difficulty_formula = ~ x,
    discrimination_formula = ~ 1,
    threshold_formula = NULL,  # No threshold covariates
    model = "GRM",
    control = list(eval.max = 2000, iter.max = 1000, trace = 1)
  )
}, error = function(e) {
  cat("\nERROR:", e$message, "\n")
  return(NULL)
})

if (!is.null(fit)) {
  cat("\n✓ Model fitted successfully!\n\n")
  cat("Results:\n")
  cat("  Converged:", fit$convergence$converged, "\n")
  cat("  LogLik:", round(fit$logLik, 2), "\n")
  cat("  AIC:", round(fit$AIC, 2), "\n")
  cat("  BIC:", round(fit$BIC, 2), "\n\n")

  cat("Difficulty regression coefficients:\n")
  print(fit$regression_coefficients$difficulty)
  cat("\n")

  cat("Discrimination regression coefficients:\n")
  print(fit$regression_coefficients$discrimination)
  cat("\n")

  cat("Residual SDs:\n")
  cat("  Difficulty:", round(fit$residual_sd$difficulty, 3), "\n")
  cat("  Discrimination:", round(fit$residual_sd$discrimination, 3), "\n")
  cat("  Threshold:", round(fit$residual_sd$threshold, 3), "\n")

} else {
  cat("\n✗ Model fitting failed\n")
}

cat("\n=========================================================\n")
