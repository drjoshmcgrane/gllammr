#!/usr/bin/env Rscript
# Standalone test script for EIRT API changes

# Load the package
library(GLLAMMR)

cat("\n=== Testing EIRT API Changes ===\n\n")

# Test 1: Check that LPCM is rejected
cat("Test 1: LPCM model name should be rejected...\n")
responses <- matrix(c(0, 1, 1, 0), 2, 2)
item_data <- data.frame(item_id = 1:2, x = c(0, 1))

tryCatch({
  fit_eirt(
    response_matrix = responses,
    item_data = item_data,
    model = "LPCM"
  )
  cat("  FAIL: LPCM was accepted (should have been rejected)\n")
}, error = function(e) {
  if (grepl("should be one of", e$message)) {
    cat("  PASS: LPCM correctly rejected\n")
  } else {
    cat("  FAIL: Wrong error message:", e$message, "\n")
  }
})

# Test 2: Check item_residuals parameter exists
cat("\nTest 2: item_residuals parameter should exist...\n")
eirt_formals <- formals(fit_eirt)
if ("item_residuals" %in% names(eirt_formals)) {
  if (identical(eirt_formals$item_residuals, TRUE)) {
    cat("  PASS: item_residuals parameter exists with default TRUE\n")
  } else {
    cat("  FAIL: item_residuals default is not TRUE\n")
  }
} else {
  cat("  FAIL: item_residuals parameter not found\n")
}

# Test 3: Test pure LLTM vs LLTM+error
cat("\nTest 3: Pure LLTM (item_residuals = FALSE) vs LLTM+error...\n")
set.seed(123)
n_persons <- 50
n_items <- 8

# True abilities
theta_true <- rnorm(n_persons, 0, 1)

# Item covariates
item_data <- data.frame(
  item_id = 1:n_items,
  difficulty_covar = rnorm(n_items, 0, 0.5)
)

# Generate difficulties from covariates only
gamma_true <- 0.8
b_true <- gamma_true * item_data$difficulty_covar

# Generate responses
responses <- matrix(NA, n_persons, n_items)
for (i in 1:n_persons) {
  for (j in 1:n_items) {
    prob <- plogis(theta_true[i] - b_true[j])
    responses[i, j] <- rbinom(1, 1, prob)
  }
}

tryCatch({
  # Fit pure LLTM
  fit_pure <- fit_eirt(
    response_matrix = responses,
    item_data = item_data,
    difficulty_formula = ~ difficulty_covar,
    model = "Rasch",
    item_residuals = FALSE
  )

  # Fit LLTM + error
  fit_error <- fit_eirt(
    response_matrix = responses,
    item_data = item_data,
    difficulty_formula = ~ difficulty_covar,
    model = "Rasch",
    item_residuals = TRUE
  )

  # Check classes
  if (inherits(fit_pure, "gllamm_eirt") && inherits(fit_error, "gllamm_eirt")) {
    cat("  PASS: Both models fit successfully\n")

    # Check item_residuals attribute
    if (!is.null(fit_pure$item_residuals) && !fit_pure$item_residuals) {
      cat("  PASS: Pure LLTM has item_residuals = FALSE\n")
    } else {
      cat("  INFO: Could not verify item_residuals attribute on pure LLTM\n")
    }

    if (!is.null(fit_error$item_residuals) && fit_error$item_residuals) {
      cat("  PASS: LLTM+error has item_residuals = TRUE\n")
    } else {
      cat("  INFO: Could not verify item_residuals attribute on LLTM+error\n")
    }

    # Compare number of parameters
    n_par_pure <- length(fit_pure$tmb_obj$par)
    n_par_error <- length(fit_error$tmb_obj$par)

    if (n_par_pure < n_par_error) {
      cat("  PASS: Pure LLTM has fewer parameters (", n_par_pure,
          ") than LLTM+error (", n_par_error, ")\n")
    } else {
      cat("  WARN: Pure LLTM does not have fewer parameters\n")
    }
  }
}, error = function(e) {
  cat("  FAIL: Error during fitting:", e$message, "\n")
})

# Test 4: Test 2PL with discrimination predictors
cat("\nTest 4: 2PL with discrimination predictors...\n")
set.seed(456)

# Item covariates
item_data <- data.frame(
  item_id = 1:n_items,
  difficulty_covar = rnorm(n_items, 0, 0.5),
  discrimination_covar = rnorm(n_items, 0, 0.3)
)

# Generate parameters
b_true <- 0.8 * item_data$difficulty_covar
a_true <- exp(0.5 * item_data$discrimination_covar)

# Generate responses
responses <- matrix(NA, n_persons, n_items)
for (i in 1:n_persons) {
  for (j in 1:n_items) {
    prob <- plogis(a_true[j] * (theta_true[i] - b_true[j]))
    responses[i, j] <- rbinom(1, 1, prob)
  }
}

tryCatch({
  fit_2pl <- fit_eirt(
    response_matrix = responses,
    item_data = item_data,
    difficulty_formula = ~ difficulty_covar,
    discrimination_formula = ~ discrimination_covar,
    model = "2PL"
  )

  if (inherits(fit_2pl, "gllamm_eirt")) {
    cat("  PASS: 2PL with discrimination predictors fitted successfully\n")

    if ("gamma" %in% names(fit_2pl$coefficients)) {
      cat("  PASS: Difficulty coefficients (gamma) present\n")
    }

    if ("delta" %in% names(fit_2pl$coefficients)) {
      cat("  PASS: Discrimination coefficients (delta) present\n")
    }
  }
}, error = function(e) {
  cat("  FAIL: Error during 2PL fitting:", e$message, "\n")
})

cat("\n=== EIRT API Tests Complete ===\n\n")
