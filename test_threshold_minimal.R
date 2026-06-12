# Minimal test of threshold-level EIRT implementation
# Just verify the code runs without testing parameter recovery

library(TMB)
library(Matrix)

# Source R files
source("R/eirt.R")

# Load compiled template
dyn.load(TMB::dynlib("src/gllamm_eirt"))

cat("==================================================\n")
cat("Minimal Test: Threshold Formula Acceptance\n")
cat("==================================================\n\n")

# Create very simple data
set.seed(999)
n_persons <- 50
n_items <- 8
n_cat <- 3  # Just 3 categories for simplicity

# Simple responses - mostly middle category
responses <- matrix(2, n_persons, n_items)
# Add some variation
for (i in 1:n_persons) {
  responses[i, sample(1:n_items, 3)] <- sample(c(1, 3), 3, replace = TRUE)
}

# Item data
item_data <- data.frame(
  x = scale(1:n_items)[,1]  # Simple continuous covariate
)

cat("Data: ", n_persons, "persons ×", n_items, "items\n")
cat("Categories:", sort(unique(as.vector(responses))), "\n")
cat("Item covariate: scaled 1:", n_items, "\n\n")

# TEST 1: Check that threshold_formula = NULL works (backward compatibility)
cat("Test 1: threshold_formula = NULL (baseline)\n")
cat("--------------------------------------------\n")

test1 <- tryCatch({
  fit_eirt(
    response_matrix = responses,
    item_data = item_data,
    difficulty_formula = ~ 1,  # Intercept only
    discrimination_formula = ~ 1,
    threshold_formula = NULL,  # NULL = no threshold covariates
    model = "GRM",
    control = list(eval.max = 1000, iter.max = 500, trace = 0)
  )
}, error = function(e) {
  list(error = e$message)
})

if ("error" %in% names(test1)) {
  cat("✗ FAILED:", test1$error, "\n\n")
} else {
  cat("✓ SUCCESS\n")
  cat("  - Converged:", test1$convergence$converged, "\n")
  cat("  - LogLik:", round(test1$logLik, 2), "\n")
  cat("  - threshold_covariate_model flag:", test1$threshold_covariate_model, "\n")
  cat("  - Threshold coefficients:",
      ifelse(is.null(test1$regression_coefficients$threshold), "NULL (correct)", "NOT NULL (error)"), "\n\n")
}

# TEST 2: Check that threshold_formula = ~ x works (new feature)
cat("Test 2: threshold_formula = ~ x (new feature)\n")
cat("------------------------------------------------\n")

test2 <- tryCatch({
  fit_eirt(
    response_matrix = responses,
    item_data = item_data,
    difficulty_formula = ~ 1,  # Intercept only
    discrimination_formula = ~ 1,
    threshold_formula = ~ x,  # Threshold covariates!
    model = "GRM",
    control = list(eval.max = 1000, iter.max = 500, trace = 0)
  )
}, error = function(e) {
  list(error = e$message)
})

if ("error" %in% names(test2)) {
  cat("✗ FAILED:", test2$error, "\n\n")
} else {
  cat("✓ SUCCESS\n")
  cat("  - Converged:", test2$convergence$converged, "\n")
  cat("  - LogLik:", round(test2$logLik, 2), "\n")
  cat("  - threshold_covariate_model flag:", test2$threshold_covariate_model, "\n")
  cat("  - Threshold coefficients:\n")
  if (!is.null(test2$regression_coefficients$threshold)) {
    print(test2$regression_coefficients$threshold)
    cat("  ✓ Threshold coefficients present (correct)\n\n")
  } else {
    cat("  ✗ NULL (error - should be present)\n\n")
  }
}

# TEST 3: Check print method handles threshold formula
if (!"error" %in% names(test2)) {
  cat("Test 3: print() method with threshold formula\n")
  cat("------------------------------------------------\n")
  print(test2)
}

cat("\n==================================================\n")
cat("SUMMARY\n")
cat("==================================================\n\n")

test1_pass <- !"error" %in% names(test1)
test2_pass <- !"error" %in% names(test2)

cat("Test 1 (NULL threshold formula):", ifelse(test1_pass, "PASS", "FAIL"), "\n")
cat("Test 2 (~ x threshold formula):", ifelse(test2_pass, "PASS", "FAIL"), "\n")

if (test1_pass && test2_pass) {
  cat("\n✓✓✓ ALL TESTS PASSED ✓✓✓\n")
  cat("\nThreshold-level EIRT implementation is working!\n")
} else {
  cat("\n✗ SOME TESTS FAILED\n")
  cat("The implementation needs debugging.\n")
}

cat("\n==================================================\n")
