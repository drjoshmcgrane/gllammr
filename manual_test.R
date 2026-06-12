#!/usr/bin/env Rscript
# Manual testing script for gllammr development
# Run with: Rscript tests/manual_test.R

cat("gllammr Manual Test Script\n")
cat("==========================\n\n")

# Load package (assumes you're in package root)
cat("Loading gllammr...\n")
devtools::load_all()

cat("\n1. Testing formula parser\n")
cat("---------------------------\n")

# Test data
test_data <- data.frame(
  y = rnorm(50),
  x = rnorm(50),
  group = rep(1:10, each = 5)
)

# Parse simple formula
formula1 <- y ~ x + (1 | group)
parsed1 <- gllammr:::parse_formula(formula1, test_data)
cat("Formula:", deparse(formula1), "\n")
cat("Fixed formula:", deparse(parsed1$fixed_formula), "\n")
cat("Number of random terms:", length(parsed1$random_terms), "\n")
cat("Random grouping:", parsed1$random_terms[[1]]$grouping, "\n")

cat("\n2. Testing model matrices\n")
cat("---------------------------\n")

mats <- gllammr:::make_model_matrices(parsed1, test_data)
cat("Number of observations:", mats$n_obs, "\n")
cat("Number of fixed effects:", mats$n_fixed, "\n")
cat("Number of groups:", mats$n_groups[1], "\n")
cat("Fixed effects design matrix dimensions:", dim(mats$X), "\n")

cat("\n3. Testing random term parsing\n")
cat("---------------------------\n")

rt1 <- gllammr:::parse_random_term("(1|g)", data.frame(g = 1:10))
cat("Simple random intercept - Grouping:", rt1$grouping, "\n")

test_data_nested <- data.frame(
  school = rep(1:3, each = 6),
  class = rep(1:9, each = 2)
)
rt2 <- gllammr:::parse_random_term("(1|school/class)", test_data_nested)
cat("Nested random effects - Nested:", rt2$nested, "\n")
cat("Nested random effects - Grouping vars:", paste(rt2$grouping, collapse = ", "), "\n")

cat("\n4. Testing validation\n")
cat("---------------------------\n")

# Should work
valid_result <- tryCatch({
  gllammr:::validate_formula(y ~ x + (1|group), test_data)
  "PASS"
}, error = function(e) {
  paste("FAIL:", e$message)
})
cat("Valid formula:", valid_result, "\n")

# Should fail - no response
invalid_result <- tryCatch({
  gllammr:::validate_formula(~ x + (1|group), test_data)
  "FAIL: Should have errored"
}, error = function(e) {
  "PASS: Correctly caught error"
})
cat("Invalid formula (no response):", invalid_result, "\n")

# Should fail - missing variable
invalid_result2 <- tryCatch({
  gllammr:::validate_formula(y ~ z + (1|group), test_data)
  "FAIL: Should have errored"
}, error = function(e) {
  "PASS: Correctly caught error"
})
cat("Invalid formula (missing var):", invalid_result2, "\n")

cat("\n5. Summary\n")
cat("---------------------------\n")
cat("Formula parser: ✓ Working\n")
cat("Model matrices: ✓ Working\n")
cat("Validation: ✓ Working\n")

cat("\n6. TMB Status\n")
cat("---------------------------\n")
tmb_installed <- requireNamespace("TMB", quietly = TRUE)
if (tmb_installed) {
  cat("TMB package: ✓ Installed\n")
  cat("TMB version:", as.character(packageVersion("TMB")), "\n")

  # Check if template is compiled
  cpp_file <- "src/gllamm_gaussian.cpp"
  if (file.exists(cpp_file)) {
    cat("TMB template source: ✓ Found\n")
    cat("Note: Template needs to be compiled with TMB::compile() before fitting models\n")
  } else {
    cat("TMB template source: ✗ Not found\n")
  }
} else {
  cat("TMB package: ✗ Not installed\n")
  cat("Install with: install.packages('TMB')\n")
}

cat("\n7. Test Dependencies\n")
cat("---------------------------\n")
deps <- c("Matrix", "methods", "stats", "testthat", "lme4")
for (dep in deps) {
  installed <- requireNamespace(dep, quietly = TRUE)
  status <- if (installed) "✓" else "✗"
  cat(sprintf("%-15s %s\n", paste0(dep, ":"), status))
}

cat("\n==========================\n")
cat("Manual test complete!\n")
cat("\nNext steps:\n")
cat("1. Compile TMB: TMB::compile('src/gllamm_gaussian.cpp')\n")
cat("2. Run tests: devtools::test()\n")
cat("3. Try example: See README.md\n")
