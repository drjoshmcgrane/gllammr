#!/usr/bin/env Rscript
# Check for response coding issues in IRT implementations

cat("\n=== Checking Response Coding Issues ===\n\n")

# Read all IRT-related files
irt_code <- readLines("R/irt.R")
eirt_code <- readLines("R/eirt.R")
tmb_irt <- readLines("src/gllamm_irt.hpp")
tmb_irt_poly <- readLines("src/gllamm_irt_poly.hpp")
tmb_eirt <- readLines("src/gllamm_eirt.hpp")

# Check 1: Dichotomous coding expectations
cat("Check 1: Dichotomous IRT coding expectations...\n")
irt_str <- paste(irt_code, collapse = "\n")

if (grepl("require binary responses \\(0/1\\)", irt_str)) {
  cat("  ✓ Dichotomous IRT expects 0/1 coding\n")
} else {
  cat("  ? Dichotomous coding expectation unclear\n")
}

# Check 2: Polytomous coding expectations
cat("\nCheck 2: Polytomous IRT coding expectations...\n")
tmb_poly_str <- paste(tmb_irt_poly, collapse = "\n")

if (grepl("Item responses \\(1, 2, \\.\\.\\., K\\)", tmb_poly_str)) {
  cat("  ✓ Polytomous IRT expects 1/2/.../K coding (comment found)\n")
} else {
  cat("  ? Polytomous coding expectation not documented\n")
}

if (grepl("CppAD::Integer\\(y\\(i\\)\\) - 1", tmb_poly_str)) {
  cat("  ✓ TMB converts to 0-indexed (subtracts 1)\n")
} else {
  cat("  ? TMB indexing unclear\n")
}

# Check 3: Polytomous validation in R
cat("\nCheck 3: Standard polytomous IRT validation...\n")

if (grepl("min\\(item_vals\\) != 1", irt_str)) {
  cat("  ✓ Validates responses start at 1\n")
} else {
  cat("  ✗ MISSING: No validation for 1-based coding\n")
}

if (grepl("max\\(item_vals\\) != n_categories_per_item", irt_str)) {
  cat("  ✓ Validates responses go to K\n")
} else {
  cat("  ✗ MISSING: No validation for max category\n")
}

# Check 4: EIRT validation
cat("\nCheck 4: EIRT validation...\n")
eirt_str <- paste(eirt_code, collapse = "\n")

if (grepl("min\\(.*item", eirt_str)) {
  cat("  ✓ EIRT has response validation\n")
} else {
  cat("  ✗ WARNING: EIRT has NO response validation!\n")
}

# Check 5: EIRT TMB template
cat("\nCheck 5: EIRT TMB template...\n")
tmb_eirt_str <- paste(tmb_eirt, collapse = "\n")

if (grepl("CppAD::Integer\\(y\\(i\\)\\) - 1", tmb_eirt_str)) {
  cat("  ✓ EIRT TMB also converts to 0-indexed\n")
  cat("  ⚠️  This means EIRT also expects 1-based coding for polytomous!\n")
} else {
  cat("  ? EIRT TMB indexing unclear\n")
}

# Check 6: Mixed items issue
cat("\n=== IDENTIFIED ISSUES ===\n\n")

cat("Issue 1: Coding incompatibility for mixed items\n")
cat("  Problem: Dichotomous (0/1) vs Polytomous (1/2/.../K) coding\n")
cat("  Impact: Cannot mix binary and polytomous items without recoding\n")
cat("  Example that FAILS:\n")
cat("    cbind(\n")
cat("      matrix(sample(0:1, 100, TRUE), 50, 2),  # Binary: 0/1\n")
cat("      matrix(sample(1:4, 150, TRUE), 50, 3)   # 4-cat: 1/2/3/4\n")
cat("    )\n")
cat("  Error: Binary items will fail validation (min != 1)\n\n")

cat("Issue 2: EIRT has no response validation\n")
cat("  Problem: fit_eirt() doesn't validate response coding\n")
cat("  Impact: Silent failures or wrong results\n")
cat("  Risk: Users could pass 0/1 data to polytomous EIRT\n\n")

cat("Issue 3: Documentation unclear\n")
cat("  Problem: Help files don't explain coding requirements clearly\n")
cat("  Impact: User confusion about how to code responses\n\n")

cat("\n=== RECOMMENDED FIXES ===\n\n")

cat("Fix 1: Add response validation to EIRT\n")
cat("  - Check that polytomous responses are 1-based\n")
cat("  - Provide clear error messages\n\n")

cat("Fix 2: Add auto-detection and recoding for mixed items\n")
cat("  - Detect if items are 0/1 coded\n")
cat("  - Auto-recode to 1/2 for polytomous models\n")
cat("  - Warn user about recoding\n\n")

cat("Fix 3: Update documentation\n")
cat("  - Clearly state coding requirements\n")
cat("  - Add examples with proper coding\n")
cat("  - Document mixed items workflow\n\n")

cat("Fix 4: Update mixed items example in NEWS.md\n")
cat("  - Current example uses 0/1 for binary items (WRONG!)\n")
cat("  - Should use 1/2 for binary items in polytomous context\n\n")

cat("=== Test Cases ===\n\n")

cat("Test 1: Standard polytomous (should work)\n")
cat("  responses <- matrix(sample(1:4, 200, TRUE), 50, 4)\n")
cat("  fit_irt(responses, model = 'PCM')  # OK\n\n")

cat("Test 2: Binary as polytomous with 1/2 coding (should work)\n")
cat("  responses <- matrix(sample(1:2, 200, TRUE), 50, 4)\n")
cat("  fit_irt(responses, model = 'PCM')  # OK (PCM with K=2)\n\n")

cat("Test 3: Binary with 0/1 coding (FAILS!)\n")
cat("  responses <- matrix(sample(0:1, 200, TRUE), 50, 4)\n")
cat("  fit_irt(responses, model = 'PCM')  # ERROR: min != 1\n\n")

cat("Test 4: Mixed items with 0/1 binary (FAILS!)\n")
cat("  responses <- cbind(\n")
cat("    matrix(sample(0:1, 100, TRUE), 50, 2),\n")
cat("    matrix(sample(1:4, 150, TRUE), 50, 3)\n")
cat("  )\n")
cat("  fit_irt(responses, model = 'PCM')  # ERROR on items 1-2\n\n")

cat("Test 5: Mixed items with 1/2 binary (SHOULD work)\n")
cat("  responses <- cbind(\n")
cat("    matrix(sample(1:2, 100, TRUE), 50, 2),  # Recoded!\n")
cat("    matrix(sample(1:4, 150, TRUE), 50, 3)\n")
cat("  )\n")
cat("  fit_irt(responses, model = 'PCM')  # Should work\n\n")

cat("\n=== Check Complete ===\n")
