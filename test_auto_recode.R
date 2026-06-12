#!/usr/bin/env Rscript
# Test auto-recoding for polytomous IRT

cat("\n=== Testing Auto-Recoding for Polytomous IRT ===\n\n")

# Read source files
irt_code <- readLines("R/irt.R")
eirt_code <- readLines("R/eirt.R")

irt_str <- paste(irt_code, collapse = "\n")
eirt_str <- paste(eirt_code, collapse = "\n")

# Test 1: Check for auto-recoding logic in fit_irt_polytomous
cat("Test 1: fit_irt_polytomous has auto-recoding...\n")

if (grepl("Auto-recoding.*binary.*from 0/1 to 1/2", irt_str)) {
  cat("  ✓ PASS: Auto-recoding message found\n")
} else {
  cat("  ✗ FAIL: Auto-recoding not implemented\n")
}

if (grepl("needs_recoding", irt_str)) {
  cat("  ✓ PASS: needs_recoding flag found\n")
} else {
  cat("  ✗ FAIL: needs_recoding flag not found\n")
}

if (grepl("min_val == 0 && max_val == n_cats - 1", irt_str)) {
  cat("  ✓ PASS: 0-based detection logic found\n")
} else {
  cat("  ✗ FAIL: 0-based detection not found\n")
}

if (grepl("response_matrix\\[, j\\] <- response_matrix\\[, j\\] \\+ 1", irt_str)) {
  cat("  ✓ PASS: Recoding (+1) logic found\n")
} else {
  cat("  ✗ FAIL: Recoding logic not found\n")
}

# Test 2: Check for auto-recoding in EIRT
cat("\nTest 2: fit_eirt has auto-recoding...\n")

if (grepl("Auto-recoding.*binary.*from 0/1 to 1/2", eirt_str)) {
  cat("  ✓ PASS: Auto-recoding message found in EIRT\n")
} else {
  cat("  ✗ FAIL: Auto-recoding not implemented in EIRT\n")
}

if (grepl("needs_recoding", eirt_str)) {
  cat("  ✓ PASS: needs_recoding flag found in EIRT\n")
} else {
  cat("  ✗ FAIL: needs_recoding flag not found in EIRT\n")
}

# Test 3: Check validation for invalid coding
cat("\nTest 3: Validation for invalid coding...\n")

if (grepl("has invalid response coding", irt_str)) {
  cat("  ✓ PASS: Invalid coding error message found\n")
} else {
  cat("  ✗ FAIL: No error for invalid coding\n")
}

if (grepl("Polytomous models require", irt_str)) {
  cat("  ✓ PASS: Helpful error message about requirements\n")
} else {
  cat("  ✗ FAIL: Error message not helpful\n")
}

# Test 4: Check NEWS.md updated
cat("\nTest 4: NEWS.md documents auto-recoding...\n")

news <- readLines("NEWS.md")
news_str <- paste(news, collapse = "\n")

if (grepl("Auto-recoding", news_str)) {
  cat("  ✓ PASS: Auto-recoding mentioned in NEWS.md\n")
} else {
  cat("  ✗ FAIL: Auto-recoding not documented\n")
}

if (grepl("auto-recoded from 0/1 to 1/2", news_str)) {
  cat("  ✓ PASS: Recoding behavior explained\n")
} else {
  cat("  ✗ FAIL: Recoding behavior not explained\n")
}

cat("\n=== Expected Behavior ===\n\n")

cat("Scenario 1: All items properly coded (1-based)\n")
cat("  Input:  matrix(sample(1:4, 200, TRUE), 50, 4)\n")
cat("  Result: No recoding, no message\n\n")

cat("Scenario 2: Binary items coded 0/1 with polytomous\n")
cat("  Input:  cbind(matrix(sample(0:1, 100, TRUE), 50, 2),\n")
cat("                matrix(sample(1:4, 150, TRUE), 50, 3))\n")
cat("  Result: Items 1-2 auto-recoded to 1/2\n")
cat("          Message: 'Auto-recoding 2 binary item(s)...'\n\n")

cat("Scenario 3: Invalid coding (e.g., 0/1/2 for 3 categories)\n")
cat("  Input:  matrix(sample(0:2, 200, TRUE), 50, 4)\n")
cat("  Result: ERROR - invalid coding, not sequential from 1\n\n")

cat("Scenario 4: All binary coded 0/1\n")
cat("  Input:  matrix(sample(0:1, 200, TRUE), 50, 4)\n")
cat("  Result: All 4 items recoded to 1/2\n")
cat("          Message displayed\n")
cat("          PCM treats as 2-category items\n\n")

cat("\n=== Test Complete ===\n")
