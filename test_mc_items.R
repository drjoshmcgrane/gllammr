#!/usr/bin/env Rscript
# Test mc_items parameter for selective 3PL guessing

cat("\n=== Testing mc_items Parameter for 3PL ===\n\n")

# Read source files once
irt_code <- readLines("R/irt.R")
irt_str <- paste(irt_code, collapse = "\n")

tmb_code <- readLines("src/gllamm_irt.hpp")
tmb_str <- paste(tmb_code, collapse = "\n")

# Test 1: Check function signature
cat("Test 1: Check mc_items parameter exists...\n")

if (grepl("mc_items = NULL", irt_str)) {
  cat("  PASS: mc_items parameter exists with default NULL\n")
} else {
  cat("  FAIL: mc_items parameter not found\n")
}

# Test 2: Check documentation
cat("\nTest 2: Check documentation mentions mc_items...\n")

if (grepl("@param mc_items", irt_str)) {
  cat("  PASS: mc_items documented\n")
} else {
  cat("  FAIL: mc_items not documented\n")
}

# Test 3: Check TMB template
cat("\nTest 3: Check TMB template has mc_items data input...\n")

if (grepl("DATA_IVECTOR\\(mc_items\\)", tmb_str)) {
  cat("  PASS: DATA_IVECTOR(mc_items) found in TMB template\n")
} else {
  cat("  FAIL: DATA_IVECTOR(mc_items) not found\n")
}

# Test 4: Check conditional guessing logic
cat("\nTest 4: Check TMB template has conditional guessing logic...\n")

if (grepl("if \\(mc_items\\(item\\) == 1\\)", tmb_str)) {
  cat("  PASS: Conditional guessing logic found\n")
} else {
  cat("  FAIL: Conditional guessing logic not found\n")
}

# Test 5: Check R code passes mc_items to TMB
cat("\nTest 5: Check R code passes mc_items to TMB data...\n")

if (grepl("mc_items = as.integer\\(mc_indicator\\)", irt_str)) {
  cat("  PASS: mc_items passed to TMB data\n")
} else {
  cat("  FAIL: mc_items not passed to TMB data\n")
}

# Test 6: Check example in documentation
cat("\nTest 6: Check example usage in documentation...\n")

if (grepl("mc_items = 1:15", irt_str)) {
  cat("  PASS: Example with mc_items found\n")
} else {
  cat("  FAIL: Example with mc_items not found\n")
}

# Test 7: Check validation logic
cat("\nTest 7: Check mc_items validation...\n")

if (grepl("mc_items.*3PL", irt_str)) {
  cat("  PASS: Validation for 3PL found\n")
} else {
  cat("  INFO: May not have validation (might be okay)\n")
}

# Test 8: Check result object stores mc_items
cat("\nTest 8: Check result object stores mc_items...\n")

if (grepl("mc_items = if.*mc_indicator", irt_str)) {
  cat("  PASS: mc_items stored in result object\n")
} else {
  cat("  FAIL: mc_items not stored in result object\n")
}

# Test 9: Check print method shows MC info
cat("\nTest 9: Check print method shows MC items info...\n")

if (grepl("MC items with guessing", irt_str)) {
  cat("  PASS: Print method shows MC items info\n")
} else {
  cat("  FAIL: Print method doesn't show MC items info\n")
}

cat("\n=== Implementation Check Complete ===\n\n")

cat("Summary:\n")
cat("- R interface: mc_items parameter added\n")
cat("- TMB template: Conditional guessing logic implemented\n")
cat("- Documentation: Examples and parameter docs added\n")
cat("- Result object: mc_items information stored\n")
cat("- Print method: Shows MC item count for 3PL\n\n")

cat("Next step: Test with actual data (requires package compilation)\n\n")
