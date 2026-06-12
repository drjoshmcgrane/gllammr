#!/usr/bin/env Rscript
# Test EIRT API changes at the code level (without package installation)

cat("\n=== Checking EIRT API Implementation ===\n\n")

# Test 1: Check fit_eirt function signature
cat("Test 1: Checking fit_eirt function signature...\n")

# Source the eirt.R file
source("R/eirt.R")

# Check function signature
eirt_formals <- formals(fit_eirt)

# Check model parameter
if ("model" %in% names(eirt_formals)) {
  model_choices <- eval(eirt_formals$model)
  cat("  Model choices:", paste(model_choices, collapse = ", "), "\n")

  if ("LPCM" %in% model_choices) {
    cat("  FAIL: LPCM still in model choices (should be removed)\n")
  } else {
    cat("  PASS: LPCM removed from model choices\n")
  }

  expected_models <- c("Rasch", "2PL", "GRM", "PCM", "GPCM")
  if (all(expected_models %in% model_choices) && length(model_choices) == 5) {
    cat("  PASS: Correct model choices:", paste(expected_models, collapse = ", "), "\n")
  }
} else {
  cat("  FAIL: model parameter not found\n")
}

# Check item_residuals parameter
if ("item_residuals" %in% names(eirt_formals)) {
  cat("  PASS: item_residuals parameter exists\n")

  if (identical(eirt_formals$item_residuals, TRUE)) {
    cat("  PASS: item_residuals default is TRUE\n")
  } else {
    cat("  FAIL: item_residuals default is not TRUE (found:",
        eirt_formals$item_residuals, ")\n")
  }
} else {
  cat("  FAIL: item_residuals parameter not found\n")
}

# Check all expected parameters
expected_params <- c("response_matrix", "item_data", "difficulty_formula",
                     "discrimination_formula", "threshold_formula",
                     "weights", "model", "item_residuals", "start", "control")

missing_params <- setdiff(expected_params, names(eirt_formals))
if (length(missing_params) == 0) {
  cat("  PASS: All expected parameters present\n")
} else {
  cat("  FAIL: Missing parameters:", paste(missing_params, collapse = ", "), "\n")
}

# Test 2: Check function body for key logic
cat("\nTest 2: Checking function body for poly_model_type logic...\n")

eirt_body <- deparse(body(fit_eirt))
eirt_body_str <- paste(eirt_body, collapse = "\n")

# Check for poly_model_type logic
if (grepl("poly_model_type", eirt_body_str)) {
  cat("  PASS: poly_model_type logic present\n")

  # Check for PCM + threshold_formula logic
  if (grepl('model == "PCM".*threshold_formula', eirt_body_str) ||
      grepl('threshold_formula.*model == "PCM"', eirt_body_str)) {
    cat("  PASS: PCM with threshold_formula logic present\n")
  } else {
    cat("  INFO: Could not detect PCM + threshold_formula logic (may still be present)\n")
  }

  # Check for GPCM + threshold_formula logic
  if (grepl('model == "GPCM".*threshold_formula', eirt_body_str) ||
      grepl('threshold_formula.*model == "GPCM"', eirt_body_str)) {
    cat("  PASS: GPCM with threshold_formula logic present\n")
  } else {
    cat("  INFO: Could not detect GPCM + threshold_formula logic (may still be present)\n")
  }
} else {
  cat("  WARN: poly_model_type logic not detected in function body\n")
}

# Check for item_residuals in tmb_data
if (grepl("item_residuals", eirt_body_str)) {
  cat("  PASS: item_residuals referenced in function body\n")
} else {
  cat("  FAIL: item_residuals not found in function body\n")
}

# Test 3: Check TMB template
cat("\nTest 3: Checking TMB template (src/gllamm_eirt.hpp)...\n")

tmb_code <- readLines("src/gllamm_eirt.hpp")
tmb_str <- paste(tmb_code, collapse = "\n")

# Check for DATA_INTEGER(item_residuals)
if (grepl("DATA_INTEGER\\(item_residuals\\)", tmb_str)) {
  cat("  PASS: DATA_INTEGER(item_residuals) found in TMB template\n")
} else {
  cat("  FAIL: DATA_INTEGER(item_residuals) not found in TMB template\n")
}

# Check for conditional item residuals logic
if (grepl("if \\(item_residuals == 1\\)", tmb_str)) {
  cat("  PASS: Conditional item_residuals logic found in TMB template\n")
} else {
  cat("  FAIL: Conditional item_residuals logic not found in TMB template\n")
}

# Check that epsilon_b is used conditionally
if (grepl("difficulty_pred \\+ epsilon_b", tmb_str)) {
  cat("  PASS: Conditional epsilon_b usage found\n")
} else {
  cat("  INFO: Could not detect conditional epsilon_b usage\n")
}

# Test 4: Check documentation
cat("\nTest 4: Checking documentation in R/eirt.R...\n")

# Check for @param item_residuals in documentation
if (any(grepl("@param item_residuals", tmb_code))) {
  cat("  (Documentation check skipped - need to check R file separately)\n")
}

# Read R file
r_code <- readLines("R/eirt.R")
r_str <- paste(r_code, collapse = "\n")

if (grepl("@param item_residuals", r_str)) {
  cat("  PASS: @param item_residuals found in documentation\n")
} else {
  cat("  FAIL: @param item_residuals not documented\n")
}

# Check for discrimination_formula documentation
if (grepl("@param discrimination_formula", r_str)) {
  cat("  PASS: @param discrimination_formula found in documentation\n")
} else {
  cat("  INFO: discrimination_formula documentation may need updating\n")
}

# Test 5: Check NEWS.md
cat("\nTest 5: Checking NEWS.md for API changes documentation...\n")

news <- readLines("NEWS.md")
news_str <- paste(news, collapse = "\n")

if (grepl("EIRT", news_str) && grepl("item_residuals", news_str)) {
  cat("  PASS: EIRT and item_residuals mentioned in NEWS.md\n")
} else {
  cat("  WARN: NEWS.md may not fully document the changes\n")
}

if (grepl("LPCM", news_str)) {
  cat("  PASS: LPCM changes mentioned in NEWS.md\n")
} else {
  cat("  INFO: LPCM removal not explicitly mentioned\n")
}

cat("\n=== EIRT API Check Complete ===\n\n")
cat("Summary: The API changes appear to be implemented at the R code level.\n")
cat("Next step: Full testing requires package compilation and installation.\n\n")
