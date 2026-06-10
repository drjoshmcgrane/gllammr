# Enhanced Ordinal Models Test Suite
# Tests for ACL, CRL, and PPO implementations

library(GLLAMMR)
dyn.load('/Users/josh/Documents/Claude_Code/GLLAMMR/src/gllamm_ordinal.so')

cat("========================================\n")
cat("Enhanced Ordinal Models Test Suite\n")
cat("========================================\n\n")

# Simulate test data
set.seed(42)
n <- 200
test_data <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  group = factor(rep(1:10, each = 20))
)

# Generate ordinal response with mild non-proportionality
beta1 <- c(0.5, 0.3, 0.2)  # Different effects per threshold
eta <- outer(test_data$x1, beta1, "*") +
       outer(test_data$x2, beta1 * 0.5, "*")
probs <- apply(eta, 1, function(e) {
  c(plogis(-1 - e[1]),
    plogis(0 - e[2]) - plogis(-1 - e[1]),
    plogis(1 - e[3]) - plogis(0 - e[2]),
    1 - plogis(1 - e[3]))
})
test_data$y <- factor(apply(probs, 2, function(p) sample(1:4, 1, prob = p)),
                      levels = 1:4, ordered = TRUE)

cat("Test data: N =", n, "| Categories = 4\n\n")

# Test 1: All link types fit successfully
cat("Test 1: Fitting all 6 link types\n")
cat("----------------------------------\n")

link_types <- c("logit", "probit", "acl", "crl_forward", "crl_backward", "ppo")
fits <- list()
test1_pass <- TRUE

for (link in link_types) {
  cat(sprintf("  %-15s: ", link))

  fit <- tryCatch({
    fit_ordinal(y ~ x1 + x2 + (1 | group), test_data, link = link)
  }, error = function(e) {
    cat("FAILED -", e$message, "\n")
    test1_pass <<- FALSE
    return(NULL)
  })

  if (!is.null(fit)) {
    if (fit$convergence$converged) {
      cat("✓ Converged | LogLik =", round(fit$logLik, 2), "\n")
      fits[[link]] <- fit
    } else {
      cat("✗ Did not converge\n")
      test1_pass <- FALSE
    }
  }
}

cat("\nTest 1:", ifelse(test1_pass, "PASS ✓", "FAIL ✗"), "\n\n")


# Test 2: PPO coefficients structure
cat("Test 2: PPO coefficient matrix structure\n")
cat("-----------------------------------------\n")

ppo_fit <- fits[["ppo"]]
test2_pass <- FALSE

if (!is.null(ppo_fit)) {
  if (!is.null(ppo_fit$coefficients$beta_ppo)) {
    beta_ppo <- ppo_fit$coefficients$beta_ppo
    expected_dims <- c(3, 3)  # 3 thresholds x 3 fixed effects (intercept + x1 + x2)

    if (all(dim(beta_ppo) == expected_dims)) {
      cat("  Matrix dimensions:", dim(beta_ppo)[1], "x", dim(beta_ppo)[2], "✓\n")
      cat("  Rownames:", paste(rownames(beta_ppo), collapse = ", "), "✓\n")
      cat("  Colnames:", paste(colnames(beta_ppo), collapse = ", "), "✓\n")
      test2_pass <- TRUE
    } else {
      cat("  ✗ Wrong dimensions:", dim(beta_ppo)[1], "x", dim(beta_ppo)[2],
          "(expected", expected_dims[1], "x", expected_dims[2], ")\n")
    }
  } else {
    cat("  ✗ beta_ppo not found in coefficients\n")
  }
}

cat("\nTest 2:", ifelse(test2_pass, "PASS ✓", "FAIL ✗"), "\n\n")


# Test 3: Proportional odds test
cat("Test 3: Proportional odds assumption test\n")
cat("------------------------------------------\n")

po_fit <- fits[["logit"]]
test3_pass <- FALSE

if (!is.null(po_fit)) {
  po_test <- tryCatch({
    test_proportional_odds(po_fit, data = test_data)
  }, error = function(e) {
    cat("  ✗ Error:", e$message, "\n")
    return(NULL)
  })

  if (!is.null(po_test)) {
    cat("  LRT statistic:", round(po_test$statistic, 3), "\n")
    cat("  Degrees of freedom:", po_test$df, "\n")
    cat("  p-value:", format.pval(po_test$p_value), "\n")
    cat("  Conclusion:", po_test$conclusion, "\n")

    # Check that PPO model has better (or equal) log-likelihood
    if (po_test$ppo_logLik >= po_test$base_logLik - 1e-6) {
      cat("  ✓ PPO logLik ≥ PO logLik\n")
      test3_pass <- TRUE
    } else {
      cat("  ✗ PPO logLik <", "PO logLik\n")
    }
  }
}

cat("\nTest 3:", ifelse(test3_pass, "PASS ✓", "FAIL ✗"), "\n\n")


# Test 4: Fit statistics
cat("Test 4: Model fit statistics\n")
cat("----------------------------\n")

test4_pass <- TRUE

for (link in c("logit", "ppo")) {
  cat("  ", link, ":\n")
  fit <- fits[[link]]

  if (!is.null(fit)) {
    fit_stats <- tryCatch({
      fit(fit, test_po = FALSE)
    }, error = function(e) {
      cat("    ✗ Error:", e$message, "\n")
      test4_pass <<- FALSE
      return(NULL)
    })

    if (!is.null(fit_stats)) {
      if (!is.null(fit_stats$pseudo_R2)) {
        cat("    Pseudo-R²:", round(fit_stats$pseudo_R2, 3), "✓\n")
      } else {
        cat("    ✗ Pseudo-R² missing\n")
        test4_pass <- FALSE
      }

      if (!is.null(fit_stats$AIC)) {
        cat("    AIC:", round(fit_stats$AIC, 2), "✓\n")
      } else {
        cat("    ✗ AIC missing\n")
        test4_pass <- FALSE
      }
    }
  }
}

cat("\nTest 4:", ifelse(test4_pass, "PASS ✓", "FAIL ✗"), "\n\n")


# Test 5: Plotting functions
cat("Test 5: Plotting functions\n")
cat("--------------------------\n")

test5_pass <- TRUE

for (link in c("logit", "ppo")) {
  cat("  ", link, ":\n")
  fit <- fits[[link]]

  if (!is.null(fit)) {
    # Suppress plotting output
    pdf(file = NULL)
    plot_result <- tryCatch({
      plot(fit, covariate = "x1", which = 1)
      cat("    plot() executed ✓\n")
      TRUE
    }, error = function(e) {
      cat("    ✗ Error:", e$message, "\n")
      FALSE
    })
    dev.off()

    test5_pass <- test5_pass && plot_result
  }
}

cat("\nTest 5:", ifelse(test5_pass, "PASS ✓", "FAIL ✗"), "\n\n")


# Test 6: Model comparison (different links should give different fits)
cat("Test 6: Model comparison across link types\n")
cat("-------------------------------------------\n")

logliks <- sapply(fits, function(f) f$logLik)
unique_logliks <- length(unique(round(logliks, 4)))

cat("  Unique log-likelihoods:", unique_logliks, "/", length(logliks), "\n")

if (unique_logliks >= 4) {
  cat("  ✓ Different links produce different fits\n")
  test6_pass <- TRUE
} else {
  cat("  ✗ Too many duplicate log-likelihoods\n")
  test6_pass <- FALSE
}

# Show AIC comparison
cat("\n  AIC comparison:\n")
aics <- sapply(fits, function(f) f$AIC)
sorted_aics <- sort(aics)
for (i in 1:length(sorted_aics)) {
  cat("    ", sprintf("%-15s: %.2f", names(sorted_aics)[i], sorted_aics[i]), "\n")
}

cat("\nTest 6:", ifelse(test6_pass, "PASS ✓", "FAIL ✗"), "\n\n")


# Overall summary
cat("========================================\n")
cat("Test Summary\n")
cat("========================================\n")

all_tests <- c(test1_pass, test2_pass, test3_pass, test4_pass, test5_pass, test6_pass)
tests_passed <- sum(all_tests)
tests_total <- length(all_tests)

cat("Tests passed:", tests_passed, "/", tests_total, "\n")

if (all(all_tests)) {
  cat("\n✅ ALL TESTS PASSED ✅\n")
  cat("\nEnhanced ordinal models (ACL, CRL, PPO) are fully functional!\n")
} else {
  cat("\n⚠️  SOME TESTS FAILED ⚠️\n")
  cat("Failed tests:", which(!all_tests), "\n")
}

cat("========================================\n")
