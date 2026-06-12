#!/usr/bin/env Rscript
# Manual test script for multinomial and survival marginal predictions

library(MASS)

# Try to load package, otherwise source files
if (!require(GLLAMMR, quietly = TRUE)) {
  source("R/formula.R")
  source("R/families.R")
  source("R/tmb_interface_v2.R")
  source("R/multinomial.R")
  source("R/survival.R")
  source("R/predict_multinomial.R")
  source("R/predict_survival.R")
  source("R/marginal_utils.R")

  # Load TMB templates if needed
  for (tmb_file in c("src/gllamm_multinomial.so", "src/gllamm_survival.so")) {
    if (file.exists(tmb_file)) {
      dyn.load(tmb_file)
    }
  }
}

set.seed(789)

cat("=== Testing Multinomial Marginal Predictions ===\n\n")

# Test 1: Multinomial Logit with Random Intercept
cat("Test 1: Multinomial Logit with Random Intercept\n")
cat("-------------------------------------------------\n")

n_groups <- 15
n_per_group <- 20
n_obs <- n_groups * n_per_group

mult_data <- data.frame(
  y = sample(0:2, n_obs, replace = TRUE, prob = c(0.4, 0.35, 0.25)),
  x1 = rnorm(n_obs),
  x2 = rnorm(n_obs),
  group = rep(1:n_groups, each = n_per_group)
)

tryCatch({
  fit_mult <- fit_multinomial(
    y ~ x1 + x2 + (1 | group),
    data = mult_data
  )

  cat("✓ Multinomial model fitted successfully\n")
  cat("  Number of categories:", fit_mult$n_categories, "\n")
  cat("  Beta matrix dimensions:", paste(dim(fit_mult$coefficients$beta_matrix), collapse=" x "), "\n\n")

  # Test conditional predictions - probabilities
  cat("Testing conditional probability predictions...\n")
  pred_cond <- predict(fit_mult, type = "probs")
  cat("✓ Conditional probabilities computed\n")
  cat("  Dimensions:", paste(dim(pred_cond), collapse=" x "), "\n")
  cat("  First obs probs:", round(pred_cond[1,], 3), "\n")
  cat("  Sum to 1?", all(abs(rowSums(pred_cond) - 1) < 1e-6), "\n\n")

  # Test conditional predictions - class
  cat("Testing class predictions...\n")
  pred_class <- predict(fit_mult, type = "class")
  cat("✓ Class predictions computed\n")
  cat("  Class distribution:", table(pred_class), "\n\n")

  # Test marginal predictions
  cat("Testing marginal predictions...\n")
  pred_marg <- predict(fit_mult, type = "marginal", n_sim = 500)
  cat("✓ Marginal probabilities computed\n")
  cat("  Dimensions:", paste(dim(pred_marg), collapse=" x "), "\n")
  cat("  First obs probs:", round(pred_marg[1,], 3), "\n")
  cat("  Sum to 1?", all(abs(rowSums(pred_marg) - 1) < 1e-6), "\n\n")

  # Compare conditional vs marginal
  cat("Comparing conditional vs marginal predictions:\n")
  cat("  Mean absolute difference:", round(mean(abs(pred_cond - pred_marg)), 4), "\n")
  cat("  Expected: marginal should differ from conditional\n\n")

  # Compare with empirical proportions
  cat("Comparing marginal predictions to empirical proportions:\n")
  emp_props <- table(mult_data$y) / length(mult_data$y)
  marg_avg <- colMeans(pred_marg)
  for (k in 0:(fit_mult$n_categories - 1)) {
    cat(sprintf("  Category %d: Predicted = %.3f, Empirical = %.3f\n",
                k, marg_avg[k+1], emp_props[as.character(k)]))
  }

  cat("\n✓ Test 1 PASSED\n\n")

}, error = function(e) {
  cat("✗ Test 1 FAILED:", conditionMessage(e), "\n\n")
})


# Test 2: Multinomial with New Data
cat("Test 2: Multinomial Predictions on New Data\n")
cat("---------------------------------------------\n")

new_data_mult <- data.frame(
  x1 = c(-2, 0, 2),
  x2 = c(-1, 0, 1),
  group = c(1, 1, 1)
)

tryCatch({
  pred_new_cond <- predict(fit_mult, newdata = new_data_mult, type = "probs")
  cat("✓ Conditional predictions on new data\n")
  cat("  Dimensions:", paste(dim(pred_new_cond), collapse=" x "), "\n")

  pred_new_marg <- predict(fit_mult, newdata = new_data_mult, type = "marginal", n_sim = 500)
  cat("✓ Marginal predictions on new data\n")
  cat("  Dimensions:", paste(dim(pred_new_marg), collapse=" x "), "\n")

  cat("\nPredictions across covariate range:\n")
  for (i in 1:nrow(new_data_mult)) {
    cat(sprintf("  x1=%.1f, x2=%.1f: P(Y=0) = %.3f, P(Y=1) = %.3f, P(Y=2) = %.3f\n",
                new_data_mult$x1[i], new_data_mult$x2[i],
                pred_new_marg[i, 1], pred_new_marg[i, 2], pred_new_marg[i, 3]))
  }

  cat("\n✓ Test 2 PASSED\n\n")

}, error = function(e) {
  cat("✗ Test 2 FAILED:", conditionMessage(e), "\n\n")
})


cat("\n=== Testing Survival Marginal Predictions ===\n\n")

# Test 3: Exponential Survival Model
cat("Test 3: Exponential Survival with Random Intercept\n")
cat("----------------------------------------------------\n")

n_groups_surv <- 10
n_per_group_surv <- 30
n_obs_surv <- n_groups_surv * n_per_group_surv

surv_data <- data.frame(
  time = rexp(n_obs_surv, rate = 0.1),
  status = rbinom(n_obs_surv, 1, 0.7),
  x1 = rnorm(n_obs_surv),
  group = rep(1:n_groups_surv, each = n_per_group_surv)
)

tryCatch({
  fit_surv_exp <- fit_survival(
    time ~ x1 + (1 | group),
    data = surv_data,
    event = surv_data$status,
    distribution = "exponential"
  )

  cat("✓ Exponential survival model fitted successfully\n")
  cat("  Fixed effects:", round(fit_surv_exp$coefficients$fixed, 3), "\n")
  cat("  Random SD:", round(sqrt(fit_surv_exp$coefficients$random_var[[1]]), 3), "\n\n")

  # Test linear predictor
  cat("Testing linear predictor predictions...\n")
  pred_lp <- predict(fit_surv_exp, type = "lp")
  cat("✓ Linear predictor computed\n")
  cat("  Length:", length(pred_lp), "\n")
  cat("  Range:", round(range(pred_lp), 3), "\n\n")

  # Test risk predictions
  cat("Testing relative risk predictions...\n")
  pred_risk <- predict(fit_surv_exp, type = "risk")
  cat("✓ Relative risk computed\n")
  cat("  Mean risk:", round(mean(pred_risk), 3), "\n\n")

  # Test conditional survival predictions
  times <- c(5, 10, 15, 20)
  cat("Testing conditional survival predictions at times:", times, "\n")
  pred_surv_cond <- predict(fit_surv_exp, type = "survival", times = times)
  cat("✓ Conditional survival computed\n")
  cat("  Dimensions:", paste(dim(pred_surv_cond), collapse=" x "), "\n")
  cat("  First obs survival at times:", round(pred_surv_cond[1,], 3), "\n\n")

  # Test marginal survival predictions
  cat("Testing marginal survival predictions...\n")
  pred_surv_marg <- predict(fit_surv_exp, type = "marginal_survival", times = times, n_sim = 500)
  cat("✓ Marginal survival computed\n")
  cat("  Dimensions:", paste(dim(pred_surv_marg), collapse=" x "), "\n")
  cat("  First obs survival at times:", round(pred_surv_marg[1,], 3), "\n\n")

  # Compare conditional vs marginal
  cat("Comparing conditional vs marginal survival:\n")
  cat("  Difference at t=5:", round(mean(pred_surv_cond[,1] - pred_surv_marg[,1]), 4), "\n")
  cat("  Difference at t=20:", round(mean(pred_surv_cond[,4] - pred_surv_marg[,4]), 4), "\n\n")

  # Test hazard predictions
  cat("Testing marginal hazard predictions...\n")
  pred_haz_marg <- predict(fit_surv_exp, type = "marginal_hazard", times = times, n_sim = 500)
  cat("✓ Marginal hazard computed\n")
  cat("  First obs hazard at times:", round(pred_haz_marg[1,], 3), "\n")
  cat("  Note: For exponential, hazard should be approximately constant\n\n")

  cat("✓ Test 3 PASSED\n\n")

}, error = function(e) {
  cat("✗ Test 3 FAILED:", conditionMessage(e), "\n\n")
})


# Test 4: Weibull Survival Model
cat("Test 4: Weibull Survival with Random Intercept\n")
cat("------------------------------------------------\n")

tryCatch({
  fit_surv_weib <- fit_survival(
    time ~ x1 + (1 | group),
    data = surv_data,
    event = surv_data$status,
    distribution = "Weibull"
  )

  cat("✓ Weibull survival model fitted successfully\n")
  cat("  Shape parameter:", round(fit_surv_weib$shape_parameter, 3), "\n\n")

  # Marginal survival
  pred_surv_weib <- predict(fit_surv_weib, type = "marginal_survival",
                            times = times, n_sim = 500)
  cat("✓ Marginal survival computed\n")
  cat("  First obs survival at times:", round(pred_surv_weib[1,], 3), "\n\n")

  # Marginal hazard (should not be constant for Weibull)
  pred_haz_weib <- predict(fit_surv_weib, type = "marginal_hazard",
                           times = times, n_sim = 500)
  cat("✓ Marginal hazard computed\n")
  cat("  First obs hazard at times:", round(pred_haz_weib[1,], 3), "\n")
  cat("  Note: For Weibull with shape ≠ 1, hazard changes over time\n\n")

  cat("✓ Test 4 PASSED\n\n")

}, error = function(e) {
  cat("✗ Test 4 FAILED:", conditionMessage(e), "\n\n")
})


# Test 5: Survival Predictions on New Data
cat("Test 5: Survival Predictions on New Data\n")
cat("------------------------------------------\n")

new_data_surv <- data.frame(
  x1 = c(-1, 0, 1),
  group = c(1, 1, 1)
)

tryCatch({
  pred_new_surv <- predict(fit_surv_exp, newdata = new_data_surv,
                          type = "marginal_survival", times = c(10, 20), n_sim = 500)
  cat("✓ Marginal survival predictions on new data\n")
  cat("  Dimensions:", paste(dim(pred_new_surv), collapse=" x "), "\n")

  cat("\nSurvival predictions across covariate range:\n")
  for (i in 1:nrow(new_data_surv)) {
    cat(sprintf("  x1=%.1f: S(10) = %.3f, S(20) = %.3f\n",
                new_data_surv$x1[i], pred_new_surv[i, 1], pred_new_surv[i, 2]))
  }

  cat("\n✓ Test 5 PASSED\n\n")

}, error = function(e) {
  cat("✗ Test 5 FAILED:", conditionMessage(e), "\n\n")
})


cat("=== Multinomial & Survival Marginal Predictions Testing Complete ===\n")
