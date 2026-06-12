#!/usr/bin/env Rscript

# Manual test for marginal predictions
library(TMB)

# Load compiled models
dyn.load(dynlib("src/gllamm_gaussian"))
dyn.load(dynlib("src/gllamm_binomial"))
dyn.load(dynlib("src/gllamm_poisson"))

# Source R functions
source("R/marginal_utils.R")
source("R/predict.R")
source("R/gllamm.R")
source("R/formula_parser.R")
source("R/model_matrices.R")
source("R/tmb_interface_v2.R")

cat("===== Testing Marginal Predictions =====\n\n")

# Test 1: Gaussian - marginal should equal conditional
cat("Test 1: Gaussian model - marginal should equal fixed effects\n")
set.seed(123)
n_groups <- 10
n_per_group <- 10
n <- n_groups * n_per_group

data_gauss <- data.frame(
  y = rnorm(n, mean = 2, sd = 1),
  x = rnorm(n),
  group = factor(rep(1:n_groups, each = n_per_group))
)

tryCatch({
  fit_gauss <- gllamm(y ~ x + (1 | group), data = data_gauss, family = gaussian())

  pred_fixed <- predict(fit_gauss, re.form = NA, type = "response")
  pred_marg <- predict(fit_gauss, type = "marginal", n_sim = 1000)

  diff <- max(abs(pred_fixed - pred_marg))
  cat("  Max difference between fixed and marginal:", diff, "\n")
  cat("  Test 1:", ifelse(diff < 1e-6, "PASSED", "FAILED"), "\n\n")
}, error = function(e) {
  cat("  Test 1 ERROR:", e$message, "\n\n")
})


# Test 2: Binomial - marginal should differ from conditional
cat("Test 2: Binomial model - marginal predictions run successfully\n")
set.seed(456)
n_groups <- 12
n_per_group <- 15
n <- n_groups * n_per_group

# Create data with strong group effects
group_effects <- rnorm(n_groups, 0, 1)
x_vals <- rnorm(n)
group_ids <- rep(1:n_groups, each = n_per_group)

eta <- 0.5 * x_vals + group_effects[group_ids]
prob <- plogis(eta)
y_vals <- rbinom(n, 1, prob)

data_binom <- data.frame(
  y = y_vals,
  x = x_vals,
  group = factor(group_ids)
)

tryCatch({
  fit_binom <- gllamm(y ~ x + (1 | group), data = data_binom, family = binomial())

  # Conditional at u=0
  pred_cond <- predict(fit_binom, re.form = NA, type = "response")

  # Marginal
  pred_marg <- predict(fit_binom, type = "marginal", n_sim = 2000)

  cat("  Conditional mean:", round(mean(pred_cond), 4), "\n")
  cat("  Marginal mean:", round(mean(pred_marg), 4), "\n")
  cat("  Max difference:", round(max(abs(pred_cond - pred_marg)), 4), "\n")

  # Marginal should be different (usually more conservative)
  mean_diff <- abs(mean(pred_cond) - mean(pred_marg))
  cat("  Mean difference:", round(mean_diff, 4), "\n")

  # Test passes if marginal predictions run and are in valid range
  valid_range <- all(pred_marg >= 0 & pred_marg <= 1)
  cat("  All predictions in [0,1]:", valid_range, "\n")
  cat("  Test 2:", ifelse(valid_range, "PASSED", "FAILED"), "\n\n")
}, error = function(e) {
  cat("  Test 2 ERROR:", e$message, "\n\n")
})


# Test 3: Standard errors
cat("Test 3: Marginal predictions with standard errors\n")
tryCatch({
  result <- predict(fit_binom, type = "marginal", se.fit = TRUE, n_sim = 1000)

  cat("  Result is a list:", is.list(result), "\n")
  cat("  Has 'fit' component:", "fit" %in% names(result), "\n")
  cat("  Has 'se.fit' component:", "se.fit" %in% names(result), "\n")

  if (is.list(result) && "se.fit" %in% names(result)) {
    cat("  Mean SE:", round(mean(result$se.fit), 4), "\n")
    cat("  Range of SE:", round(min(result$se.fit), 4), "to", round(max(result$se.fit), 4), "\n")
    all_positive <- all(result$se.fit >= 0)
    cat("  All SE non-negative:", all_positive, "\n")
    cat("  Test 3:", ifelse(all_positive, "PASSED", "FAILED"), "\n\n")
  } else {
    cat("  Test 3: FAILED (wrong structure)\n\n")
  }
}, error = function(e) {
  cat("  Test 3 ERROR:", e$message, "\n\n")
})


# Test 4: Poisson model
cat("Test 4: Poisson model - marginal predictions\n")
set.seed(789)
n_groups <- 10
n_per_group <- 12
n <- n_groups * n_per_group

data_pois <- data.frame(
  y = rpois(n, lambda = 5),
  x = rnorm(n, 0, 0.5),
  group = factor(rep(1:n_groups, each = n_per_group))
)

tryCatch({
  fit_pois <- gllamm(y ~ x + (1 | group), data = data_pois, family = poisson())

  pred_marg <- predict(fit_pois, type = "marginal", n_sim = 1000)

  cat("  Mean marginal prediction:", round(mean(pred_marg), 3), "\n")
  cat("  Mean observed y:", round(mean(data_pois$y), 3), "\n")

  all_nonneg <- all(pred_marg >= 0)
  cat("  All predictions non-negative:", all_nonneg, "\n")
  cat("  Test 4:", ifelse(all_nonneg, "PASSED", "FAILED"), "\n\n")
}, error = function(e) {
  cat("  Test 4 ERROR:", e$message, "\n\n")
})


# Test 5: Newdata predictions
cat("Test 5: Marginal predictions with newdata\n")
tryCatch({
  newdata <- data.frame(
    x = c(-2, -1, 0, 1, 2),
    group = factor(c(1, 1, 1, 1, 1))  # Use first group
  )

  pred_new <- predict(fit_binom, newdata = newdata, type = "marginal", n_sim = 1000)

  cat("  Predictions for x = -2 to 2:\n")
  for (i in 1:5) {
    cat("    x =", newdata$x[i], "-> p =", round(pred_new[i], 4), "\n")
  }

  # Should show monotonic relationship with x (approximately)
  is_ordered <- all(diff(pred_new) > -0.1)  # Allow small violations due to MC noise
  cat("  Roughly monotonic with x:", is_ordered, "\n")
  cat("  Test 5:", ifelse(is_ordered, "PASSED", "FAILED"), "\n\n")
}, error = function(e) {
  cat("  Test 5 ERROR:", e$message, "\n\n")
})


# Test 6: Convergence with sample size
cat("Test 6: MC convergence - more samples = more stable\n")
tryCatch({
  set.seed(111)
  pred_500 <- predict(fit_binom, type = "marginal", n_sim = 500)

  set.seed(111)
  pred_5000 <- predict(fit_binom, type = "marginal", n_sim = 5000)

  # Run twice to check stability
  set.seed(222)
  pred_500b <- predict(fit_binom, type = "marginal", n_sim = 500)

  set.seed(222)
  pred_5000b <- predict(fit_binom, type = "marginal", n_sim = 5000)

  var_500 <- var(pred_500 - pred_500b)
  var_5000 <- var(pred_5000 - pred_5000b)

  cat("  Variance (n=500):", format(var_500, scientific = TRUE), "\n")
  cat("  Variance (n=5000):", format(var_5000, scientific = TRUE), "\n")
  cat("  Improvement ratio:", round(var_500 / var_5000, 2), "x\n")

  improved <- var_5000 < var_500
  cat("  More samples = less variance:", improved, "\n")
  cat("  Test 6:", ifelse(improved, "PASSED", "FAILED"), "\n\n")
}, error = function(e) {
  cat("  Test 6 ERROR:", e$message, "\n\n")
})


cat("===== All Manual Tests Complete =====\n")
