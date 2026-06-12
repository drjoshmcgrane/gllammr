# Test EIRT API Changes
# Tests for:
# 1. Pure LLTM (item_residuals = FALSE)
# 2. PCM with threshold_formula (former LPCM)
# 3. GPCM with threshold_formula (new capability)
# 4. Discrimination predictors for 2PL/GPCM

test_that("Pure LLTM works (item_residuals = FALSE)", {
  skip_if_not_installed("MASS")

  # Generate simple Rasch data
  set.seed(123)
  n_persons <- 100
  n_items <- 10

  # True abilities
  theta_true <- rnorm(n_persons, 0, 1)

  # Item covariates
  item_data <- data.frame(
    item_id = 1:n_items,
    difficulty_covar = rnorm(n_items, 0, 0.5)
  )

  # Generate difficulties from covariates only (pure LLTM)
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

  # Fit pure LLTM (no residuals)
  fit_pure <- fit_eirt(
    response_matrix = responses,
    item_data = item_data,
    difficulty_formula = ~ difficulty_covar,
    model = "Rasch",
    item_residuals = FALSE
  )

  expect_s3_class(fit_pure, "gllamm_eirt")
  expect_false(fit_pure$item_residuals)

  # Fit LLTM + error (with residuals)
  fit_error <- fit_eirt(
    response_matrix = responses,
    item_data = item_data,
    difficulty_formula = ~ difficulty_covar,
    model = "Rasch",
    item_residuals = TRUE
  )

  expect_s3_class(fit_error, "gllamm_eirt")
  expect_true(fit_error$item_residuals)

  # Pure LLTM should have fewer parameters
  expect_lt(length(fit_pure$tmb_obj$par), length(fit_error$tmb_obj$par))
})

test_that("PCM with threshold_formula works (former LPCM)", {
  skip("Requires polytomous data - implement when testing framework ready")

  # This will test:
  # fit_eirt(..., model = "PCM", threshold_formula = ~ x)
  # Should work identically to old model = "LPCM"
})

test_that("GPCM with threshold_formula works (new capability)", {
  skip("Requires polytomous data - implement when testing framework ready")

  # This will test:
  # fit_eirt(..., model = "GPCM", threshold_formula = ~ x)
  # This is a NEW capability (GPCM couldn't have threshold predictors before)
})

test_that("2PL with discrimination predictors works", {
  skip_if_not_installed("MASS")

  # Generate 2PL data
  set.seed(456)
  n_persons <- 100
  n_items <- 10

  theta_true <- rnorm(n_persons, 0, 1)

  # Item covariates
  item_data <- data.frame(
    item_id = 1:n_items,
    difficulty_covar = rnorm(n_items, 0, 0.5),
    discrimination_covar = rnorm(n_items, 0, 0.3)
  )

  # Generate parameters from covariates
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

  # Fit 2PL with discrimination predictors
  fit_2pl <- fit_eirt(
    response_matrix = responses,
    item_data = item_data,
    difficulty_formula = ~ difficulty_covar,
    discrimination_formula = ~ discrimination_covar,
    model = "2PL"
  )

  expect_s3_class(fit_2pl, "gllamm_eirt")
  expect_equal(fit_2pl$model, "2PL")

  # Should have both difficulty (gamma) and discrimination (delta)
  # regression coefficients
  gamma_hat <- fit_2pl$regression_coefficients$difficulty
  delta_hat <- fit_2pl$regression_coefficients$discrimination
  expect_false(is.null(gamma_hat))
  expect_false(is.null(delta_hat))
  expect_true(all(is.finite(gamma_hat)))
  expect_true(all(is.finite(delta_hat)))
  expect_true("difficulty_covar" %in% names(gamma_hat))
  expect_true("discrimination_covar" %in% names(delta_hat))
})

test_that("LPCM model name is rejected", {
  # Generate minimal data
  responses <- matrix(c(0, 1, 1, 0), 2, 2)
  item_data <- data.frame(item_id = 1:2, x = c(0, 1))

  # Should error because "LPCM" is no longer a valid model choice
  expect_error(
    fit_eirt(
      response_matrix = responses,
      item_data = item_data,
      model = "LPCM"
    ),
    "should be one of"
  )
})

test_that("item_residuals parameter is documented and accessible", {
  # Check that the parameter exists in function signature
  eirt_formals <- formals(fit_eirt)
  expect_true("item_residuals" %in% names(eirt_formals))
  expect_equal(eirt_formals$item_residuals, TRUE)  # Default should be TRUE
})
