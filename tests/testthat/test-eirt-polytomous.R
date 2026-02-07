test_that("EIRT accepts polytomous data with GRM", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(123)
  n_persons <- 100
  n_items <- 10
  n_categories <- 4

  # Create item covariates
  item_data <- data.frame(
    complexity = rnorm(n_items, 0, 1)
  )

  # Simulate GRM responses
  theta <- rnorm(n_persons, 0, 1)
  discrimination <- exp(rnorm(n_items, 0, 0.3))

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      # Simple threshold structure
      thresholds <- c(-1.5, -0.5, 0.5)

      probs <- numeric(n_categories)
      probs[1] <- plogis(discrimination[j] * (theta[i] - thresholds[1]))

      for (k in 2:(n_categories - 1)) {
        p_le_k <- plogis(discrimination[j] * (theta[i] - thresholds[k]))
        p_le_k_minus_1 <- plogis(discrimination[j] * (theta[i] - thresholds[k-1]))
        probs[k] <- p_le_k - p_le_k_minus_1
      }

      probs[n_categories] <- 1 - plogis(discrimination[j] * (theta[i] - thresholds[n_categories-1]))

      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Fit EIRT with GRM
  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ complexity,
    discrimination_formula = ~ 1,
    model = "GRM"
  )

  expect_s3_class(fit, "gllamm_eirt")
  expect_equal(fit$model, "GRM")
  expect_equal(fit$n_persons, n_persons)
  expect_equal(fit$n_items, n_items)
})


test_that("EIRT recovers item covariate effects for GRM", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(456)
  n_persons <- 150
  n_items <- 15
  n_categories <- 5

  # Create item covariate
  item_data <- data.frame(
    item_difficulty_predictor = rnorm(n_items, 0, 1)
  )

  # True parameters
  gamma_0 <- 0.0
  gamma_1 <- 0.6  # Positive effect makes items harder

  theta <- rnorm(n_persons, 0, 1)
  discrimination <- rep(1.5, n_items)

  # Generate base difficulties from covariates
  base_difficulty <- gamma_0 + gamma_1 * item_data$item_difficulty_predictor + rnorm(n_items, 0, 0.2)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      # Thresholds centered around base difficulty
      thresholds <- base_difficulty[j] + seq(-2, 2, length.out = n_categories - 1)

      probs <- numeric(n_categories)
      probs[1] <- plogis(discrimination[j] * (theta[i] - thresholds[1]))

      for (k in 2:(n_categories - 1)) {
        p_le_k <- plogis(discrimination[j] * (theta[i] - thresholds[k]))
        p_le_k_minus_1 <- plogis(discrimination[j] * (theta[i] - thresholds[k-1]))
        probs[k] <- p_le_k - p_le_k_minus_1
      }

      probs[n_categories] <- 1 - plogis(discrimination[j] * (theta[i] - thresholds[n_categories-1]))

      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Fit EIRT
  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ item_difficulty_predictor,
    discrimination_formula = ~ 1,
    model = "GRM"
  )

  gamma_hat <- fit$regression_coefficients$difficulty

  # Check recovery (with reasonable tolerance for polytomous)
  expect_equal(gamma_hat[["item_difficulty_predictor"]], gamma_1, tolerance = 0.4)
})


test_that("EIRT with PCM model", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(789)
  n_persons <- 100
  n_items <- 10
  n_categories <- 4

  # Item covariates
  item_data <- data.frame(
    item_type = factor(rep(c("TypeA", "TypeB"), each = 5))
  )

  # Simulate PCM responses (discrimination = 1 for all items)
  theta <- rnorm(n_persons, 0, 1)
  discrimination <- rep(1, n_items)  # Fixed for PCM

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      thresholds <- sort(rnorm(n_categories - 1, 0, 1))

      probs <- numeric(n_categories)
      probs[1] <- plogis(discrimination[j] * (theta[i] - thresholds[1]))

      for (k in 2:(n_categories - 1)) {
        p_le_k <- plogis(discrimination[j] * (theta[i] - thresholds[k]))
        p_le_k_minus_1 <- plogis(discrimination[j] * (theta[i] - thresholds[k-1]))
        probs[k] <- p_le_k - p_le_k_minus_1
      }

      probs[n_categories] <- 1 - plogis(discrimination[j] * (theta[i] - thresholds[n_categories-1]))

      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Fit EIRT with PCM
  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ item_type,
    discrimination_formula = ~ 1,
    model = "PCM"
  )

  expect_s3_class(fit, "gllamm_eirt")
  expect_equal(fit$model, "PCM")

  # Check that discrimination coefficients exist
  delta_hat <- fit$regression_coefficients$discrimination
  expect_equal(length(delta_hat), 1)
})


test_that("EIRT with GPCM model", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(111)
  n_persons <- 120
  n_items <- 12
  n_categories <- 3

  # Item covariates
  item_data <- data.frame(
    complexity = rnorm(n_items, 0, 1)
  )

  # Simulate GPCM responses
  theta <- rnorm(n_persons, 0, 1)
  discrimination <- exp(0.3 * item_data$complexity + rnorm(n_items, 0, 0.2))

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      thresholds <- c(-0.5, 0.5)

      probs <- numeric(n_categories)
      probs[1] <- plogis(discrimination[j] * (theta[i] - thresholds[1]))
      probs[2] <- plogis(discrimination[j] * (theta[i] - thresholds[2])) -
                  plogis(discrimination[j] * (theta[i] - thresholds[1]))
      probs[3] <- 1 - plogis(discrimination[j] * (theta[i] - thresholds[2]))

      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Fit EIRT with GPCM
  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ 1,
    discrimination_formula = ~ complexity,
    model = "GPCM"
  )

  expect_s3_class(fit, "gllamm_eirt")
  expect_equal(fit$model, "GPCM")

  # Check discrimination regression coefficients
  delta_hat <- fit$regression_coefficients$discrimination
  expect_equal(length(delta_hat), 2)
  expect_named(delta_hat, c("(Intercept)", "complexity"))
})


test_that("EIRT with mixed number of categories", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(222)
  n_persons <- 100
  n_items <- 8

  # Item covariates
  item_data <- data.frame(
    type = factor(rep(c("Short", "Long"), each = 4))
  )

  theta <- rnorm(n_persons, 0, 1)

  # Create mixed category responses
  responses <- matrix(NA, n_persons, n_items)

  # Items 1-4: 3 categories
  for (i in 1:n_persons) {
    for (j in 1:4) {
      probs <- c(0.3, 0.4, 0.3)
      responses[i, j] <- sample(1:3, 1, prob = probs)
    }
  }

  # Items 5-8: 5 categories
  for (i in 1:n_persons) {
    for (j in 5:8) {
      probs <- c(0.1, 0.2, 0.4, 0.2, 0.1)
      responses[i, j] <- sample(1:5, 1, prob = probs)
    }
  }

  # Fit EIRT with mixed categories
  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ type,
    discrimination_formula = ~ 1,
    model = "GRM"
  )

  expect_s3_class(fit, "gllamm_eirt")
  expect_equal(fit$n_items, n_items)
})


test_that("EIRT polytomous with multiple covariates", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(333)
  n_persons <- 150
  n_items <- 15
  n_categories <- 4

  # Multiple covariates
  item_data <- data.frame(
    difficulty_pred = rnorm(n_items, 0, 1),
    discrimination_pred = rnorm(n_items, 0, 1)
  )

  # Generate data
  theta <- rnorm(n_persons, 0, 1)
  difficulty <- 0.5 * item_data$difficulty_pred + rnorm(n_items, 0, 0.3)
  discrimination <- exp(0.3 * item_data$discrimination_pred + rnorm(n_items, 0, 0.2))

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      thresholds <- difficulty[j] + c(-1.5, -0.5, 0.5)

      probs <- numeric(n_categories)
      probs[1] <- plogis(discrimination[j] * (theta[i] - thresholds[1]))

      for (k in 2:(n_categories - 1)) {
        p_le_k <- plogis(discrimination[j] * (theta[i] - thresholds[k]))
        p_le_k_minus_1 <- plogis(discrimination[j] * (theta[i] - thresholds[k-1]))
        probs[k] <- p_le_k - p_le_k_minus_1
      }

      probs[n_categories] <- 1 - plogis(discrimination[j] * (theta[i] - thresholds[n_categories-1]))

      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Fit with multiple covariates
  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ difficulty_pred,
    discrimination_formula = ~ discrimination_pred,
    model = "GRM"
  )

  gamma_hat <- fit$regression_coefficients$difficulty
  delta_hat <- fit$regression_coefficients$discrimination

  expect_equal(length(gamma_hat), 2)
  expect_equal(length(delta_hat), 2)
})


test_that("EIRT polytomous print method", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(444)
  n_persons <- 80
  n_items <- 8
  n_categories <- 4

  item_data <- data.frame(
    x = rnorm(n_items)
  )

  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ x,
    discrimination_formula = ~ 1,
    model = "GRM"
  )

  expect_output(print(fit), "Explanatory IRT Model")
  expect_output(print(fit), "GRM")
  expect_output(print(fit), "Difficulty regression")
})


test_that("EIRT polytomous summary method", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(555)
  n_persons <- 100
  n_items <- 10
  n_categories <- 3

  item_data <- data.frame(
    covariate = rnorm(n_items)
  )

  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ covariate,
    discrimination_formula = ~ 1,
    model = "GRM"
  )

  expect_output(summary(fit), "Fitted Item Parameters")
})


test_that("EIRT polytomous with missing data", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(666)
  n_persons <- 100
  n_items <- 10
  n_categories <- 4

  item_data <- data.frame(
    x = rnorm(n_items)
  )

  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  # Introduce missing data
  missing_idx <- sample(1:length(responses), size = 0.15 * length(responses))
  responses[missing_idx] <- NA

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ x,
    discrimination_formula = ~ 1,
    model = "GRM"
  )

  expect_s3_class(fit, "gllamm_eirt")
  expect_equal(fit$n_persons, n_persons)
  expect_equal(fit$n_items, n_items)
})


test_that("EIRT polytomous residual SDs are positive", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(777)
  n_persons <- 100
  n_items <- 10
  n_categories <- 4

  item_data <- data.frame(
    x = rnorm(n_items)
  )

  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ x,
    discrimination_formula = ~ 1,
    model = "GRM"
  )

  expect_gt(fit$residual_sd$difficulty, 0)
  expect_gt(fit$residual_sd$discrimination, 0)
  expect_gt(fit$ability_sd, 0)
})


test_that("EIRT polytomous convergence check", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(888)
  n_persons <- 80
  n_items <- 8
  n_categories <- 3

  item_data <- data.frame(
    x = rnorm(n_items)
  )

  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ x,
    discrimination_formula = ~ 1,
    model = "GRM"
  )

  expect_true("convergence" %in% names(fit))
  expect_true("converged" %in% names(fit$convergence))
})


test_that("EIRT polytomous AIC/BIC computation", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(999)
  n_persons <- 100
  n_items <- 10
  n_categories <- 4

  item_data <- data.frame(
    x = rnorm(n_items)
  )

  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ x,
    discrimination_formula = ~ 1,
    model = "GRM"
  )

  expect_true("AIC" %in% names(fit))
  expect_true("BIC" %in% names(fit))
  expect_true("logLik" %in% names(fit))

  expect_gt(fit$AIC, 0)
  expect_gt(fit$BIC, 0)
})


test_that("EIRT polytomous item parameter extraction", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(1010)
  n_persons <- 100
  n_items <- 10
  n_categories <- 4

  item_data <- data.frame(
    x = rnorm(n_items)
  )

  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ x,
    discrimination_formula = ~ 1,
    model = "GRM"
  )

  # Check item parameters are returned
  expect_true("item_parameters" %in% names(fit))
  expect_true("difficulty" %in% names(fit$item_parameters))
  expect_true("discrimination" %in% names(fit$item_parameters))

  expect_equal(length(fit$item_parameters$difficulty), n_items)
  expect_equal(length(fit$item_parameters$discrimination), n_items)
})
