test_that("EIRT accepts valid input for Rasch model", {
  skip_if_not_installed("TMB")

  set.seed(123)
  n_persons <- 100
  n_items <- 10

  # Create item-level covariates
  item_data <- data.frame(
    word_freq = rnorm(n_items, 0, 1),
    item_length = rpois(n_items, 5)
  )

  # Simulate responses
  theta <- rnorm(n_persons, 0, 1)
  difficulty <- rnorm(n_items, 0, 1)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      p <- plogis(theta[i] - difficulty[j])
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  # Fit EIRT model
  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ word_freq + item_length,
    discrimination_formula = ~ 1,
    model = "Rasch"
  )

  expect_s3_class(fit, "gllamm_eirt")
  expect_s3_class(fit, "gllamm")
  expect_equal(fit$model, "Rasch")
  expect_equal(fit$n_persons, n_persons)
  expect_equal(fit$n_items, n_items)
})


test_that("EIRT recovers known item covariate effects", {
  skip_if_not_installed("TMB")

  set.seed(456)
  n_persons <- 200
  n_items <- 20

  # Create item covariates
  item_data <- data.frame(
    word_freq = rnorm(n_items, 0, 1)
  )

  # True parameters
  gamma_0 <- 0.5   # Intercept
  gamma_1 <- -0.8  # Word frequency effect (higher freq = easier)
  sigma_epsilon_b <- 0.3

  theta <- rnorm(n_persons, 0, 1)

  # Generate difficulty from covariates
  difficulty <- gamma_0 + gamma_1 * item_data$word_freq + rnorm(n_items, 0, sigma_epsilon_b)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      p <- plogis(theta[i] - difficulty[j])
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  # Fit EIRT model
  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ word_freq,
    discrimination_formula = ~ 1,
    model = "Rasch"
  )

  # Check recovery of regression coefficients
  gamma_hat <- fit$regression_coefficients$difficulty

  expect_equal(length(gamma_hat), 2)
  expect_named(gamma_hat, c("(Intercept)", "word_freq"))

  # Should recover within reasonable tolerance
  expect_equal(gamma_hat[["(Intercept)"]], gamma_0, tolerance = 0.3)
  expect_equal(gamma_hat[["word_freq"]], gamma_1, tolerance = 0.3)
})


test_that("EIRT with 2PL model and discrimination covariates", {
  skip_if_not_installed("TMB")

  set.seed(789)
  n_persons <- 150
  n_items <- 15

  # Create item covariates
  item_data <- data.frame(
    item_type = factor(rep(c("easy", "hard"), length.out = n_items)),
    complexity = rnorm(n_items, 0, 1)
  )

  # True parameters
  gamma_0 <- 0.0
  delta_0 <- 0.5  # Log-discrimination intercept
  delta_1 <- 0.3  # Complexity effect on discrimination

  theta <- rnorm(n_persons, 0, 1)

  difficulty <- gamma_0 + rnorm(n_items, 0, 0.5)
  log_discrimination <- delta_0 + delta_1 * item_data$complexity + rnorm(n_items, 0, 0.2)
  discrimination <- exp(log_discrimination)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      p <- plogis(discrimination[j] * (theta[i] - difficulty[j]))
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  # Fit EIRT model with discrimination covariates
  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ 1,
    discrimination_formula = ~ complexity,
    model = "2PL"
  )

  expect_s3_class(fit, "gllamm_eirt")
  expect_equal(fit$model, "2PL")

  # Check that discrimination coefficients are estimated
  delta_hat <- fit$regression_coefficients$discrimination
  expect_equal(length(delta_hat), 2)
  expect_named(delta_hat, c("(Intercept)", "complexity"))
})


test_that("EIRT with multiple covariates", {
  skip_if_not_installed("TMB")

  set.seed(111)
  n_persons <- 120
  n_items <- 12

  # Multiple item covariates
  item_data <- data.frame(
    word_freq = rnorm(n_items, 0, 1),
    word_length = rpois(n_items, 6),
    is_abstract = rbinom(n_items, 1, 0.5)
  )

  # Generate responses
  theta <- rnorm(n_persons, 0, 1)
  difficulty <- 0.5 - 0.6 * item_data$word_freq +
                0.3 * item_data$word_length +
                0.8 * item_data$is_abstract +
                rnorm(n_items, 0, 0.4)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      p <- plogis(theta[i] - difficulty[j])
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  # Fit with multiple covariates
  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ word_freq + word_length + is_abstract,
    discrimination_formula = ~ 1,
    model = "Rasch"
  )

  gamma_hat <- fit$regression_coefficients$difficulty

  expect_equal(length(gamma_hat), 4)
  expect_named(gamma_hat, c("(Intercept)", "word_freq", "word_length", "is_abstract"))

  # Check signs of effects (word_freq should be negative, length/abstract positive)
  expect_lt(gamma_hat[["word_freq"]], 0)
  expect_gt(gamma_hat[["word_length"]], 0)
})


test_that("EIRT print method works", {
  skip_if_not_installed("TMB")
  skip_on_cran()  # extra EIRT fit; Rasch/2PL smoke fits run on CRAN

  set.seed(222)
  n_persons <- 100
  n_items <- 10

  item_data <- data.frame(
    covariate1 = rnorm(n_items)
  )

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ covariate1,
    discrimination_formula = ~ 1,
    model = "Rasch"
  )

  # Print should not error
  expect_output(print(fit), "Explanatory IRT Model")
  expect_output(print(fit), "Rasch")
  expect_output(print(fit), "Difficulty regression")
  expect_output(print(fit), "covariate1")
})


test_that("EIRT summary method works", {
  skip_if_not_installed("TMB")
  skip_on_cran()  # extra EIRT fit; Rasch/2PL smoke fits run on CRAN

  set.seed(333)
  n_persons <- 80
  n_items <- 8

  item_data <- data.frame(
    x1 = rnorm(n_items)
  )

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ x1,
    discrimination_formula = ~ 1,
    model = "Rasch"
  )

  expect_output(summary(fit), "Fitted Item Parameters")
  expect_output(summary(fit), "first 10 items")
})


test_that("EIRT validates item_data dimensions", {
  skip_if_not_installed("TMB")

  set.seed(444)
  n_persons <- 50
  n_items <- 10

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)

  # Wrong number of rows in item_data
  item_data_wrong <- data.frame(
    covariate = rnorm(5)  # Only 5 rows, should be 10
  )

  expect_error(
    fit_eirt(
      responses,
      item_data = item_data_wrong,
      difficulty_formula = ~ covariate,
      model = "Rasch"
    ),
    "item_data must have 10 rows"
  )
})


test_that("EIRT with intercept-only formulas", {
  skip_if_not_installed("TMB")

  set.seed(555)
  n_persons <- 100
  n_items <- 10

  item_data <- data.frame(
    dummy = rep(1, n_items)  # Just a placeholder
  )

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)

  # Intercept-only for both difficulty and discrimination
  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ 1,
    discrimination_formula = ~ 1,
    model = "2PL"
  )

  expect_s3_class(fit, "gllamm_eirt")

  # Should have only intercepts
  gamma_hat <- fit$regression_coefficients$difficulty
  delta_hat <- fit$regression_coefficients$discrimination

  expect_equal(length(gamma_hat), 1)
  expect_equal(length(delta_hat), 1)
  expect_named(gamma_hat, "(Intercept)")
  expect_named(delta_hat, "(Intercept)")
})


test_that("EIRT handles missing data", {
  skip_if_not_installed("TMB")

  set.seed(666)
  n_persons <- 100
  n_items <- 10

  item_data <- data.frame(
    covariate = rnorm(n_items)
  )

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)

  # Introduce 20% missing data
  missing_idx <- sample(1:length(responses), size = 0.2 * length(responses))
  responses[missing_idx] <- NA

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ covariate,
    discrimination_formula = ~ 1,
    model = "Rasch"
  )

  expect_s3_class(fit, "gllamm_eirt")
  expect_equal(fit$n_persons, n_persons)
  expect_equal(fit$n_items, n_items)
})


test_that("EIRT residual standard deviations are positive", {
  skip_if_not_installed("TMB")
  skip_on_cran()  # extra EIRT fit; Rasch/2PL smoke fits run on CRAN

  set.seed(777)
  n_persons <- 100
  n_items <- 10

  item_data <- data.frame(
    x = rnorm(n_items)
  )

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ x,
    discrimination_formula = ~ 1,
    model = "2PL"
  )

  # Check residual SDs are positive
  expect_gt(fit$residual_sd$difficulty, 0)
  expect_gt(fit$residual_sd$discrimination, 0)
  expect_gt(fit$ability_sd, 0)
})


test_that("EIRT comparison with two-stage approach", {
  skip_if_not_installed("TMB")
  skip_on_cran()  # multiple EIRT/IRT fits for the two-stage comparison; CI-only

  set.seed(888)
  n_persons <- 150
  n_items <- 15

  # Create strong covariate effect
  item_data <- data.frame(
    word_freq = rnorm(n_items, 0, 1)
  )

  gamma_0 <- 0.0
  gamma_1 <- -1.0  # Strong negative effect

  theta <- rnorm(n_persons, 0, 1)
  difficulty <- gamma_0 + gamma_1 * item_data$word_freq + rnorm(n_items, 0, 0.2)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      p <- plogis(theta[i] - difficulty[j])
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  # Direct EIRT fit
  fit_eirt <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ word_freq,
    discrimination_formula = ~ 1,
    model = "Rasch"
  )

  # Two-stage: fit_irt() then regress difficulties
  fit_irt <- fit_irt(responses, model = "Rasch", se = FALSE)
  difficulties <- fit_irt$item_parameters$difficulty

  two_stage_lm <- lm(difficulties ~ word_freq, data = item_data)

  # Compare coefficients (should be similar)
  gamma_direct <- fit_eirt$regression_coefficients$difficulty
  gamma_two_stage <- coef(two_stage_lm)

  # Correlation should be high
  expect_gt(cor(gamma_direct, gamma_two_stage), 0.95)
})


test_that("EIRT with categorical item covariates", {
  skip_if_not_installed("TMB")

  set.seed(999)
  n_persons <- 120
  n_items <- 12

  # Categorical covariate
  item_data <- data.frame(
    item_type = factor(rep(c("Type_A", "Type_B", "Type_C"), each = 4))
  )

  theta <- rnorm(n_persons, 0, 1)

  # Different difficulties for each type
  difficulty <- numeric(n_items)
  difficulty[item_data$item_type == "Type_A"] <- rnorm(4, -0.5, 0.3)
  difficulty[item_data$item_type == "Type_B"] <- rnorm(4, 0.0, 0.3)
  difficulty[item_data$item_type == "Type_C"] <- rnorm(4, 0.5, 0.3)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      p <- plogis(theta[i] - difficulty[j])
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ item_type,
    discrimination_formula = ~ 1,
    model = "Rasch"
  )

  gamma_hat <- fit$regression_coefficients$difficulty

  # Should have 3 coefficients (intercept + 2 dummy variables)
  expect_equal(length(gamma_hat), 3)
  expect_true("(Intercept)" %in% names(gamma_hat))
})


test_that("EIRT convergence information is returned", {
  skip_if_not_installed("TMB")
  skip_on_cran()  # extra EIRT fit; Rasch/2PL smoke fits run on CRAN

  set.seed(1111)
  n_persons <- 80
  n_items <- 8

  item_data <- data.frame(
    x = rnorm(n_items)
  )

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ x,
    discrimination_formula = ~ 1,
    model = "Rasch"
  )

  expect_true("convergence" %in% names(fit))
  expect_true("converged" %in% names(fit$convergence))
  expect_true("message" %in% names(fit$convergence))
})


test_that("EIRT AIC and BIC are computed", {
  skip_if_not_installed("TMB")
  skip_on_cran()  # extra EIRT fit; Rasch/2PL smoke fits run on CRAN

  set.seed(1212)
  n_persons <- 100
  n_items <- 10

  item_data <- data.frame(
    x = rnorm(n_items)
  )

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ x,
    discrimination_formula = ~ 1,
    model = "Rasch"
  )

  expect_true("AIC" %in% names(fit))
  expect_true("BIC" %in% names(fit))
  expect_true("logLik" %in% names(fit))

  # AIC and BIC should be positive (larger is worse)
  expect_gt(fit$AIC, 0)
  expect_gt(fit$BIC, 0)
})
