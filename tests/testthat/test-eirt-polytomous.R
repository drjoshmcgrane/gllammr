# Helper: simulate GRM responses using standard IRT parameterization
# P(Y >= k) = plogis(a * (theta - tau_{k-1}))
# P(Y = k) = P(Y >= k) - P(Y >= k+1)
sim_grm <- function(n_persons, n_items, theta, discrimination, thresholds_list) {
  K <- length(thresholds_list[[1]]) + 1L
  responses <- matrix(NA_integer_, n_persons, n_items)
  for (i in seq_len(n_persons)) {
    for (j in seq_len(n_items)) {
      tau <- thresholds_list[[j]]
      p_exceed <- c(1, plogis(discrimination[j] * (theta[i] - tau)), 0)
      probs <- p_exceed[-length(p_exceed)] - p_exceed[-1]
      responses[i, j] <- sample.int(K, 1, prob = pmax(probs, 0))
    }
  }
  responses
}

# Helper: simulate PCM/GPCM responses using adjacent-categories logit
# P(Y = k) proportional to exp(sum_{m=1}^k a*(theta - delta_m))
sim_pcm <- function(n_persons, n_items, theta, discrimination, thresholds_list) {
  K <- length(thresholds_list[[1]]) + 1L
  responses <- matrix(NA_integer_, n_persons, n_items)
  for (i in seq_len(n_persons)) {
    for (j in seq_len(n_items)) {
      tau <- thresholds_list[[j]]
      cumsums <- c(0, cumsum(discrimination[j] * (theta[i] - tau)))
      probs <- exp(cumsums - max(cumsums))
      probs <- probs / sum(probs)
      responses[i, j] <- sample.int(K, 1, prob = probs)
    }
  }
  responses
}


test_that("EIRT accepts polytomous data with GRM", {
  skip_if_not_installed("TMB")

  set.seed(123)
  n_persons <- 100
  n_items <- 10
  n_categories <- 4

  item_data <- data.frame(complexity = rnorm(n_items))
  theta <- rnorm(n_persons)
  discrimination <- exp(rnorm(n_items, 0, 0.3))
  thresholds_list <- replicate(n_items, c(-1.5, -0.5, 0.5), simplify = FALSE)

  responses <- sim_grm(n_persons, n_items, theta, discrimination, thresholds_list)

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


test_that("EIRT regression coefficients have correct structure for GRM", {
  skip_if_not_installed("TMB")
  skip_on_cran()  # extra EIRT fit; GRM/PCM/GPCM smoke fits run on CRAN

  set.seed(456)
  n_persons <- 120
  n_items <- 12
  n_categories <- 4

  item_data <- data.frame(item_difficulty_predictor = rnorm(n_items))
  theta <- rnorm(n_persons)
  thresholds_list <- replicate(n_items, c(-1.5, -0.5, 0.5), simplify = FALSE)
  responses <- sim_grm(n_persons, n_items, theta, rep(1.2, n_items), thresholds_list)

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ item_difficulty_predictor,
    discrimination_formula = ~ 1,
    model = "GRM"
  )

  gamma_hat <- fit$regression_coefficients$difficulty
  expect_equal(length(gamma_hat), 2)
  expect_named(gamma_hat, c("(Intercept)", "item_difficulty_predictor"))
  expect_true(is.numeric(gamma_hat))
  expect_false(any(is.na(gamma_hat)))
})


test_that("EIRT with PCM model (adjacent-categories logit)", {
  skip_if_not_installed("TMB")

  set.seed(789)
  n_persons <- 100
  n_items <- 10
  n_categories <- 4

  item_data <- data.frame(
    item_type = factor(rep(c("TypeA", "TypeB"), each = 5))
  )

  theta <- rnorm(n_persons)
  thresholds_list <- replicate(n_items, sort(rnorm(n_categories - 1)), simplify = FALSE)
  responses <- sim_pcm(n_persons, n_items, theta, rep(1, n_items), thresholds_list)

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ item_type,
    discrimination_formula = ~ 1,
    model = "PCM"
  )

  expect_s3_class(fit, "gllamm_eirt")
  expect_equal(fit$model, "PCM")

  delta_hat <- fit$regression_coefficients$discrimination
  expect_equal(length(delta_hat), 1)
})


test_that("EIRT with GPCM model (adjacent-categories + discrimination)", {
  skip_if_not_installed("TMB")

  set.seed(111)
  n_persons <- 120
  n_items <- 12
  n_categories <- 3

  item_data <- data.frame(complexity = rnorm(n_items))
  theta <- rnorm(n_persons)
  discrimination <- exp(0.3 * item_data$complexity + rnorm(n_items, 0, 0.2))
  thresholds_list <- replicate(n_items, c(-0.5, 0.5), simplify = FALSE)
  responses <- sim_pcm(n_persons, n_items, theta, discrimination, thresholds_list)

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ 1,
    discrimination_formula = ~ complexity,
    model = "GPCM"
  )

  expect_s3_class(fit, "gllamm_eirt")
  expect_equal(fit$model, "GPCM")

  delta_hat <- fit$regression_coefficients$discrimination
  expect_equal(length(delta_hat), 2)
  expect_named(delta_hat, c("(Intercept)", "complexity"))
})


test_that("EIRT PCM with threshold_formula (threshold-difficulty regression)", {
  skip_if_not_installed("TMB")

  set.seed(202)
  n_persons <- 120
  n_items <- 12
  n_categories <- 4

  item_data <- data.frame(
    abstractness = rnorm(n_items),
    cognitive_level = rnorm(n_items)
  )

  theta <- rnorm(n_persons)
  thresholds_list <- replicate(n_items, sort(rnorm(n_categories - 1)), simplify = FALSE)
  responses <- sim_pcm(n_persons, n_items, theta, rep(1, n_items), thresholds_list)

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ abstractness,
    threshold_formula = ~ cognitive_level,
    model = "PCM"
  )

  expect_s3_class(fit, "gllamm_eirt")
  expect_equal(fit$model, "PCM")

  # PCM with threshold predictors should have threshold regression coefficients
  xi_hat <- fit$regression_coefficients$threshold
  expect_false(is.null(xi_hat))
  expect_equal(nrow(xi_hat), 2)  # 2 predictors (intercept + cognitive_level)
  expect_equal(ncol(xi_hat), n_categories - 1L)  # 3 thresholds
})


test_that("EIRT with mixed number of categories", {
  skip_if_not_installed("TMB")
  skip_on_cran()  # extra EIRT fit; GRM/PCM/GPCM smoke fits run on CRAN

  set.seed(222)
  n_persons <- 100
  n_items <- 8

  item_data <- data.frame(type = factor(rep(c("Short", "Long"), each = 4)))

  responses <- matrix(NA_integer_, n_persons, n_items)
  for (i in seq_len(n_persons)) {
    responses[i, 1:4] <- sample.int(3, 4, replace = TRUE)
    responses[i, 5:8] <- sample.int(5, 4, replace = TRUE)
  }

  fit <- fit_eirt(
    responses,
    item_data = item_data,
    difficulty_formula = ~ type,
    discrimination_formula = ~ 1,
    model = "GRM"
  )

  expect_s3_class(fit, "gllamm_eirt")
  expect_equal(fit$n_items, n_items)
  expect_equal(fit$max_categories, 5L)
})


test_that("EIRT polytomous with multiple covariates", {
  skip_if_not_installed("TMB")
  skip_on_cran()  # extra EIRT fit; GRM/PCM/GPCM smoke fits run on CRAN

  set.seed(333)
  n_persons <- 120
  n_items <- 12
  n_categories <- 4

  item_data <- data.frame(
    difficulty_pred = rnorm(n_items),
    discrimination_pred = rnorm(n_items)
  )

  theta <- rnorm(n_persons)
  difficulty <- 0.5 * item_data$difficulty_pred + rnorm(n_items, 0, 0.3)
  discrimination <- exp(0.3 * item_data$discrimination_pred + rnorm(n_items, 0, 0.2))
  thresholds_list <- lapply(difficulty, function(b) b + c(-1.5, -0.5, 0.5))
  responses <- sim_grm(n_persons, n_items, theta, discrimination, thresholds_list)

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
  skip_on_cran()  # extra EIRT fit; GRM/PCM/GPCM smoke fits run on CRAN

  set.seed(444)
  n_persons <- 80; n_items <- 8; n_categories <- 4

  item_data <- data.frame(x = rnorm(n_items))
  responses <- matrix(sample.int(n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit <- fit_eirt(responses, item_data = item_data,
                  difficulty_formula = ~ x, model = "GRM")

  expect_output(print(fit), "Explanatory IRT Model")
  expect_output(print(fit), "GRM")
  expect_output(print(fit), "Difficulty regression")
})


test_that("EIRT polytomous summary method", {
  skip_if_not_installed("TMB")
  skip_on_cran()  # extra EIRT fit; GRM/PCM/GPCM smoke fits run on CRAN

  set.seed(555)
  n_persons <- 100; n_items <- 10; n_categories <- 3

  item_data <- data.frame(covariate = rnorm(n_items))
  responses <- matrix(sample.int(n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit <- fit_eirt(responses, item_data = item_data,
                  difficulty_formula = ~ covariate, model = "GRM")

  expect_output(summary(fit), "Fitted Item Parameters")
})


test_that("EIRT polytomous with missing data", {
  skip_if_not_installed("TMB")
  skip_on_cran()  # extra EIRT fit; GRM/PCM/GPCM smoke fits run on CRAN

  set.seed(666)
  n_persons <- 100; n_items <- 10; n_categories <- 4

  item_data <- data.frame(x = rnorm(n_items))
  responses <- matrix(sample.int(n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)
  responses[sample.int(length(responses), size = floor(0.15 * length(responses)))] <- NA

  fit <- fit_eirt(responses, item_data = item_data,
                  difficulty_formula = ~ x, model = "GRM")

  expect_s3_class(fit, "gllamm_eirt")
  expect_equal(fit$n_persons, n_persons)
  expect_equal(fit$n_items, n_items)
})


test_that("EIRT polytomous residual SDs are positive", {
  skip_if_not_installed("TMB")
  skip_on_cran()  # extra EIRT fit; GRM/PCM/GPCM smoke fits run on CRAN

  set.seed(777)
  n_persons <- 100; n_items <- 10; n_categories <- 4

  item_data <- data.frame(x = rnorm(n_items))
  responses <- matrix(sample.int(n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit <- fit_eirt(responses, item_data = item_data,
                  difficulty_formula = ~ x, model = "GRM")

  expect_gt(fit$residual_sd$difficulty, 0)
  expect_gt(fit$residual_sd$discrimination, 0)
  expect_gte(fit$ability_sd, 0)
})


test_that("EIRT polytomous convergence structure", {
  skip_if_not_installed("TMB")
  skip_on_cran()  # extra EIRT fit; GRM/PCM/GPCM smoke fits run on CRAN

  set.seed(888)
  n_persons <- 80; n_items <- 8; n_categories <- 3

  item_data <- data.frame(x = rnorm(n_items))
  responses <- matrix(sample.int(n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit <- fit_eirt(responses, item_data = item_data,
                  difficulty_formula = ~ x, model = "GRM")

  expect_true("convergence" %in% names(fit))
  expect_true("converged" %in% names(fit$convergence))
})


test_that("EIRT polytomous AIC/BIC computation", {
  skip_if_not_installed("TMB")
  skip_on_cran()  # extra EIRT fit; GRM/PCM/GPCM smoke fits run on CRAN

  set.seed(999)
  n_persons <- 100; n_items <- 10; n_categories <- 4

  item_data <- data.frame(x = rnorm(n_items))
  responses <- matrix(sample.int(n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit <- fit_eirt(responses, item_data = item_data,
                  difficulty_formula = ~ x, model = "GRM")

  expect_true("AIC" %in% names(fit))
  expect_true("BIC" %in% names(fit))
  expect_true("logLik" %in% names(fit))
  expect_gt(fit$AIC, 0)
  expect_gt(fit$BIC, 0)
})


test_that("EIRT polytomous item parameter extraction", {
  skip_if_not_installed("TMB")
  skip_on_cran()  # extra EIRT fit; GRM/PCM/GPCM smoke fits run on CRAN

  set.seed(1010)
  n_persons <- 100; n_items <- 10; n_categories <- 4

  item_data <- data.frame(x = rnorm(n_items))
  responses <- matrix(sample.int(n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit <- fit_eirt(responses, item_data = item_data,
                  difficulty_formula = ~ x, model = "GRM")

  expect_true("item_parameters" %in% names(fit))
  expect_true("difficulty" %in% names(fit$item_parameters))
  expect_true("discrimination" %in% names(fit$item_parameters))
  expect_equal(length(fit$item_parameters$difficulty), n_items)
  expect_equal(length(fit$item_parameters$discrimination), n_items)
})
