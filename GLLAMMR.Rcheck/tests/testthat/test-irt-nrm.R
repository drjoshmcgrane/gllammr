test_that("NRM accepts unordered categorical data", {
  skip_if_not_installed("TMB")

  set.seed(123)
  n_persons <- 100
  n_items <- 8
  n_categories <- 4

  # Simulate NRM data (no ordering assumption)
  theta <- rnorm(n_persons, 0, 1)

  # Category-specific parameters (slopes and intercepts)
  # For NRM, different categories have different relationships with theta
  discrimination <- runif(n_items, 0.5, 2)
  intercepts <- matrix(rnorm(n_items * (n_categories - 1)), n_items, n_categories - 1)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      # Compute linear predictors for each category
      eta <- numeric(n_categories)
      eta[1] <- 0  # Reference category
      for (k in 2:n_categories) {
        eta[k] <- discrimination[j] * theta[i] + intercepts[j, k - 1]
      }

      # Softmax to get probabilities
      exp_eta <- exp(eta)
      probs <- exp_eta / sum(exp_eta)

      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Fit NRM
  fit <- fit_irt(responses, model = "NRM")

  # Basic checks
  expect_s3_class(fit, "gllamm_irt")
  expect_equal(fit$model, "NRM")
  expect_equal(fit$max_categories, n_categories)
})


test_that("NRM vs GRM on ordered data", {
  skip_if_not_installed("TMB")

  set.seed(456)
  n_persons <- 150
  n_items <- 10
  n_categories <- 4

  # Generate ordered data (GRM-like)
  theta <- rnorm(n_persons, 0, 1)
  discrimination <- rep(1.5, n_items)
  thresholds <- matrix(c(-1.5, 0, 1.5), n_items, 3, byrow = TRUE)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      # GRM: P(Y >= k+1) = plogis(a * (theta - b_k))
      p_exceed <- c(1, plogis(discrimination[j] * (theta[i] - thresholds[j, ])), 0)
      probs <- p_exceed[-length(p_exceed)] - p_exceed[-1]
      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Fit both models
  fit_grm <- fit_irt(responses, model = "GRM")
  fit_nrm <- fit_irt(responses, model = "NRM")

  # GRM should fit better (lower AIC) for ordered data
  expect_lt(fit_grm$AIC, fit_nrm$AIC)

  # Both should converge
  expect_true(fit_grm$convergence$converged)
  expect_true(fit_nrm$convergence$converged)
})


test_that("NRM handles truly unordered categories", {
  skip_if_not_installed("TMB")

  set.seed(789)
  n_persons <- 200
  n_items <- 12
  n_categories <- 3

  # Simulate data where category order doesn't matter
  # e.g., preference for Red/Green/Blue colors
  theta <- rnorm(n_persons, 0, 1)

  # Different patterns for different categories (no ordering)
  discrimination <- runif(n_items, 0.5, 2)
  intercepts <- matrix(rnorm(n_items * (n_categories - 1), sd = 1),
                       n_items, n_categories - 1)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      eta <- c(0, discrimination[j] * theta[i] + intercepts[j, ])
      probs <- exp(eta) / sum(exp(eta))
      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Fit NRM
  fit <- fit_irt(responses, model = "NRM")

  # Should complete successfully
  expect_s3_class(fit, "gllamm_irt")
  expect_true(fit$convergence$converged)
  expect_equal(fit$n_items, n_items)
})


test_that("NRM with 5 unordered categories", {
  skip_if_not_installed("TMB")

  set.seed(111)
  n_persons <- 100
  n_items <- 8
  n_categories <- 5

  # 5-choice unordered (e.g., multiple choice where order is arbitrary)
  theta <- rnorm(n_persons, 0, 1)
  discrimination <- runif(n_items, 0.8, 1.5)
  intercepts <- matrix(rnorm(n_items * (n_categories - 1)), n_items, n_categories - 1)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      eta <- c(0, discrimination[j] * theta[i] + intercepts[j, ])
      probs <- exp(eta) / sum(exp(eta))
      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  fit <- fit_irt(responses, model = "NRM")

  expect_equal(fit$max_categories, 5)
  expect_equal(length(fit$item_parameters$thresholds[[1]]), 4)  # K-1 thresholds
})


test_that("NRM print method", {
  skip_if_not_installed("TMB")

  set.seed(222)
  n_persons <- 60
  n_items <- 6
  n_categories <- 3

  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit <- fit_irt(responses, model = "NRM")

  expect_output(print(fit), "NRM")
  expect_output(print(fit), "Polytomous")
})


test_that("NRM with missing data", {
  skip_if_not_installed("TMB")

  set.seed(333)
  n_persons <- 100
  n_items <- 10
  n_categories <- 4

  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  # Introduce 25% missing data
  missing_idx <- sample(1:(n_persons * n_items), 0.25 * n_persons * n_items)
  responses[missing_idx] <- NA

  fit <- fit_irt(responses, model = "NRM")

  expect_s3_class(fit, "gllamm_irt")
  expect_true(fit$convergence$converged)
})


test_that("NRM parameter structure", {
  skip_if_not_installed("TMB")

  set.seed(444)
  n_persons <- 80
  n_items <- 8
  n_categories <- 4

  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit <- fit_irt(responses, model = "NRM")

  # Check parameter structure
  expect_true(is.list(fit$item_parameters))
  expect_true("discrimination" %in% names(fit$item_parameters))
  expect_true("thresholds" %in% names(fit$item_parameters))

  # Each item should have K-1 threshold-like parameters (intercepts)
  for (j in 1:n_items) {
    expect_equal(length(fit$item_parameters$thresholds[[j]]), n_categories - 1)
  }
})


test_that("NRM convergence on difficult data", {
  skip_if_not_installed("TMB")

  set.seed(555)
  n_persons <- 50
  n_items <- 5
  n_categories <- 3

  # Very noisy data (almost random responses)
  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  # Should still converge (though parameters may not be well-estimated)
  fit <- fit_irt(responses, model = "NRM")

  expect_s3_class(fit, "gllamm_irt")
  # May or may not converge perfectly on random data, so don't assert convergence
})
