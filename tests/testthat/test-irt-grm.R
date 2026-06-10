test_that("GRM accepts valid polytomous input", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(123)
  n_persons <- 100
  n_items <- 10
  n_categories <- 5

  # Simulate GRM data
  theta <- rnorm(n_persons, 0, 1)
  discrimination <- runif(n_items, 0.5, 2)

  # Generate threshold parameters (ordered)
  thresholds <- matrix(NA, n_items, n_categories - 1)
  for (j in 1:n_items) {
    thresholds[j, ] <- sort(rnorm(n_categories - 1, 0, 1))
  }

  # Generate responses
  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      # Compute probabilities for each category
      probs <- numeric(n_categories)

      # P(Y = 1)
      probs[1] <- plogis(discrimination[j] * (theta[i] - thresholds[j, 1]))

      # P(Y = k) for k = 2, ..., K-1
      for (k in 2:(n_categories-1)) {
        probs[k] <- plogis(discrimination[j] * (theta[i] - thresholds[j, k])) -
                    plogis(discrimination[j] * (theta[i] - thresholds[j, k-1]))
      }

      # P(Y = K)
      probs[n_categories] <- 1 - plogis(discrimination[j] * (theta[i] - thresholds[j, n_categories-1]))

      # Sample response
      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Fit GRM model
  fit <- fit_irt(responses, model = "GRM")

  # Basic checks
  expect_s3_class(fit, "gllamm_irt")
  expect_equal(fit$model, "GRM")
  expect_equal(fit$n_persons, n_persons)
  expect_equal(fit$n_items, n_items)
  expect_equal(fit$max_categories, n_categories)
})


test_that("GRM parameter recovery with known parameters", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(456)
  n_persons <- 500
  n_items <- 15
  n_categories <- 4

  # True parameters
  true_theta <- rnorm(n_persons, 0, 1)
  true_discrimination <- rep(1.5, n_items)  # Constant discrimination
  true_thresholds <- matrix(c(-1.5, 0, 1.5), n_items, 3, byrow = TRUE)

  # Generate responses
  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      probs <- numeric(n_categories)
      probs[1] <- plogis(true_discrimination[j] * (true_theta[i] - true_thresholds[j, 1]))
      for (k in 2:(n_categories-1)) {
        probs[k] <- plogis(true_discrimination[j] * (true_theta[i] - true_thresholds[j, k])) -
                    plogis(true_discrimination[j] * (true_theta[i] - true_thresholds[j, k-1]))
      }
      probs[n_categories] <- 1 - plogis(true_discrimination[j] * (true_theta[i] - true_thresholds[j, n_categories-1]))
      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Fit model
  fit <- fit_irt(responses, model = "GRM")

  # Check parameter recovery (with reasonable tolerance)
  # Discrimination
  expect_equal(mean(fit$item_parameters$discrimination),
               mean(true_discrimination),
               tolerance = 0.3)

  # Threshold correlation (accounting for identification/scaling)
  for (j in 1:min(5, n_items)) {
    recovered_thresh <- fit$item_parameters$thresholds[[j]]
    true_thresh <- true_thresholds[j, ]
    expect_equal(cor(recovered_thresh, true_thresh), 1, tolerance = 0.2)
  }

  # Ability correlation
  expect_gt(cor(fit$person_abilities, true_theta), 0.7)
})


test_that("GRM handles different numbers of categories", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(789)
  n_persons <- 80
  n_items <- 8

  # Items with different numbers of categories
  responses <- matrix(NA, n_persons, n_items)

  # Items 1-4: 3 categories
  for (i in 1:n_persons) {
    for (j in 1:4) {
      responses[i, j] <- sample(1:3, 1)
    }
  }

  # Items 5-8: 5 categories
  for (i in 1:n_persons) {
    for (j in 5:8) {
      responses[i, j] <- sample(1:5, 1)
    }
  }

  # Fit model
  fit <- fit_irt(responses, model = "GRM")

  # Check that different categories are recognized
  expect_equal(fit$n_categories[1:4], rep(3, 4))
  expect_equal(fit$n_categories[5:8], rep(5, 4))
  expect_equal(fit$max_categories, 5)
})


test_that("GRM handles missing data", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(111)
  n_persons <- 100
  n_items <- 10
  n_categories <- 4

  # Generate complete data
  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  # Introduce 20% missing data
  missing_idx <- sample(1:(n_persons * n_items), 0.2 * n_persons * n_items)
  responses[missing_idx] <- NA

  # Fit model
  fit <- fit_irt(responses, model = "GRM")

  # Should complete without error
  expect_s3_class(fit, "gllamm_irt")
  expect_true(fit$convergence$converged)
})


test_that("GRM print method works for polytomous models", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(222)
  n_persons <- 50
  n_items <- 5
  n_categories <- 4

  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit <- fit_irt(responses, model = "GRM")

  # Print should work without error
  expect_output(print(fit), "GRM")
  expect_output(print(fit), "Polytomous")
  expect_output(print(fit), "Item discriminations")
})


test_that("GRM validates response coding", {
  skip_if_not_installed("TMB")

  set.seed(333)
  n_persons <- 50
  n_items <- 5

  # 0-based coding is auto-recoded to 1-based with a message
  responses_zero_based <- matrix(sample(0:4, n_persons * n_items, replace = TRUE),
                                 n_persons, n_items)
  expect_message(fit_irt(responses_zero_based, model = "GRM"),
                 "Auto-recoding")

  # Truly invalid coding (categories starting at 2) must error
  responses_invalid <- matrix(sample(2:6, n_persons * n_items, replace = TRUE),
                              n_persons, n_items)
  expect_error(fit_irt(responses_invalid, model = "GRM"),
               "invalid response coding")
})


test_that("GRM model type dispatch works correctly", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(444)
  n_persons <- 50
  n_items <- 5

  # Binary data
  responses_binary <- matrix(sample(0:1, n_persons * n_items, replace = TRUE),
                              n_persons, n_items)

  # Should error when trying to fit GRM to binary data
  expect_warning(fit_irt(responses_binary, model = "GRM"),
                 "data appears dichotomous")

  # Polytomous data
  responses_poly <- matrix(sample(1:5, n_persons * n_items, replace = TRUE),
                           n_persons, n_items)

  # Should error when trying to fit Rasch to polytomous data
  expect_error(fit_irt(responses_poly, model = "Rasch"),
               "Dichotomous models.*require binary responses")
})
