test_that("PCM constrains discriminations to 1", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(123)
  n_persons <- 100
  n_items <- 10
  n_categories <- 4

  # Simulate PCM data (all discriminations = 1)
  theta <- rnorm(n_persons, 0, 1)
  step_difficulties <- matrix(rnorm(n_items * (n_categories - 1)), n_items, n_categories - 1)

  # Generate responses using PCM logic
  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      # Compute cumulative sums
      cumsum_vec <- c(0, cumsum(theta[i] - step_difficulties[j, ]))

      # Compute probabilities
      exp_cumsum <- exp(cumsum_vec)
      probs <- exp_cumsum / sum(exp_cumsum)

      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Fit PCM model
  fit <- fit_irt(responses, model = "PCM")

  # Check that model type is correct
  expect_equal(fit$model, "PCM")

  # For PCM, discriminations should be approximately 1
  # (may vary slightly due to identification/scaling)
  expect_true(all(abs(fit$item_parameters$discrimination - 1) < 0.5))
})


test_that("GPCM estimates varying discriminations", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(456)
  n_persons <- 200
  n_items <- 15
  n_categories <- 5

  # Simulate GPCM data with varying discriminations
  theta <- rnorm(n_persons, 0, 1)
  discrimination <- runif(n_items, 0.5, 2.5)
  step_difficulties <- matrix(rnorm(n_items * (n_categories - 1)), n_items, n_categories - 1)

  # Generate responses
  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      cumsum_vec <- c(0, cumsum(discrimination[j] * (theta[i] - step_difficulties[j, ])))
      exp_cumsum <- exp(cumsum_vec)
      probs <- exp_cumsum / sum(exp_cumsum)
      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Fit GPCM model
  fit <- fit_irt(responses, model = "GPCM")

  # Check that discriminations vary
  expect_gt(sd(fit$item_parameters$discrimination), 0.1)

  # Check correlation with true discriminations
  expect_gt(cor(fit$item_parameters$discrimination, discrimination), 0.5)
})


test_that("PCM vs GPCM model comparison", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(789)
  n_persons <- 150
  n_items <- 12
  n_categories <- 4

  # Simulate data that could be either model
  theta <- rnorm(n_persons, 0, 1)
  discrimination <- rep(1, n_items)  # Equal discriminations
  step_difficulties <- matrix(rnorm(n_items * (n_categories - 1)), n_items, n_categories - 1)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      cumsum_vec <- c(0, cumsum(discrimination[j] * (theta[i] - step_difficulties[j, ])))
      exp_cumsum <- exp(cumsum_vec)
      probs <- exp_cumsum / sum(exp_cumsum)
      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Fit both models
  fit_pcm <- fit_irt(responses, model = "PCM")
  fit_gpcm <- fit_irt(responses, model = "GPCM")

  # Both should converge
  expect_true(fit_pcm$convergence$converged)
  expect_true(fit_gpcm$convergence$converged)

  # PCM should have lower AIC (fewer parameters) when data is truly PCM
  # GPCM has more parameters but similar fit
  expect_lt(fit_pcm$AIC, fit_gpcm$AIC + 5)  # Within 5 AIC units
})


test_that("PCM parameter recovery", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(111)
  n_persons <- 500
  n_items <- 20
  n_categories <- 3

  # True parameters
  true_theta <- rnorm(n_persons, 0, 1)
  true_steps <- matrix(c(-1, 1), n_items, 2, byrow = TRUE)

  # Generate data
  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      cumsum_vec <- c(0, cumsum(true_theta[i] - true_steps[j, ]))
      exp_cumsum <- exp(cumsum_vec)
      probs <- exp_cumsum / sum(exp_cumsum)
      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Fit model
  fit <- fit_irt(responses, model = "PCM")

  # Check ability correlation
  expect_gt(cor(fit$person_abilities, true_theta), 0.8)

  # Check step difficulty recovery (within reasonable tolerance)
  for (j in 1:min(5, n_items)) {
    recovered_steps <- fit$item_parameters$thresholds[[j]]
    true_steps_j <- true_steps[j, ]
    # Allow for identification/scaling differences
    expect_equal(cor(recovered_steps, true_steps_j), 1, tolerance = 0.3)
  }
})


test_that("GPCM with 7 categories", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(222)
  n_persons <- 100
  n_items <- 8
  n_categories <- 7

  # Large number of categories
  theta <- rnorm(n_persons, 0, 1)
  discrimination <- runif(n_items, 0.8, 1.5)
  step_difficulties <- matrix(rnorm(n_items * (n_categories - 1), sd = 0.5),
                              n_items, n_categories - 1)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      cumsum_vec <- c(0, cumsum(discrimination[j] * (theta[i] - step_difficulties[j, ])))
      exp_cumsum <- exp(cumsum_vec)
      probs <- exp_cumsum / sum(exp_cumsum)
      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Fit GPCM
  fit <- fit_irt(responses, model = "GPCM")

  # Should handle 7 categories
  expect_equal(fit$max_categories, 7)
  expect_true(all(fit$n_categories == 7))

  # Each item should have 6 thresholds
  for (j in 1:n_items) {
    expect_equal(length(fit$item_parameters$thresholds[[j]]), 6)
  }
})


test_that("PCM with missing data patterns", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(333)
  n_persons <- 120
  n_items <- 10
  n_categories <- 4

  # Generate complete data
  theta <- rnorm(n_persons, 0, 1)
  step_difficulties <- matrix(rnorm(n_items * (n_categories - 1)), n_items, n_categories - 1)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      cumsum_vec <- c(0, cumsum(theta[i] - step_difficulties[j, ]))
      exp_cumsum <- exp(cumsum_vec)
      probs <- exp_cumsum / sum(exp_cumsum)
      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Introduce structured missingness (some items missing for some persons)
  responses[1:20, 1:3] <- NA  # First 20 persons missing first 3 items
  responses[50:70, 7:10] <- NA  # Middle group missing last 4 items

  # Fit model
  fit <- fit_irt(responses, model = "PCM")

  # Should complete successfully
  expect_s3_class(fit, "gllamm_irt")
  expect_true(fit$convergence$converged)
})


test_that("GPCM print and summary methods", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(444)
  n_persons <- 80
  n_items <- 8
  n_categories <- 4

  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit <- fit_irt(responses, model = "GPCM")

  # Print should show polytomous model info
  expect_output(print(fit), "GPCM")
  expect_output(print(fit), "Polytomous")
  expect_output(print(fit), "Item discriminations")

  # Summary should show ability quartiles
  expect_output(summary(fit), "Ability quartiles")
})


test_that("PCM vs GPCM discrimination constraints", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(555)
  n_persons <- 100
  n_items <- 10
  n_categories <- 3

  # Same data for both models
  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  fit_pcm <- fit_irt(responses, model = "PCM")
  fit_gpcm <- fit_irt(responses, model = "GPCM")

  # PCM discriminations should be close to 1
  expect_true(all(abs(fit_pcm$item_parameters$discrimination - 1) < 0.3))

  # GPCM discriminations should be free to vary
  # (though with random data, may not vary much)
  expect_true(is.numeric(fit_gpcm$item_parameters$discrimination))
  expect_equal(length(fit_gpcm$item_parameters$discrimination), n_items)
})


test_that("PCM/GPCM handle extreme responses", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(666)
  n_persons <- 100
  n_items <- 10
  n_categories <- 5

  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  # Force some extreme responses (all 1s or all 5s for some persons)
  responses[1:5, ] <- 1   # Low ability persons
  responses[96:100, ] <- 5  # High ability persons

  # Both models should handle this
  fit_pcm <- fit_irt(responses, model = "PCM")
  fit_gpcm <- fit_irt(responses, model = "GPCM")

  expect_s3_class(fit_pcm, "gllamm_irt")
  expect_s3_class(fit_gpcm, "gllamm_irt")

  # Check that extreme persons have extreme ability estimates
  expect_lt(mean(fit_pcm$person_abilities[1:5]), -1)
  expect_gt(mean(fit_pcm$person_abilities[96:100]), 1)
})
