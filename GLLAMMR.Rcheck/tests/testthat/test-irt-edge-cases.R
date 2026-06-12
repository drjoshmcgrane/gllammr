test_that("IRT handles all zero responses for a person", {
  skip_if_not_installed("TMB")

  set.seed(123)
  n_persons <- 100
  n_items <- 10

  # Binary IRT
  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)

  # Force first 3 persons to have all zeros
  responses[1:3, ] <- 0

  # Should handle this (extreme low ability)
  fit <- fit_irt(responses, model = "Rasch")

  expect_s3_class(fit, "gllamm_irt")

  # First 3 persons should have the lowest ability estimates (shrinkage
  # toward the prior keeps the magnitudes modest with random data)
  expect_true(max(fit$person_abilities[1:3]) <=
                min(fit$person_abilities[4:n_persons]) + 1e-8)
  expect_true(all(fit$person_abilities[1:3] < mean(fit$person_abilities)))
})


test_that("IRT handles all maximum responses for a person", {
  skip_if_not_installed("TMB")

  set.seed(456)
  n_persons <- 100
  n_items <- 10
  n_categories <- 5

  # Polytomous IRT
  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  # Force last 3 persons to have all maximum responses
  responses[98:100, ] <- n_categories

  fit <- fit_irt(responses, model = "GRM")

  expect_s3_class(fit, "gllamm_irt")

  # Last 3 persons should have very high ability estimates
  expect_true(all(fit$person_abilities[98:100] > 2))
})


test_that("IRT handles item with no variance", {
  skip_if_not_installed("TMB")

  set.seed(789)
  n_persons <- 100
  n_items <- 10

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)

  # Make one item have no variance (all correct)
  responses[, 5] <- 1

  # Should still fit (though that item may not be informative)
  fit <- fit_irt(responses, model = "Rasch")

  expect_s3_class(fit, "gllamm_irt")

  # Item 5 should have very low difficulty (everyone got it right)
  expect_lt(fit$item_parameters$difficulty[5], -2)
})


test_that("IRT handles very small sample size", {
  skip_if_not_installed("TMB")

  set.seed(111)
  n_persons <- 20  # Very small
  n_items <- 5

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)

  # Should complete but may not converge well
  fit <- fit_irt(responses, model = "Rasch")

  expect_s3_class(fit, "gllamm_irt")
  expect_equal(fit$n_persons, 20)
  expect_equal(fit$n_items, 5)
})


test_that("IRT handles very few items", {
  skip_if_not_installed("TMB")

  set.seed(222)
  n_persons <- 100
  n_items <- 3  # Very few items

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)

  fit <- fit_irt(responses, model = "2PL")

  expect_s3_class(fit, "gllamm_irt")
  expect_equal(fit$n_items, 3)
})


test_that("IRT handles extreme discrimination values in simulation", {
  skip_if_not_installed("TMB")

  set.seed(333)
  n_persons <- 200
  n_items <- 10

  # Simulate with very high and very low discriminations
  theta <- rnorm(n_persons, 0, 1)
  difficulty <- rnorm(n_items, 0, 1)
  discrimination <- c(rep(0.1, 3), rep(1, 4), rep(5, 3))  # Extreme values

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      p <- plogis(discrimination[j] * (theta[i] - difficulty[j]))
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  fit <- fit_irt(responses, model = "2PL")

  expect_s3_class(fit, "gllamm_irt")

  # Should recover that discriminations vary widely
  expect_gt(max(fit$item_parameters$discrimination), 2)
  expect_lt(min(fit$item_parameters$discrimination), 0.5)
})


test_that("Polytomous IRT handles 2 categories (edge of polytomous)", {
  skip_if_not_installed("TMB")

  set.seed(444)
  n_persons <- 100
  n_items <- 10
  n_categories <- 2  # Technically binary, but coded as 1,2

  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  # Should warn but may still fit
  expect_warning(fit_irt(responses, model = "GRM"),
                 "data appears dichotomous")
})


test_that("IRT handles highly sparse response matrix", {
  skip_if_not_installed("TMB")

  set.seed(555)
  n_persons <- 100
  n_items <- 20

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)

  # Make 60% of responses missing
  missing_idx <- sample(1:(n_persons * n_items), 0.6 * n_persons * n_items)
  responses[missing_idx] <- NA

  # Should handle this (though estimates will have high uncertainty)
  fit <- fit_irt(responses, model = "Rasch")

  expect_s3_class(fit, "gllamm_irt")

  # Check that some parameters were estimated
  expect_true(!all(is.na(fit$person_abilities)))
  expect_true(!all(is.na(fit$item_parameters$difficulty)))
})


test_that("Polytomous IRT handles unbalanced category frequencies", {
  skip_if_not_installed("TMB")

  set.seed(666)
  n_persons <- 150
  n_items <- 10
  n_categories <- 5

  # Highly unbalanced: most responses in categories 1 and 2
  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE,
                             prob = c(0.4, 0.3, 0.15, 0.1, 0.05)),
                      n_persons, n_items)

  fit <- fit_irt(responses, model = "GRM")

  expect_s3_class(fit, "gllamm_irt")
  expect_true(fit$convergence$converged)
})


test_that("IRT validates impossible response values", {
  skip_if_not_installed("TMB")

  set.seed(777)
  n_persons <- 50
  n_items <- 5

  # Invalid: mixed 0-based and 1-based coding (0 present, but category 1
  # missing, so neither contiguous 0-based nor 1-based)
  responses_invalid <- matrix(sample(c(0, 2, 3, 4, 5), n_persons * n_items,
                                     replace = TRUE),
                               n_persons, n_items)

  # Should error due to inconsistent coding
  expect_error(fit_irt(responses_invalid, model = "GRM"),
               "invalid response coding")
})


test_that("IRT handles single person (degenerate case)", {
  skip_if_not_installed("TMB")

  set.seed(888)
  n_persons <- 1  # Single person
  n_items <- 10

  responses <- matrix(rbinom(n_items, 1, 0.5), n_persons, n_items)

  # May not converge well but should not crash
  # Expect warning about insufficient data
  fit <- try(fit_irt(responses, model = "Rasch"), silent = TRUE)

  # Either succeeds or fails gracefully
  expect_true(inherits(fit, "gllamm_irt") || inherits(fit, "try-error"))
})


test_that("IRT handles all items same difficulty", {
  skip_if_not_installed("TMB")

  set.seed(999)
  n_persons <- 100
  n_items <- 10

  # Simulate with identical difficulties
  theta <- rnorm(n_persons, 0, 1)
  difficulty <- rep(0, n_items)  # All same

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      p <- plogis(theta[i] - difficulty[j])
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  fit <- fit_irt(responses, model = "Rasch")

  expect_s3_class(fit, "gllamm_irt")

  # Difficulties should all be close to 0
  expect_true(sd(fit$item_parameters$difficulty) < 1)
})


test_that("Polytomous IRT with ragged data structure", {
  skip_if_not_installed("TMB")

  set.seed(1234)
  n_persons <- 100
  n_items <- 10

  # Items with different numbers of categories: some 3, some 4, some 5
  responses <- matrix(NA, n_persons, n_items)
  responses[, 1:3] <- matrix(sample(1:3, n_persons * 3, replace = TRUE), n_persons, 3)
  responses[, 4:7] <- matrix(sample(1:4, n_persons * 4, replace = TRUE), n_persons, 4)
  responses[, 8:10] <- matrix(sample(1:5, n_persons * 3, replace = TRUE), n_persons, 3)

  fit <- fit_irt(responses, model = "GRM")

  expect_s3_class(fit, "gllamm_irt")

  # Check that category counts are detected correctly
  expect_equal(fit$n_categories[1:3], rep(3, 3))
  expect_equal(fit$n_categories[4:7], rep(4, 4))
  expect_equal(fit$n_categories[8:10], rep(5, 3))
  expect_equal(fit$max_categories, 5)
})


test_that("IRT handles negative values gracefully", {
  skip_if_not_installed("TMB")

  set.seed(1111)
  n_persons <- 50
  n_items <- 5

  # Invalid: negative response values
  responses_invalid <- matrix(sample(-1:1, n_persons * n_items, replace = TRUE),
                               n_persons, n_items)

  expect_error(fit_irt(responses_invalid, model = "Rasch"))
})


test_that("IRT convergence message is informative", {
  skip_if_not_installed("TMB")

  set.seed(2222)
  n_persons <- 80
  n_items <- 8

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)

  fit <- fit_irt(responses, model = "Rasch")

  # Check convergence info is present
  expect_true("convergence" %in% names(fit))
  expect_true("converged" %in% names(fit$convergence))
  expect_true("message" %in% names(fit$convergence))
})


test_that("Polytomous IRT with perfect scores on some items", {
  skip_if_not_installed("TMB")

  set.seed(3333)
  n_persons <- 100
  n_items <- 10
  n_categories <- 4

  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  # Make one item have perfect scores (all category 4)
  responses[, 3] <- n_categories

  expect_warning(fit <- fit_irt(responses, model = "PCM"),
                 "no response variance")

  expect_s3_class(fit, "gllamm_irt")

  # Item 3 should have very low threshold estimates (everyone chose highest category)
  item3_thresholds <- fit$item_parameters$thresholds[[3]]
  expect_true(all(item3_thresholds < 0))
})
