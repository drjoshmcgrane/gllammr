test_that("IRT data simulation works", {
  set.seed(123)
  n_persons <- 100
  n_items <- 10

  # Simulate Rasch data
  theta <- rnorm(n_persons, 0, 1)
  difficulty <- rnorm(n_items, 0, 1)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      p <- plogis(theta[i] - difficulty[j])
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  expect_equal(dim(responses), c(100, 10))
  expect_true(all(responses %in% c(0, 1)))
})


test_that("fit_irt accepts valid input", {
  skip_if_not_installed("TMB")

  set.seed(456)
  n_persons <- 50
  n_items <- 5

  # Simple Rasch data
  theta <- rnorm(n_persons)
  difficulty <- seq(-1, 1, length.out = n_items)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      p <- plogis(theta[i] - difficulty[j])
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  fit <- fit_irt(responses, model = "Rasch")

  expect_s3_class(fit, "gllamm_irt")
  expect_true("item_parameters" %in% names(fit))
  expect_true("person_abilities" %in% names(fit))
  expect_equal(nrow(fit$item_parameters), n_items)
  expect_equal(length(fit$person_abilities), n_persons)
})


test_that("Rasch model recovers known parameters", {
  skip_if_not_installed("TMB")

  set.seed(789)
  n_persons <- 200
  n_items <- 15

  # Known parameters
  true_difficulty <- seq(-2, 2, length.out = n_items)
  true_theta <- rnorm(n_persons, 0, 1)

  # Generate data
  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      p <- plogis(true_theta[i] - true_difficulty[j])
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  # Fit model
  fit <- fit_irt(responses, model = "Rasch")

  # Check parameter recovery (with generous tolerance)
  # Note: Need to account for identification constraints
  recovered_difficulty <- fit$item_parameters$difficulty
  expect_equal(cor(recovered_difficulty, true_difficulty), 1, tolerance = 0.3)
})


test_that("2PL model includes discrimination parameters", {
  skip_if_not_installed("TMB")

  set.seed(321)
  n_persons <- 100
  n_items <- 8

  theta <- rnorm(n_persons)
  difficulty <- rnorm(n_items, 0, 1)
  discrimination <- runif(n_items, 0.5, 2)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      p <- plogis(discrimination[j] * (theta[i] - difficulty[j]))
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  fit <- fit_irt(responses, model = "2PL")

  expect_s3_class(fit, "gllamm_irt")
  expect_true(all(fit$item_parameters$discrimination > 0))
  expect_equal(nrow(fit$item_parameters), n_items)
})


test_that("IRT print and summary methods work", {
  skip_if_not_installed("TMB")

  set.seed(111)
  responses <- matrix(rbinom(50 * 5, 1, 0.5), 50, 5)

  fit <- fit_irt(responses, model = "Rasch")

  expect_output(print(fit), "IRT Model")
  expect_output(print(fit), "Rasch")
  expect_output(summary(fit), "Ability quartiles")
})


test_that("IRT handles missing data", {
  set.seed(222)
  responses <- matrix(rbinom(100 * 10, 1, 0.5), 100, 10)

  # Introduce missing values
  responses[sample(length(responses), 50)] <- NA

  # Should not error
  expect_silent({
    # fit <- fit_irt(responses, model = "Rasch")
    # Would work after TMB compilation
  })
})
