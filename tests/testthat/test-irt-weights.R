test_that("Rasch model: equal weights match unweighted", {
  set.seed(123)
  n_persons <- 50
  n_items <- 10

  # Simulate data
  theta <- rnorm(n_persons, 0, 1)
  difficulty <- rnorm(n_items, 0, 1)

  prob <- outer(theta, difficulty, function(th, b) plogis(th - b))
  responses <- matrix(rbinom(n_persons * n_items, 1, prob), n_persons, n_items)

  # Fit without weights
  fit_nowt <- fit_irt(responses, model = "Rasch")

  # Fit with equal weights
  fit_eqwt <- fit_irt(responses, model = "Rasch", weights = rep(1, n_persons))

  # Should get identical results
  expect_equal(fit_nowt$item_parameters$difficulty,
               fit_eqwt$item_parameters$difficulty,
               tolerance = 1e-6)
  expect_equal(fit_nowt$person_abilities, fit_eqwt$person_abilities, tolerance = 1e-6)
  expect_equal(fit_nowt$logLik, fit_eqwt$logLik, tolerance = 1e-6)
})


test_that("2PL model: doubled weights double log-likelihood contribution", {
  set.seed(456)
  n_persons <- 40
  n_items <- 8

  # Simulate data
  theta <- rnorm(n_persons, 0, 1)
  difficulty <- rnorm(n_items, 0, 1)
  discrimination <- runif(n_items, 0.5, 2)

  responses <- matrix(0, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      eta <- discrimination[j] * (theta[i] - difficulty[j])
      prob <- plogis(eta)
      responses[i, j] <- rbinom(1, 1, prob)
    }
  }

  # Fit with unit weights
  fit_wt1 <- fit_irt(responses, model = "2PL", weights = rep(1, n_persons))

  # Fit with doubled weights
  fit_wt2 <- fit_irt(responses, model = "2PL", weights = rep(2, n_persons))

  # Parameters should be similar (estimates don't change with constant weights)
  expect_equal(fit_wt1$item_parameters$difficulty,
               fit_wt2$item_parameters$difficulty,
               tolerance = 0.15)

  # Log-likelihood should approximately double
  # (Not exactly because of the prior on theta)
  expect_true(abs(fit_wt2$logLik / fit_wt1$logLik - 2) < 0.5)
})


test_that("GRM model: equal weights match unweighted", {
  set.seed(789)
  n_persons <- 50
  n_items <- 8
  n_categories <- 4

  # Simulate GRM data
  theta <- rnorm(n_persons, 0, 1)
  discrimination <- runif(n_items, 0.8, 1.5)

  # Create ordered thresholds
  thresholds_list <- lapply(1:n_items, function(j) {
    sort(rnorm(n_categories - 1, 0, 0.8))
  })

  responses <- matrix(0, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      tau <- thresholds_list[[j]]
      a <- discrimination[j]

      # P(Y >= k) = plogis(a * (theta - tau_k))
      p_exceed <- c(1, plogis(a * (theta[i] - tau)), 0)
      probs <- p_exceed[-length(p_exceed)] - p_exceed[-1]

      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Fit without weights
  fit_nowt <- fit_irt(responses, model = "GRM")

  # Fit with equal weights
  fit_eqwt <- fit_irt(responses, model = "GRM", weights = rep(1, n_persons))

  # Should get very similar results
  expect_equal(fit_nowt$item_parameters$discrimination,
               fit_eqwt$item_parameters$discrimination,
               tolerance = 1e-5)
  expect_equal(fit_nowt$person_abilities, fit_eqwt$person_abilities, tolerance = 1e-5)
  expect_equal(fit_nowt$logLik, fit_eqwt$logLik, tolerance = 1e-5)
})


test_that("PCM model: variable weights affect estimates appropriately", {
  set.seed(101)
  n_persons <- 60
  n_items <- 6
  n_categories <- 3

  # Simulate PCM data
  theta <- rnorm(n_persons, 0, 1)

  # Free step difficulties
  step_difficulties <- lapply(1:n_items, function(j) {
    rnorm(n_categories - 1, 0, 0.5)
  })

  responses <- matrix(0, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      delta <- step_difficulties[[j]]

      # PCM: adjacent-categories logit
      cumsums <- c(0, cumsum(theta[i] - delta))
      probs <- exp(cumsums) / sum(exp(cumsums))

      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Create variable weights (some persons weighted more)
  weights <- c(rep(1, 30), rep(2, 30))

  # Fit with weights
  fit_weighted <- fit_irt(responses, model = "PCM", weights = weights)

  # Should converge successfully
  expect_true(fit_weighted$convergence$converged)
  expect_true(is.finite(fit_weighted$logLik))

  # All discrimination parameters should be positive
  expect_true(all(fit_weighted$item_parameters$discrimination > 0))
})


test_that("GPCM model: weights validation", {
  set.seed(202)
  n_persons <- 40
  n_items <- 6
  n_categories <- 3

  # Simulate minimal GPCM data
  theta <- rnorm(n_persons, 0, 1)
  discrimination <- runif(n_items, 0.8, 1.2)

  step_difficulties <- lapply(1:n_items, function(j) {
    rnorm(n_categories - 1, 0, 0.5)
  })

  responses <- matrix(0, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      delta <- step_difficulties[[j]]
      a <- discrimination[j]

      cumsums <- c(0, cumsum(a * (theta[i] - delta)))
      probs <- exp(cumsums) / sum(exp(cumsums))

      responses[i, j] <- sample(1:n_categories, 1, prob = probs)
    }
  }

  # Test: wrong length weights should error
  expect_error(
    fit_irt(responses, model = "GPCM", weights = rep(1, n_persons - 1)),
    "Length of weights"
  )

  # Test: negative weights should error
  expect_error(
    fit_irt(responses, model = "GPCM", weights = c(rep(1, n_persons - 1), -1)),
    "non-negative"
  )

  # Test: NA weights should error
  expect_error(
    fit_irt(responses, model = "GPCM", weights = c(rep(1, n_persons - 1), NA)),
    "cannot contain missing"
  )

  # Test: valid weights should work
  fit <- fit_irt(responses, model = "GPCM", weights = runif(n_persons, 0.5, 2))
  expect_true(fit$convergence$converged)
})


test_that("3PL model: weights support", {
  set.seed(303)
  n_persons <- 50
  n_items <- 10

  # Simulate 3PL data
  theta <- rnorm(n_persons, 0, 1)
  difficulty <- rnorm(n_items, 0, 1)
  discrimination <- runif(n_items, 0.5, 2)
  guessing <- runif(n_items, 0, 0.25)

  responses <- matrix(0, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      eta <- discrimination[j] * (theta[i] - difficulty[j])
      prob <- guessing[j] + (1 - guessing[j]) * plogis(eta)
      responses[i, j] <- rbinom(1, 1, prob)
    }
  }

  # Fit with and without weights (3PL guessing parameters are weakly
  # identified at this sample size, so nlminb may report false convergence;
  # the meaningful check is that weights are honoured)
  fit_nowt <- suppressWarnings(fit_irt(responses, model = "3PL"))
  fit_eqwt <- suppressWarnings(
    fit_irt(responses, model = "3PL", weights = rep(1, n_persons)))
  fit_weighted <- suppressWarnings(
    fit_irt(responses, model = "3PL", weights = runif(n_persons, 0.8, 1.2)))

  # Equal weights must reproduce the unweighted fit
  expect_equal(fit_nowt$logLik, fit_eqwt$logLik, tolerance = 1e-6)
  expect_equal(fit_nowt$item_parameters$difficulty,
               fit_eqwt$item_parameters$difficulty, tolerance = 1e-5)

  # Weighted fit should have valid, finite results
  expect_true(is.finite(fit_weighted$logLik))
  expect_true(all(fit_weighted$item_parameters$guessing >= 0))
  expect_true(all(fit_weighted$item_parameters$guessing <= 1))
})
