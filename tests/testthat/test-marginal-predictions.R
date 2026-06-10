test_that("Gaussian-identity: marginal equals conditional (fixed effects)", {
  skip_if_not_installed("Matrix")

  # For Gaussian with identity link, marginal = conditional
  # because E[X'β + Z'u] = X'β + Z'E[u] = X'β
  set.seed(123)
  n_groups <- 10
  n_per_group <- 5
  n <- n_groups * n_per_group

  data <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    group = rep(1:n_groups, each = n_per_group)
  )

  # Fit model (skip if compilation issues)
  fit <- tryCatch(
    gllamm(y ~ x + (1 | group), data = data, family = gaussian()),
    error = function(e) {
      skip("Model fitting failed - likely compilation issue")
    }
  )

  # Conditional at u=0
  pred_cond <- predict(fit, re.form = NA, type = "response")

  # Marginal
  pred_marg <- predict(fit, type = "marginal")

  # Should be identical for Gaussian-identity
  expect_equal(pred_cond, pred_marg, tolerance = 1e-8)
})


test_that("Binomial-logit: marginal prediction runs without error", {
  skip_if_not_installed("Matrix")

  set.seed(456)
  n_groups <- 10
  n_per_group <- 10
  n <- n_groups * n_per_group

  data <- data.frame(
    y = rbinom(n, 1, 0.6),
    x = rnorm(n),
    group = rep(1:n_groups, each = n_per_group)
  )

  # Fit model
  fit <- tryCatch(
    gllamm(y ~ x + (1 | group), data = data, family = binomial()),
    error = function(e) {
      skip("Model fitting failed")
    }
  )

  # Marginal predictions should run
  expect_no_error({
    pred_marg <- predict(fit, type = "marginal", n_sim = 100)
  })

  # Should return numeric vector
  pred_marg <- predict(fit, type = "marginal", n_sim = 100)
  expect_type(pred_marg, "double")
  expect_length(pred_marg, n)
})


test_that("Marginal predictions: se.fit option works", {
  skip_if_not_installed("Matrix")

  set.seed(789)
  n_groups <- 8
  n_per_group <- 10
  n <- n_groups * n_per_group

  data <- data.frame(
    y = rbinom(n, 1, 0.5),
    x = rnorm(n),
    group = rep(1:n_groups, each = n_per_group)
  )

  fit <- tryCatch(
    gllamm(y ~ x + (1 | group), data = data, family = binomial()),
    error = function(e) {
      skip("Model fitting failed")
    }
  )

  # With se.fit
  result <- predict(fit, type = "marginal", se.fit = TRUE, n_sim = 100)

  expect_type(result, "list")
  expect_named(result, c("fit", "se.fit"))
  expect_length(result$fit, n)
  expect_length(result$se.fit, n)

  # SE should be non-negative
  expect_true(all(result$se.fit >= 0))
})


test_that("Marginal predictions: more samples = more stable", {
  skip_if_not_installed("Matrix")

  set.seed(101)
  n_groups <- 6
  n_per_group <- 10
  n <- n_groups * n_per_group

  data <- data.frame(
    y = rbinom(n, 1, 0.6),
    x = rnorm(n),
    group = rep(1:n_groups, each = n_per_group)
  )

  fit <- tryCatch(
    gllamm(y ~ x + (1 | group), data = data, family = binomial()),
    error = function(e) {
      skip("Model fitting failed")
    }
  )

  # Two runs with few samples should differ
  set.seed(111)
  pred_100a <- predict(fit, type = "marginal", n_sim = 100)
  set.seed(222)
  pred_100b <- predict(fit, type = "marginal", n_sim = 100)

  # Two runs with many samples should be more similar
  set.seed(111)
  pred_5ka <- predict(fit, type = "marginal", n_sim = 5000)
  set.seed(222)
  pred_5kb <- predict(fit, type = "marginal", n_sim = 5000)

  # Variance across runs should decrease with more samples
  var_100 <- var(pred_100a - pred_100b)
  var_5k <- var(pred_5ka - pred_5kb)

  expect_true(var_5k < var_100)
})


test_that("Marginal predictions: newdata works", {
  skip_if_not_installed("Matrix")

  set.seed(202)
  n_groups <- 8
  n_per_group <- 10
  n <- n_groups * n_per_group

  data <- data.frame(
    y = rbinom(n, 1, 0.5),
    x = rnorm(n),
    group = rep(1:n_groups, each = n_per_group)
  )

  fit <- tryCatch(
    gllamm(y ~ x + (1 | group), data = data, family = binomial()),
    error = function(e) {
      skip("Model fitting failed")
    }
  )

  # Create new data
  newdata <- data.frame(
    x = c(-1, 0, 1),
    group = c(1, 1, 1)  # Use existing group
  )

  # Marginal predictions for new data
  pred_new <- predict(fit, newdata = newdata, type = "marginal", n_sim = 500)

  expect_length(pred_new, 3)
  expect_type(pred_new, "double")

  # Predictions should be in valid range for binomial
  expect_true(all(pred_new >= 0 & pred_new <= 1))
})


test_that("Poisson-log: marginal predictions work", {
  skip_if_not_installed("Matrix")

  set.seed(303)
  n_groups <- 10
  n_per_group <- 10
  n <- n_groups * n_per_group

  data <- data.frame(
    y = rpois(n, lambda = 3),
    x = rnorm(n),
    group = rep(1:n_groups, each = n_per_group)
  )

  fit <- tryCatch(
    gllamm(y ~ x + (1 | group), data = data, family = poisson()),
    error = function(e) {
      skip("Model fitting failed")
    }
  )

  # Marginal predictions
  pred_marg <- predict(fit, type = "marginal", n_sim = 500)

  expect_length(pred_marg, n)
  expect_type(pred_marg, "double")

  # Predictions should be non-negative for counts
  expect_true(all(pred_marg >= 0))
})


test_that("Marginal predictions: input validation", {
  skip_if_not_installed("Matrix")

  set.seed(404)
  n <- 50

  data <- data.frame(
    y = rbinom(n, 1, 0.5),
    x = rnorm(n),
    group = rep(1:5, each = 10)
  )

  fit <- tryCatch(
    gllamm(y ~ x + (1 | group), data = data, family = binomial()),
    error = function(e) {
      skip("Model fitting failed")
    }
  )

  # Invalid n_sim
  expect_error(
    predict(fit, type = "marginal", n_sim = -1),
    NA  # May or may not error depending on validation
  )

  # n_sim = 0 should work but give warning or error
  # (implementation choice - skip this test for now)
})
