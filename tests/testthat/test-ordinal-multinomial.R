test_that("Ordinal model accepts valid input", {
  skip_if_not_installed("TMB")

  set.seed(123)
  n <- 100
  data <- data.frame(
    y = factor(sample(1:5, n, replace = TRUE), ordered = TRUE),
    x = rnorm(n),
    group = rep(1:10, each = 10)
  )

  fit <- fit_ordinal(y ~ x + (1 | group), data = data)

  expect_s3_class(fit, "gllamm_ordinal")
  expect_equal(fit$n_categories, 5)
})


test_that("Ordinal threshold parameters are ordered", {
  skip_if_not_installed("TMB")

  set.seed(456)
  n <- 200
  data <- data.frame(
    y = ordered(sample(1:4, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:20, each = 10)
  )

  fit <- fit_ordinal(y ~ x + (1 | group), data = data, link = "logit")

  # Thresholds should be strictly increasing
  thresholds <- fit$coefficients$thresholds
  expect_true(all(diff(thresholds) > 0))
})


test_that("Proportional odds vs probit link", {
  skip_if_not_installed("TMB")

  set.seed(789)
  n <- 150
  data <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:15, each = 10)
  )

  fit_logit <- fit_ordinal(y ~ x + (1 | group), data = data, link = "logit")
  fit_probit <- fit_ordinal(y ~ x + (1 | group), data = data, link = "probit")

  expect_equal(fit_logit$link, "logit")
  expect_equal(fit_probit$link, "probit")
})


test_that("Multinomial model with 3 categories", {
  skip_if_not_installed("TMB")

  set.seed(111)
  n <- 120
  data <- data.frame(
    y = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:12, each = 10)
  )

  fit <- fit_multinomial(y ~ x + (1 | group), data = data)

  expect_s3_class(fit, "gllamm_multinomial")
  expect_equal(fit$n_categories, 3)
  expect_equal(nrow(fit$coefficients$beta), 2)  # K-1 rows
})


test_that("Multinomial reference category", {
  skip_if_not_installed("TMB")

  set.seed(222)
  data <- data.frame(
    y = factor(rep(c("Low", "Med", "High"), each = 30)),
    x = rnorm(90),
    group = rep(1:9, each = 10)
  )

  fit <- fit_multinomial(y ~ x + (1 | group), data = data)

  # First level should be reference
  expect_equal(fit$reference, "High")  # Alphabetically first
})


test_that("Ordinal print and summary methods", {
  skip_if_not_installed("TMB")

  set.seed(333)
  data <- data.frame(
    y = ordered(sample(1:4, 80, replace = TRUE)),
    x = rnorm(80),
    group = rep(1:8, each = 10)
  )

  fit <- fit_ordinal(y ~ x + (1 | group), data = data)

  expect_output(print(fit), "Ordinal Regression")
  expect_output(print(fit), "Threshold parameters")
})


test_that("Multinomial print method", {
  skip_if_not_installed("TMB")

  set.seed(444)
  data <- data.frame(
    y = factor(sample(c("A", "B", "C"), 90, replace = TRUE)),
    x = rnorm(90),
    group = rep(1:9, each = 10)
  )

  fit <- fit_multinomial(y ~ x + (1 | group), data = data)

  expect_output(print(fit), "Multinomial Regression")
  expect_output(print(fit), "Reference category")
})
