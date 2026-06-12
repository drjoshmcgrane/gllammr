test_that("Diagnostic plots work for GLMM", {
  skip_if_not_installed("TMB")

  set.seed(123)
  data <- data.frame(
    y = rnorm(100),
    x = rnorm(100),
    group = rep(1:10, each = 10)
  )

  fit <- gllamm(y ~ x + (1 | group), data = data)

  # Should not error
  expect_silent(plot(fit, which = 1))
  expect_silent(plot(fit, which = 2))
})


test_that("Influence diagnostics return data frame", {
  skip_if_not_installed("TMB")

  set.seed(456)
  data <- data.frame(
    y = rnorm(50),
    x = rnorm(50),
    group = rep(1:5, each = 10)
  )

  fit <- gllamm(y ~ x + (1 | group), data = data)

  infl <- influence(fit)

  expect_s3_class(infl, "data.frame")
  expect_equal(nrow(infl), 50)
  expect_true("std_residual" %in% names(infl))
})


test_that("Outlier detection finds extreme values", {
  skip_if_not_installed("TMB")

  set.seed(789)
  data <- data.frame(
    y = rnorm(100),
    x = rnorm(100),
    group = rep(1:10, each = 10)
  )

  # Add outliers
  data$y[c(1, 50, 99)] <- c(10, -10, 8)

  fit <- gllamm(y ~ x + (1 | group), data = data)

  outliers <- find_outliers(fit, threshold = 2.5)

  expect_true(!is.null(outliers))
  expect_true(nrow(outliers) >= 3)
})


test_that("Goodness of fit produces output", {
  skip_if_not_installed("TMB")

  set.seed(111)
  data <- data.frame(
    y = rnorm(80),
    x = rnorm(80),
    group = rep(1:8, each = 10)
  )

  fit <- gllamm(y ~ x + (1 | group), data = data)

  expect_output(gof.gllamm(fit), "Goodness of Fit")
  expect_output(gof.gllamm(fit), "Log-likelihood")
})


test_that("ICC calculation for Gaussian model", {
  skip_if_not_installed("TMB")

  set.seed(222)
  data <- data.frame(
    y = rnorm(100),
    x = rnorm(100),
    group = rep(1:10, each = 10)
  )

  fit <- gllamm(y ~ x + (1 | group), data = data)

  icc_vals <- icc(fit)

  expect_type(icc_vals, "double")
  expect_true(all(icc_vals >= 0 & icc_vals <= 1))
})


test_that("ICC for binomial model uses approximation", {
  skip_if_not_installed("TMB")

  set.seed(333)
  data <- data.frame(
    y = rbinom(100, 1, 0.5),
    x = rnorm(100),
    group = rep(1:10, each = 10)
  )

  fit <- gllamm(y ~ x + (1 | group), data = data, family = binomial())

  expect_message(icc(fit), "approximation")
})


test_that("Enhanced prediction on new data", {
  skip_if_not_installed("TMB")

  set.seed(444)
  train_data <- data.frame(
    y = rnorm(80),
    x = rnorm(80),
    group = rep(1:8, each = 10)
  )

  test_data <- data.frame(
    x = rnorm(20),
    group = rep(1:4, each = 5)  # Mix of old and new groups
  )

  fit <- gllamm(y ~ x + (1 | group), data = train_data)

  pred <- predict(fit, newdata = test_data)

  expect_equal(length(pred), 20)
  expect_type(pred, "double")
})
