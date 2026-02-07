test_that("Basic gllamm structure is correct", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(123)
  n_groups <- 10
  n_per_group <- 5
  data <- data.frame(
    y = rnorm(n_groups * n_per_group),
    x = rnorm(n_groups * n_per_group),
    group = rep(1:n_groups, each = n_per_group)
  )

  fit <- gllamm(y ~ x + (1 | group), data = data)

  expect_s3_class(fit, "gllamm")
  expect_true(all(c("coefficients", "vcov", "logLik", "AIC", "BIC") %in% names(fit)))
})


test_that("gllamm validates inputs", {
  data <- data.frame(y = 1:10, x = 1:10, group = rep(1:2, 5))

  expect_error(gllamm(data = data), "formula.*required")
  expect_error(gllamm(y ~ x + (1 | group)), "data.*required")
  expect_error(gllamm(y ~ x + (1 | group), data = "not a dataframe"), "data frame")
})


test_that("print.gllamm doesn't error", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(123)
  data <- data.frame(
    y = rnorm(50),
    x = rnorm(50),
    group = rep(1:10, each = 5)
  )

  fit <- gllamm(y ~ x + (1 | group), data = data)

  expect_output(print(fit), "Generalized Linear Latent and Mixed Model")
  expect_output(print(fit), "Random effects")
  expect_output(print(fit), "Fixed effects")
})


test_that("summary.gllamm provides detailed output", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(123)
  data <- data.frame(
    y = rnorm(50),
    x = rnorm(50),
    group = rep(1:10, each = 5)
  )

  fit <- gllamm(y ~ x + (1 | group), data = data)

  expect_output(summary(fit), "Std. Error")
  expect_output(summary(fit), "z value")
  expect_output(summary(fit), "Pr\\(>\\|z\\|\\)")
})


test_that("Extractor functions work", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(123)
  data <- data.frame(
    y = rnorm(50),
    x = rnorm(50),
    group = rep(1:10, each = 5)
  )

  fit <- gllamm(y ~ x + (1 | group), data = data)

  # fixef
  fe <- fixef(fit)
  expect_type(fe, "double")
  expect_equal(length(fe), 2)  # Intercept + x

  # ranef
  re <- ranef(fit)
  expect_type(re, "list")
  expect_equal(length(re), 10)  # 10 groups

  # VarCorr
  vc <- VarCorr(fit)
  expect_s3_class(vc, "VarCorr.gllamm")

  # coef
  cf <- coef(fit)
  expect_type(cf, "list")
  expect_true("fixed" %in% names(cf))

  # vcov
  v <- vcov(fit)
  expect_true(is.matrix(v))

  # logLik
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")

  # fitted
  fitted_vals <- fitted(fit)
  expect_equal(length(fitted_vals), 50)

  # residuals
  resids <- residuals(fit)
  expect_equal(length(resids), 50)
})


test_that("AIC and BIC are calculated", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(123)
  data <- data.frame(
    y = rnorm(50),
    x = rnorm(50),
    group = rep(1:10, each = 5)
  )

  fit <- gllamm(y ~ x + (1 | group), data = data)

  expect_true(!is.na(fit$AIC))
  expect_true(!is.na(fit$BIC))
  expect_true(fit$BIC > fit$AIC)  # BIC penalizes complexity more
})
