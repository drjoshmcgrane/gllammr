test_that("binomial() family constructor works", {
  # Test logit link (default)
  fam1 <- binomial()
  expect_s3_class(fam1, "binomial_family")
  expect_equal(fam1$family, "binomial")
  expect_equal(fam1$link, "logit")
  expect_equal(fam1$link_code, 1L)

  # Test probit link
  fam2 <- binomial(link = "probit")
  expect_equal(fam2$link, "probit")
  expect_equal(fam2$link_code, 2L)

  # Test cloglog link
  fam3 <- binomial(link = "cloglog")
  expect_equal(fam3$link, "cloglog")
  expect_equal(fam3$link_code, 3L)

  # Test invalid link
  expect_error(binomial(link = "invalid"))
})


test_that("fit_binomial() works with logit link", {
  skip_on_cran()

  set.seed(123)
  n_groups <- 10
  n_per_group <- 20
  n <- n_groups * n_per_group

  # Simulate data
  group <- rep(1:n_groups, each = n_per_group)
  x <- rnorm(n)
  group_effect <- rep(rnorm(n_groups, 0, 0.5), each = n_per_group)
  eta <- -0.5 + 1.0 * x + group_effect
  p <- plogis(eta)
  y <- rbinom(n, 1, p)

  data <- data.frame(y = y, x = x, group = group)

  # Fit model
  fit <- fit_binomial(y ~ x + (1 | group), data = data, link = "logit")

  # Check class
  expect_s3_class(fit, "gllamm_binomial")
  expect_s3_class(fit, "gllamm")

  # Check components
  expect_true(!is.null(fit$coefficients))
  expect_true(!is.null(fit$logLik))
  expect_true(!is.null(fit$AIC))
  expect_true(!is.null(fit$BIC))
  expect_equal(fit$link, "logit")
  expect_equal(fit$n_obs, n)
  expect_equal(fit$n_groups, n_groups)

  # Check convergence
  expect_true(fit$convergence$converged)

  # Check parameter recovery (approximately)
  expect_equal(unname(fit$coefficients$fixed[2]), 1.0, tolerance = 0.3)  # Slope coefficient
})


test_that("fit_binomial() works with probit link", {
  skip_on_cran()

  set.seed(456)
  n_groups <- 10
  n_per_group <- 20
  n <- n_groups * n_per_group

  # Simulate data
  group <- rep(1:n_groups, each = n_per_group)
  x <- rnorm(n)
  group_effect <- rep(rnorm(n_groups, 0, 0.5), each = n_per_group)
  eta <- -0.3 + 0.8 * x + group_effect
  p <- pnorm(eta)
  y <- rbinom(n, 1, p)

  data <- data.frame(y = y, x = x, group = group)

  # Fit model
  fit <- fit_binomial(y ~ x + (1 | group), data = data, link = "probit")

  # Check link
  expect_equal(fit$link, "probit")

  # Check convergence
  expect_true(fit$convergence$converged)

  # Check fitted values are in (0,1)
  expect_true(all(fit$fitted_values >= 0 & fit$fitted_values <= 1))
})


test_that("fit_binomial() works with cloglog link", {
  skip_on_cran()

  set.seed(789)
  n_groups <- 10
  n_per_group <- 20
  n <- n_groups * n_per_group

  # Simulate data with cloglog link
  group <- rep(1:n_groups, each = n_per_group)
  x <- rnorm(n)
  group_effect <- rep(rnorm(n_groups, 0, 0.3), each = n_per_group)
  eta <- -1.5 + 0.7 * x + group_effect
  p <- 1 - exp(-exp(eta))  # Complementary log-log
  y <- rbinom(n, 1, p)

  data <- data.frame(y = y, x = x, group = group)

  # Fit model
  fit <- fit_binomial(y ~ x + (1 | group), data = data, link = "cloglog")

  # Check link
  expect_equal(fit$link, "cloglog")

  # Check convergence
  expect_true(fit$convergence$converged)

  # Check fitted values are in (0,1)
  expect_true(all(fit$fitted_values >= 0 & fit$fitted_values <= 1))

  # cloglog should recover parameters reasonably
  expect_equal(unname(fit$coefficients$fixed[2]), 0.7, tolerance = 0.4)
})


test_that("gllamm() works with binomial() family", {
  skip_on_cran()

  set.seed(111)
  n_groups <- 8
  n_per_group <- 15
  n <- n_groups * n_per_group

  # Simulate data
  group <- rep(1:n_groups, each = n_per_group)
  x <- rnorm(n)
  group_effect <- rep(rnorm(n_groups, 0, 0.5), each = n_per_group)
  eta <- 0.2 + 0.9 * x + group_effect
  p <- plogis(eta)
  y <- rbinom(n, 1, p)

  data <- data.frame(y = y, x = x, group = group)

  # Fit via gllamm() with binomial() family
  fit1 <- gllamm(y ~ x + (1 | group),
                 data = data,
                 family = binomial(link = "logit"))

  # Check class
  expect_s3_class(fit1, "gllamm_binomial")

  # Should be equivalent to fit_binomial()
  fit2 <- fit_binomial(y ~ x + (1 | group), data = data, link = "logit")

  expect_equal(fit1$logLik, fit2$logLik, tolerance = 1e-6)
  expect_equal(fit1$coefficients$fixed, fit2$coefficients$fixed, tolerance = 1e-6)
})


test_that("gllamm() works with binomial(link = 'cloglog')", {
  skip_on_cran()

  set.seed(222)
  n <- 100

  # Simulate rare event data (cloglog appropriate)
  group <- rep(1:10, each = 10)
  x <- rnorm(n)
  group_effect <- rep(rnorm(10, 0, 0.2), each = 10)
  eta <- -2.0 + 1.0 * x + group_effect
  p <- 1 - exp(-exp(eta))
  y <- rbinom(n, 1, p)

  data <- data.frame(y = y, x = x, group = group)

  # Fit with cloglog
  fit <- gllamm(y ~ x + (1 | group),
                data = data,
                family = binomial(link = "cloglog"))

  # Check link
  expect_equal(fit$link, "cloglog")
  expect_true(fit$convergence$converged)
})


test_that("print.gllamm_binomial() works", {
  skip_on_cran()

  set.seed(333)
  n <- 80
  group <- rep(1:8, each = 10)
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(x))

  data <- data.frame(y = y, x = x, group = group)
  fit <- fit_binomial(y ~ x + (1 | group), data = data)

  # Should not error
  expect_output(print(fit), "Binomial Regression Model")
  expect_output(print(fit), "Link function: logit")
})


test_that("summary.gllamm_binomial() works", {
  skip_on_cran()

  set.seed(444)
  n <- 80
  group <- rep(1:8, each = 10)
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(x))

  data <- data.frame(y = y, x = x, group = group)
  fit <- fit_binomial(y ~ x + (1 | group), data = data)

  # Should not error
  expect_output(summary(fit), "Binomial Regression Model")
  expect_output(summary(fit), "Fixed Effects")
  expect_output(summary(fit), "Random Effects")
})


test_that("binomial model rejects non-binary response", {
  data <- data.frame(
    y = c(0, 1, 2, 3),  # Not binary
    x = rnorm(4),
    group = c(1, 1, 2, 2)
  )

  expect_error(
    fit_binomial(y ~ x + (1 | group), data = data),
    "binary"
  )
})
