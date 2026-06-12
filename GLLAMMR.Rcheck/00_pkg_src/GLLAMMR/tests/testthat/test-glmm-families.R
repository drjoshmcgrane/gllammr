test_that("Binomial GLMM with logit link", {
  skip_if_not_installed("TMB")

  set.seed(123)
  n_groups <- 20
  n_per_group <- 10

  # Generate binary data
  group <- rep(1:n_groups, each = n_per_group)
  x <- rnorm(n_groups * n_per_group)
  u <- rnorm(n_groups, sd = 0.8)

  eta <- -0.5 + 0.7 * x + u[group]
  p <- plogis(eta)
  y <- rbinom(n_groups * n_per_group, 1, p)

  data <- data.frame(y = y, x = x, group = group)

  fit <- gllamm(y ~ x + (1 | group),
                data = data,
                family = binomial(link = "logit"))

  expect_s3_class(fit, "gllamm")
  expect_equal(fit$family$family, "binomial")
  expect_equal(fit$family$link, "logit")
})


test_that("Binomial GLMM recovers parameters", {
  skip_if_not_installed("TMB")

  set.seed(456)
  n_groups <- 60
  n_per_group <- 15

  # Known parameters
  beta0 <- 0.3
  beta1 <- -0.6
  sigma_u <- 1.0

  group <- rep(1:n_groups, each = n_per_group)
  x <- rnorm(n_groups * n_per_group)
  u <- rnorm(n_groups, sd = sigma_u)

  eta <- beta0 + beta1 * x + u[group]
  p <- plogis(eta)
  y <- rbinom(n_groups * n_per_group, 1, p)

  data <- data.frame(y = y, x = x, group = group)

  fit <- gllamm(y ~ x + (1 | group),
                data = data,
                family = binomial())

  # Check parameter recovery (generous tolerance)
  expect_equal(unname(fixef(fit)["(Intercept)"]), beta0, tolerance = 0.3)
  expect_equal(unname(fixef(fit)["x"]), beta1, tolerance = 0.3)
  expect_equal(sqrt(VarCorr(fit)[[1]][1, 1]), sigma_u, tolerance = 0.4)
})


test_that("Binomial GLMM with probit link", {
  skip_if_not_installed("TMB")

  set.seed(789)
  n_groups <- 15
  n_per_group <- 20

  group <- rep(1:n_groups, each = n_per_group)
  x <- rnorm(n_groups * n_per_group)
  u <- rnorm(n_groups, sd = 0.5)

  eta <- 0.2 + 0.5 * x + u[group]
  p <- pnorm(eta)
  y <- rbinom(n_groups * n_per_group, 1, p)

  data <- data.frame(y = y, x = x, group = group)

  fit <- gllamm(y ~ x + (1 | group),
                data = data,
                family = binomial(link = "probit"))

  expect_s3_class(fit, "gllamm")
  expect_equal(fit$family$link, "probit")
})


test_that("Poisson GLMM with log link", {
  skip_if_not_installed("TMB")

  set.seed(111)
  n_groups <- 20
  n_per_group <- 12

  # Generate count data
  group <- rep(1:n_groups, each = n_per_group)
  x <- rnorm(n_groups * n_per_group)
  u <- rnorm(n_groups, sd = 0.3)

  eta <- 1.0 + 0.4 * x + u[group]
  lambda <- exp(eta)
  y <- rpois(n_groups * n_per_group, lambda)

  data <- data.frame(y = y, x = x, group = group)

  fit <- gllamm(y ~ x + (1 | group),
                data = data,
                family = poisson())

  expect_s3_class(fit, "gllamm")
  expect_equal(fit$family$family, "poisson")
  expect_equal(fit$family$link, "log")
})


test_that("Poisson GLMM recovers parameters", {
  skip_if_not_installed("TMB")

  set.seed(222)
  n_groups <- 25
  n_per_group <- 20

  # Known parameters
  beta0 <- 0.5
  beta1 <- 0.3
  sigma_u <- 0.4

  group <- rep(1:n_groups, each = n_per_group)
  x <- rnorm(n_groups * n_per_group)
  u <- rnorm(n_groups, sd = sigma_u)

  eta <- beta0 + beta1 * x + u[group]
  lambda <- exp(eta)
  y <- rpois(n_groups * n_per_group, lambda)

  data <- data.frame(y = y, x = x, group = group)

  fit <- gllamm(y ~ x + (1 | group),
                data = data,
                family = poisson())

  # Check parameter recovery
  expect_equal(unname(fixef(fit)["(Intercept)"]), beta0, tolerance = 0.2)
  expect_equal(unname(fixef(fit)["x"]), beta1, tolerance = 0.2)
  expect_equal(sqrt(VarCorr(fit)[[1]][1, 1]), sigma_u, tolerance = 0.3)
})


test_that("GLMMs handle overdispersion", {
  skip_if_not_installed("TMB")

  set.seed(333)
  n_groups <- 30
  n_per_group <- 15

  # Create overdispersed Poisson data
  group <- rep(1:n_groups, each = n_per_group)
  x <- rnorm(n_groups * n_per_group)
  u <- rnorm(n_groups, sd = 0.6)  # Large random effect = overdispersion

  eta <- 2.0 + 0.2 * x + u[group]
  lambda <- exp(eta)
  y <- rpois(n_groups * n_per_group, lambda)

  data <- data.frame(y = y, x = x, group = group)

  fit <- gllamm(y ~ x + (1 | group),
                data = data,
                family = poisson())

  # Random effect variance should be positive
  expect_true(VarCorr(fit)[[1]][1, 1] > 0)
})


test_that("Binomial fitted values are probabilities", {
  skip_if_not_installed("TMB")

  set.seed(444)
  n_groups <- 10
  n_per_group <- 20

  group <- rep(1:n_groups, each = n_per_group)
  x <- rnorm(n_groups * n_per_group)
  u <- rnorm(n_groups, sd = 0.5)

  eta <- 0 + 0.5 * x + u[group]
  p <- plogis(eta)
  y <- rbinom(n_groups * n_per_group, 1, p)

  data <- data.frame(y = y, x = x, group = group)

  fit <- gllamm(y ~ x + (1 | group),
                data = data,
                family = binomial())

  fitted_vals <- fitted(fit)

  # Fitted values should be probabilities
  expect_true(all(fitted_vals >= 0))
  expect_true(all(fitted_vals <= 1))
})


test_that("Poisson fitted values are positive", {
  skip_if_not_installed("TMB")

  set.seed(555)
  n_groups <- 10
  n_per_group <- 15

  group <- rep(1:n_groups, each = n_per_group)
  x <- rnorm(n_groups * n_per_group)
  u <- rnorm(n_groups, sd = 0.3)

  eta <- 1 + 0.3 * x + u[group]
  lambda <- exp(eta)
  y <- rpois(n_groups * n_per_group, lambda)

  data <- data.frame(y = y, x = x, group = group)

  fit <- gllamm(y ~ x + (1 | group),
                data = data,
                family = poisson())

  fitted_vals <- fitted(fit)

  # Fitted values should be positive
  expect_true(all(fitted_vals > 0))
})
