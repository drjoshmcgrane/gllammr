# E2 (joint mixed-response) and E3 (SEM) models.

test_that("SEM with structural path matches lavaan", {
  skip_if_not_installed("lavaan")
  set.seed(71)
  n <- 600
  f1 <- rnorm(n); f2 <- 0.6 * f1 + rnorm(n, 0, 0.8)
  d <- data.frame(
    x1 = 1.0 * f1 + rnorm(n, 0, 0.6), x2 = 0.8 * f1 + rnorm(n, 0, 0.6),
    x3 = 1.2 * f1 + rnorm(n, 0, 0.6),
    y1 = 1.0 * f2 + rnorm(n, 0, 0.5), y2 = 0.9 * f2 + rnorm(n, 0, 0.5),
    y3 = 1.1 * f2 + rnorm(n, 0, 0.5))

  fit <- fit_sem(measurement = list(f1 = ~ x1 + x2 + x3, f2 = ~ y1 + y2 + y3),
                 structural = list(f2 ~ f1), data = d)
  lav <- lavaan::sem("f1 =~ x1 + x2 + x3\nf2 =~ y1 + y2 + y3\nf2 ~ f1",
                     data = d)
  pe <- lavaan::parameterEstimates(lav)

  expect_true(fit$convergence$converged)
  expect_equal(fit$loadings["x2", "f1"],
               pe$est[pe$lhs == "f1" & pe$op == "=~" & pe$rhs == "x2"],
               tolerance = 5e-3)
  expect_equal(fit$structural["f2", "f1"], pe$est[pe$op == "~"],
               tolerance = 5e-3)
})

test_that("SEM rejects cyclic structural models", {
  d <- data.frame(x1 = rnorm(50), x2 = rnorm(50), y1 = rnorm(50), y2 = rnorm(50))
  expect_error(
    fit_sem(measurement = list(f1 = ~ x1 + x2, f2 = ~ y1 + y2),
            structural = list(f2 ~ f1, f1 ~ f2), data = d),
    "cycle")
})

test_that("joint mixed-response model recovers parameters", {
  set.seed(61)
  n <- 1500; g <- 60
  grp <- factor(rep(1:g, each = n %/% g))
  x <- rnorm(n)
  u <- rnorm(g, 0, 0.8)
  d <- data.frame(
    x = x, grp = grp,
    yc = 1 + 0.5 * x + u[as.integer(grp)] + rnorm(n, 0, 1.2),
    yb = rbinom(n, 1, plogis(-0.5 + 0.7 * x + u[as.integer(grp)])))

  fit <- fit_mixed(formulas = list(gaussian = yc ~ x, binomial = yb ~ x),
                   random = ~ (1 | grp), data = d)

  expect_true(fit$convergence$converged)
  expect_equal(unname(fit$coefficients$gaussian["x"]), 0.5, tolerance = 0.1)
  expect_equal(unname(fit$coefficients$binomial["x"]), 0.7, tolerance = 0.15)
  expect_equal(fit$residual_sd, 1.2, tolerance = 0.1)
  expect_equal(fit$random_sd, 0.8, tolerance = 0.25)
})

test_that("mixed-response input validation errors clearly", {
  d <- data.frame(y = rnorm(50), g = factor(rep(1:5, 10)))
  expect_error(fit_mixed(list(weibull = y ~ 1), ~ (1 | g), d), "names among")
})
