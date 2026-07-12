# E7: NPML mass-point estimation.

test_that("NPML matches npmlreg on two-point binomial heterogeneity", {
  skip_if_not_installed("npmlreg")
  skip_on_cran()  # cross-package agreement is CI-only
  set.seed(121)
  g <- 80; n_per <- 10; n <- g * n_per
  grp <- factor(rep(1:g, each = n_per))
  x <- rnorm(n)
  cls <- sample(1:2, g, TRUE, prob = c(0.6, 0.4))
  d <- data.frame(x = x, grp = grp,
                  yb = rbinom(n, 1, plogis(c(-1, 1.5)[cls[as.integer(grp)]] +
                                             0.5 * x)))

  fit <- fit_npml(yb ~ x + (1 | grp), data = d, k = 2,
                  family = stats::binomial())
  ref <- suppressMessages(
    npmlreg::allvc(yb ~ x, random = ~ 1 | grp, data = d, k = 2,
                   family = binomial(), verbose = FALSE, plot.opt = 0))

  expect_equal(unname(coef(fit)$fixed["x"]), unname(coef(ref)["x"]),
               tolerance = 5e-3)
  expect_equal(fit$locations, unname(sort(ref$mass.points)), tolerance = 1e-2)
  expect_equal(sum(fit$masses), 1, tolerance = 1e-10)
})

test_that("NPML gaussian recovers mixture structure", {
  set.seed(122)
  g <- 60; n <- 600
  grp <- factor(rep(1:g, each = n %/% g))
  x <- rnorm(n)
  cls <- sample(1:2, g, TRUE)
  d <- data.frame(x = x, grp = grp,
                  y = c(-2, 2)[cls[as.integer(grp)]] + 0.5 * x + rnorm(n))

  fit <- fit_npml(y ~ x + (1 | grp), data = d, k = 2)
  expect_true(fit$convergence$converged)
  expect_equal(fit$locations, c(-2, 2), tolerance = 0.3)
  expect_equal(unname(coef(fit)$fixed["x"]), 0.5, tolerance = 0.1)
})

test_that("NPML input validation errors clearly", {
  d <- data.frame(y = rnorm(40), x = rnorm(40), g = factor(rep(1:4, 10)))
  expect_error(fit_npml(y ~ x + (x | g), data = d, k = 2),
               "random intercepts only")
})
