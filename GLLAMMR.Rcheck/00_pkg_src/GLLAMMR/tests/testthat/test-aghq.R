# E8: adaptive Gauss-Hermite quadrature.

simulate_hard_binary <- function(seed = 131, g = 100, n_per = 6) {
  set.seed(seed)
  n <- g * n_per
  grp <- factor(rep(1:g, each = n_per))
  x <- rnorm(n)
  u <- rnorm(g, 0, 2)
  data.frame(x = x, grp = grp,
             yb = rbinom(n, 1, plogis(-0.5 + 0.8 * x + u[as.integer(grp)])))
}

test_that("aghq(15) matches glmer nAGQ=15 where Laplace is biased", {
  skip_if_not_installed("lme4")
  d <- simulate_hard_binary()

  fit <- gllamm(yb ~ x + (1 | grp), data = d, family = stats::binomial(),
                integration = aghq(15))
  ref <- lme4::glmer(yb ~ x + (1 | grp), data = d,
                     family = stats::binomial(), nAGQ = 15)

  expect_equal(unname(coef(fit)$fixed), unname(lme4::fixef(ref)),
               tolerance = 2e-3)
  expect_equal(sqrt(fit$coefficients$random_var[[1]][1, 1]),
               unname(attr(lme4::VarCorr(ref)$grp, "stddev")),
               tolerance = 5e-3)
  expect_equal(fit$logLik, as.numeric(logLik(ref)), tolerance = 0.05)
})

test_that("aghq improves on Laplace toward the quadrature reference", {
  skip_if_not_installed("lme4")
  d <- simulate_hard_binary(seed = 132)

  fit_lap <- gllamm(yb ~ x + (1 | grp), data = d, family = stats::binomial())
  fit_aghq <- gllamm(yb ~ x + (1 | grp), data = d, family = stats::binomial(),
                     integration = aghq(15))
  ref <- lme4::glmer(yb ~ x + (1 | grp), data = d,
                     family = stats::binomial(), nAGQ = 25)

  s_lap <- sqrt(fit_lap$coefficients$random_var[[1]][1, 1])
  s_aghq <- sqrt(fit_aghq$coefficients$random_var[[1]][1, 1])
  s_ref <- unname(attr(lme4::VarCorr(ref)$grp, "stddev"))

  expect_lt(abs(s_aghq - s_ref), abs(s_lap - s_ref))
})

test_that("aghq gaussian agrees with Laplace (both exact for gaussian)", {
  set.seed(133)
  g <- 40; n <- 800
  grp <- factor(rep(1:g, each = n %/% g))
  x <- rnorm(n)
  u <- rnorm(g, 0, 1)
  d <- data.frame(y = 1 + 0.5 * x + u[as.integer(grp)] + rnorm(n),
                  x = x, grp = grp)
  fit_lap <- gllamm(y ~ x + (1 | grp), data = d)
  fit_aghq <- gllamm(y ~ x + (1 | grp), data = d, integration = aghq(7))
  expect_equal(unname(coef(fit_aghq)$fixed), unname(coef(fit_lap)$fixed),
               tolerance = 1e-3)
  expect_equal(fit_aghq$logLik, fit_lap$logLik, tolerance = 0.01)
})

test_that("aghq input validation errors clearly", {
  expect_error(aghq(1), ">= 2")
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
  expect_error(gllamm(y ~ x + (x | g), data = d, integration = aghq(7)),
               "single random intercept")
})
