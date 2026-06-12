# Random slopes in GLMMs: (x | g) and (x || g) for gaussian/binomial/poisson.
# lme4 with nAGQ = 1 uses the same Laplace approximation, so agreement is tight.

simulate_slopes_data <- function(seed = 11, g = 30, n_per = 20) {
  set.seed(seed)
  n <- g * n_per
  grp <- factor(rep(1:g, each = n_per))
  Sig <- matrix(c(1.0, 0.3, 0.3, 0.49), 2)
  U <- matrix(rnorm(g * 2), g) %*% chol(Sig)
  x <- rnorm(n)
  eta <- (1 + U[as.integer(grp), 1]) + (0.5 + U[as.integer(grp), 2]) * x
  d <- data.frame(
    x = x, grp = grp, eta = eta,
    y  = eta + rnorm(n, 0, 0.8),
    yb = rbinom(n, 1, plogis(eta - 1))
  )
  # Own seed so the poisson response does not depend on the draws above
  set.seed(seed + 1)
  d$yp <- rpois(n, exp(0.2 + 0.4 * U[as.integer(grp), 1] +
                         (0.3 + 0.4 * U[as.integer(grp), 2]) * x))
  d
}

test_that("gaussian random slopes fit and return full covariance", {
  d <- simulate_slopes_data()
  fit <- gllamm(y ~ x + (x | grp), data = d)

  expect_true(fit$convergence$converged)
  expect_length(coef(fit)$fixed, 2)
  Sigma <- fit$coefficients$random_var[[1]]
  expect_equal(dim(Sigma), c(2, 2))
  expect_true(Sigma[1, 2] != 0)          # correlation estimated
  expect_equal(Sigma[1, 2], Sigma[2, 1]) # symmetric
})

test_that("gaussian random slopes match lme4 (same Laplace approximation)", {
  skip_if_not_installed("lme4")
  d <- simulate_slopes_data()

  fit <- gllamm(y ~ x + (x | grp), data = d)
  ref <- lme4::lmer(y ~ x + (x | grp), data = d, REML = FALSE)

  expect_equal(unname(coef(fit)$fixed), unname(lme4::fixef(ref)), tolerance = 1e-4)
  expect_equal(fit$logLik, as.numeric(logLik(ref)), tolerance = 1e-4)
  Sigma <- fit$coefficients$random_var[[1]]
  Sigma_ref <- as.matrix(Matrix::bdiag(lme4::VarCorr(ref)$grp))
  expect_equal(unname(Sigma), unname(Sigma_ref), tolerance = 1e-3)
})

test_that("binomial random slopes match lme4", {
  skip_if_not_installed("lme4")
  d <- simulate_slopes_data()

  fit <- gllamm(yb ~ x + (x | grp), data = d, family = stats::binomial())
  ref <- lme4::glmer(yb ~ x + (x | grp), data = d, family = stats::binomial())

  expect_equal(unname(coef(fit)$fixed), unname(lme4::fixef(ref)), tolerance = 1e-3)
  expect_equal(fit$logLik, as.numeric(logLik(ref)), tolerance = 1e-3)
  Sigma <- fit$coefficients$random_var[[1]]
  Sigma_ref <- as.matrix(Matrix::bdiag(lme4::VarCorr(ref)$grp))
  expect_equal(unname(Sigma), unname(Sigma_ref), tolerance = 5e-2)
})

test_that("poisson random slopes match lme4", {
  skip_if_not_installed("lme4")
  d <- simulate_slopes_data()

  fit <- gllamm(yp ~ x + (x | grp), data = d, family = stats::poisson())
  # bobyqa with a generous budget: the default optimizer can stop short on
  # some datasets, which would test the optimizers rather than the likelihood
  ref <- lme4::glmer(yp ~ x + (x | grp), data = d, family = stats::poisson(),
                     control = lme4::glmerControl(optimizer = "bobyqa",
                                                  optCtrl = list(maxfun = 1e5)))

  expect_equal(unname(coef(fit)$fixed), unname(lme4::fixef(ref)), tolerance = 1e-3)
  expect_equal(fit$logLik, as.numeric(logLik(ref)), tolerance = 1e-3)
})

test_that("uncorrelated slopes (x || g) force zero covariance", {
  d <- simulate_slopes_data()
  fit <- gllamm(yb ~ x + (x || grp), data = d, family = stats::binomial())

  Sigma <- fit$coefficients$random_var[[1]]
  expect_equal(Sigma[1, 2], 0)
  expect_true(all(diag(Sigma) > 0))
})

test_that("binomial constructor path returns full covariance for slopes", {
  d <- simulate_slopes_data()
  fit <- gllamm(yb ~ x + (x | grp), data = d, family = binomial())

  Sigma <- fit$coefficients$random_var[[1]]
  expect_equal(dim(Sigma), c(2, 2))
})
