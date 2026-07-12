# Parametric survival (frailty) models: exponential and Weibull.

simulate_survival <- function(seed = 51, n = 1000, g = 40, shape = 1) {
  set.seed(seed)
  grp <- factor(rep(1:g, each = n %/% g))
  x <- rnorm(n)
  u <- rnorm(g, 0, 0.6)
  lam <- exp(-1 + 0.5 * x + u[as.integer(grp)])
  t_true <- (rexp(n))^(1 / shape) / lam
  cens <- quantile(t_true, 0.85)
  data.frame(time = pmin(t_true, cens), status = as.integer(t_true <= cens),
             x = x, grp = grp)
}

test_that("exponential frailty equals Poisson GLMM with log-time offset", {
  skip_if_not_installed("lme4")
  d <- simulate_survival()
  fit <- fit_survival(Surv(time, status) ~ x + (1 | grp), data = d,
                      distribution = "exponential")
  ref <- ref_fit(lme4::glmer(status ~ x + offset(log(time)) + (1 | grp),
                             data = d, family = stats::poisson(), nAGQ = 1,
                             control = lme4::glmerControl(optimizer = "bobyqa")))

  expect_equal(unname(coef(fit)$fixed), unname(lme4::fixef(ref)), tolerance = 1e-3)
  expect_equal(unname(fit$coefficients$random_sd),
               unname(attr(lme4::VarCorr(ref)$grp, "stddev")), tolerance = 1e-3)
  expect_equal(fit$logLik,
               as.numeric(logLik(ref)) - sum(d$status * log(d$time)),
               tolerance = 1e-2)
})

test_that("Weibull frailty recovers shape and coefficients", {
  d <- simulate_survival(seed = 52, shape = 1.6)
  fit <- fit_survival(Surv(time, status) ~ x + (1 | grp), data = d)

  expect_true(fit$convergence$converged)
  expect_equal(unname(fit$shape), 1.6, tolerance = 0.15)
  expect_equal(unname(coef(fit)$fixed[2]), 0.5, tolerance = 0.15)
  expect_equal(unname(fit$coefficients$random_sd), 0.6, tolerance = 0.2)
  expect_s3_class(fit, "gllamm_survival")
})

test_that("survival input validation errors clearly", {
  d <- simulate_survival(n = 100, g = 5)
  expect_error(fit_survival(time ~ x + (1 | grp), data = d), "Surv")
  d$bad <- -d$time
  expect_error(fit_survival(Surv(bad, status) ~ x + (1 | grp), data = d),
               "positive")
})
