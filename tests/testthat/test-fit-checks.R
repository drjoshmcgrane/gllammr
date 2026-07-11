# Shared robustness helpers: check_sdreport, warn_not_converged,
# safe_nlminb, safe_solve, and safe_chol (via rmvnorm_chol).

test_that("check_sdreport flags a try-error with a warning and se_ok FALSE", {
  fake_err <- structure("boom", class = "try-error",
                        condition = simpleError("boom"))
  expect_warning(res <- check_sdreport(fake_err, "test model"),
                 "Standard error computation failed for test model")
  expect_null(res$sdr)
  expect_false(res$se_ok)
})

test_that("check_sdreport flags a non-positive-definite Hessian", {
  fake <- list(pdHess = FALSE)
  expect_warning(res <- check_sdreport(fake, "test model"),
                 "Hessian not positive definite for test model")
  expect_identical(res$sdr, fake)
  expect_false(res$se_ok)
})

test_that("check_sdreport is silent on a good sdreport", {
  fake <- list(pdHess = TRUE)
  expect_silent(res <- check_sdreport(fake, "test model"))
  expect_identical(res$sdr, fake)
  expect_true(res$se_ok)
})

test_that("check_sdreport is silent when SEs were not requested (NULL)", {
  expect_silent(res <- check_sdreport(NULL, "test model"))
  expect_null(res$sdr)
  expect_false(res$se_ok)
})

test_that("warn_not_converged warns only when convergence is FALSE", {
  # Nested convergence field (the common package layout)
  expect_warning(warn_not_converged(list(convergence = list(converged = FALSE))),
                 "did not converge")
  # Top-level field
  expect_warning(warn_not_converged(list(converged = FALSE)),
                 "did not converge")
  # Converged, missing, or NULL: silent
  expect_silent(warn_not_converged(list(convergence = list(converged = TRUE))))
  expect_silent(warn_not_converged(list()))
  expect_silent(warn_not_converged(list(converged = NULL)))
})

test_that("safe_nlminb returns nlminb shape on success", {
  opt <- safe_nlminb(c(1, 1), function(p) sum((p - 3)^2), context = "test")
  expect_equal(opt$par, c(3, 3), tolerance = 1e-5)
  expect_true(is.finite(opt$objective))
})

test_that("safe_nlminb errors informatively on a non-finite objective", {
  # nlminb itself emits an "NA/NaN function evaluation" warning; suppress it
  # so only the informative error is exercised.
  expect_error(
    suppressWarnings(safe_nlminb(1, function(p) NaN, context = "test model")),
    "Optimization failed for test model"
  )
})

test_that("safe_nlminb wraps optimizer errors", {
  expect_error(
    safe_nlminb(1, function(p) stop("kaboom"), context = "test model"),
    "Optimization failed for test model"
  )
})

test_that("safe_solve returns the inverse or NULL with a warning", {
  M <- matrix(c(2, 0, 0, 4), 2, 2)
  expect_equal(safe_solve(M), solve(M))
  singular <- matrix(c(1, 1, 1, 1), 2, 2)
  expect_warning(inv <- safe_solve(singular, "test matrix"),
                 "Could not invert test matrix")
  expect_null(inv)
})

test_that("safe_chol warns and repairs a non-positive-definite covariance", {
  # A rank-deficient (PSD but singular) matrix: chol() fails.
  Sigma <- matrix(c(1, 1, 1, 1), 2, 2)
  expect_warning(draws <- rmvnorm_chol(50, Sigma),
                 "not positive definite")
  expect_equal(dim(draws), c(50, 2))
  expect_true(all(is.finite(draws)))

  # A genuinely non-positive-definite matrix (a negative eigenvalue).
  Sigma2 <- matrix(c(1, 2, 2, 1), 2, 2)
  expect_warning(L <- safe_chol(Sigma2), "not positive definite")
  expect_true(all(is.finite(L)))
})

test_that("rmvnorm_chol handles a negative scalar variance", {
  expect_warning(draws <- rmvnorm_chol(10, matrix(-1, 1, 1)),
                 "variance is negative")
  expect_equal(dim(draws), c(10, 1))
  expect_true(all(draws == 0))
})

test_that("coef()/predict() warn on a non-converged fit but not a good one", {
  set.seed(101)
  n <- 400
  g <- factor(rep(1:20, each = 20))
  x <- rnorm(n)
  u <- rnorm(20, 0, 0.8)
  eta <- -0.3 + 0.5 * x + u[as.integer(g)]
  yb <- rbinom(n, 1, plogis(eta))
  d <- data.frame(yb = yb, x = x, g = g)

  # Force non-convergence with a one-iteration cap.
  bad <- suppressWarnings(
    gllamm(yb ~ x + (1 | g), data = d, family = stats::binomial(),
           control = list(iter.max = 1, eval.max = 1))
  )
  expect_false(isTRUE(bad$convergence$converged))
  expect_warning(coef(bad), "did not converge")

  good <- gllamm(yb ~ x + (1 | g), data = d, family = stats::binomial())
  expect_true(isTRUE(good$convergence$converged))
  expect_no_warning(coef(good))
})
