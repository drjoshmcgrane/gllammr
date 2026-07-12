# E11: person-level random slopes in multilevel IRT.

test_that("multilevel IRT recovers class-level intercept and slope variance", {
  set.seed(141)
  n_class <- 30; per_class <- 30; np <- n_class * per_class; ni <- 15
  class_id <- rep(1:n_class, each = per_class)
  ses <- rnorm(np)
  u0 <- rnorm(n_class, 0, 0.6)
  u1 <- rnorm(n_class, 0, 0.4)
  theta <- rnorm(np) + u0[class_id] + u1[class_id] * ses
  b <- seq(-1.5, 1.5, length.out = ni)
  resp <- matrix(rbinom(np * ni, 1, plogis(outer(theta, b, "-"))), np, ni)
  pdata <- data.frame(class = factor(class_id), ses = ses)

  fit <- fit_irt(resp, model = "Rasch", person_data = pdata,
                 random = ~ (ses | class), se = FALSE)

  expect_true(fit$convergence$converged)
  expect_setequal(fit$random_effects$group_names, c("class", "class:ses"))
  s <- fit$random_effects$sigma_random
  expect_equal(unname(s["class"]), 0.6, tolerance = 0.2)
  expect_equal(unname(s["class:ses"]), 0.4, tolerance = 0.2)

  # The slope model must dominate the intercept-only fit
  fit0 <- fit_irt(resp, model = "Rasch", person_data = pdata,
                  random = ~ (1 | class), se = FALSE)
  expect_gt(fit$logLik, fit0$logLik)
})

test_that("slope-only specification (0 + x | g) works", {
  skip_on_cran()  # random-slope IRT fit; random-slope smoke kept above
  set.seed(142)
  n_class <- 20; per_class <- 25; np <- n_class * per_class; ni <- 10
  class_id <- rep(1:n_class, each = per_class)
  ses <- rnorm(np)
  u1 <- rnorm(n_class, 0, 0.5)
  theta <- rnorm(np) + u1[class_id] * ses
  resp <- matrix(rbinom(np * ni, 1,
                        plogis(outer(theta, rnorm(ni), "-"))), np, ni)
  pdata <- data.frame(class = factor(class_id), ses = ses)

  fit <- fit_irt(resp, model = "Rasch", person_data = pdata,
                 random = ~ (0 + ses | class), se = FALSE)
  expect_true(fit$convergence$converged)
  expect_equal(fit$random_effects$group_names, "class:ses")
})

test_that("unknown slope covariate errors clearly", {
  resp <- matrix(rbinom(200, 1, 0.5), 20, 10)
  pdata <- data.frame(class = factor(rep(1:4, each = 5)))
  expect_error(fit_irt(resp, model = "Rasch", person_data = pdata,
                       random = ~ (nope | class), se = FALSE), "not found")
})
