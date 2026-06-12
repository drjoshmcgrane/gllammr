# Cluster-robust (sandwich) standard errors (E9).

test_that("sandwich vcov approximates model vcov under correct specification", {
  set.seed(91)
  g <- 60; n_per <- 12; n <- g * n_per
  grp <- factor(rep(1:g, each = n_per))
  x <- rnorm(n)
  u <- rnorm(g, 0, 0.7)
  d <- data.frame(y = 1 + 0.5 * x + u[as.integer(grp)] + rnorm(n),
                  x = x, grp = grp)

  fit <- gllamm(y ~ x + (1 | grp), data = d)
  V_m <- vcov(fit)
  V_s <- vcov(fit, type = "sandwich")

  expect_equal(dim(V_s), dim(V_m))
  expect_true(all(eigen(V_s, only.values = TRUE)$values > 0))
  ratio <- sqrt(diag(V_s)) / sqrt(diag(V_m))
  expect_true(all(ratio > 0.7 & ratio < 1.5))
})

test_that("sandwich vcov works for binomial fits", {
  set.seed(92)
  g <- 50; n <- 600
  grp <- factor(rep(1:g, each = n %/% g))
  x <- rnorm(n)
  u <- rnorm(g, 0, 0.7)
  d <- data.frame(x = x, grp = grp,
                  yb = rbinom(n, 1, plogis(0.5 * x + u[as.integer(grp)])))

  fit <- gllamm(yb ~ x + (1 | grp), data = d, family = stats::binomial())
  V_s <- vcov(fit, type = "sandwich")
  ratio <- sqrt(diag(V_s)) / sqrt(diag(vcov(fit)))
  expect_true(all(ratio > 0.7 & ratio < 1.5))
})

test_that("sandwich requires a gllamm()-fitted model", {
  fake <- structure(list(formula = NULL), class = "gllamm")
  expect_error(vcov(fake, type = "sandwich"), "gllamm")
})
