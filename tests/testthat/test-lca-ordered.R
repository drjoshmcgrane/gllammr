# Order-restricted (Croon) latent class models

test_that("pava_weighted solves weighted isotonic regression", {
  # Already monotone: unchanged
  expect_equal(.pava_weighted(c(1, 2, 3), c(1, 1, 1)), c(1, 2, 3))
  # Single violation pools to the weighted mean
  expect_equal(.pava_weighted(c(2, 1), c(1, 1)), c(1.5, 1.5))
  expect_equal(.pava_weighted(c(2, 1), c(3, 1)), c(1.75, 1.75))
  # Cascading violations
  out <- .pava_weighted(c(3, 2, 1), c(1, 1, 1))
  expect_equal(out, c(2, 2, 2))
  # Result is always nondecreasing
  set.seed(1)
  for (i in 1:20) {
    y <- rnorm(6); w <- runif(6, 0.1, 2)
    expect_true(all(diff(.pava_weighted(y, w)) >= -1e-12))
  }
})

test_that("ordered LCA matches unrestricted fit when classes are separated", {
  set.seed(42)
  n <- 600
  cls <- sample(1:2, n, TRUE, prob = c(0.55, 0.45))
  pmat <- rbind(c(0.15, 0.85), c(0.2, 0.8), c(0.1, 0.9), c(0.25, 0.75))
  Y <- sapply(1:4, function(j) rbinom(n, 1, pmat[j, cls]))

  fit_u <- fit_lca(Y, nclass = 2, control = list(n_starts = 5))
  fit_o <- fit_lca(Y, nclass = 2, ordering = "increasing",
                   control = list(n_starts = 5))

  # The unrestricted optimum is interior to the order constraint here,
  # so both reach the same maximum
  expect_equal(fit_o$logLik, fit_u$logLik, tolerance = 1e-5)
  # Ordered fit pins the labels: probabilities nondecreasing across classes
  expect_true(all(apply(fit_o$item_probs, 1, diff) >= -1e-8))
})

test_that("ordering constrains every binary item across 3 classes", {
  set.seed(7)
  n <- 900
  cls <- sample(1:3, n, TRUE, prob = c(0.4, 0.35, 0.25))
  # 8 items keep the 3-class model well-identified (EM converges briskly;
  # with ~5 items 3-class binary LCA is near-flat regardless of ordering)
  pmat <- rbind(c(0.1, 0.5, 0.9), c(0.2, 0.4, 0.8),
                c(0.15, 0.55, 0.85), c(0.3, 0.5, 0.7),
                c(0.1, 0.3, 0.6), c(0.05, 0.45, 0.9),
                c(0.25, 0.6, 0.95), c(0.1, 0.4, 0.75))
  Y <- sapply(1:8, function(j) rbinom(n, 1, pmat[j, cls]))

  fit <- fit_lca(Y, nclass = 3, ordering = "increasing",
                 control = list(n_starts = 5))
  expect_true(fit$convergence$converged)
  for (j in 1:8) {
    expect_true(all(diff(fit$item_probs[j, ]) >= -1e-8))
  }
  # Recovers the generating probabilities
  expect_lt(max(abs(fit$item_probs - pmat)), 0.15)
})

test_that("ordering binds when items disagree about class order", {
  set.seed(11)
  n <- 800
  cls <- sample(1:2, n, TRUE)
  # Item 4 is anti-monotone relative to items 1-3: the constraint binds
  pmat <- rbind(c(0.2, 0.8), c(0.25, 0.75), c(0.2, 0.8), c(0.8, 0.2))
  Y <- sapply(1:4, function(j) rbinom(n, 1, pmat[j, cls]))

  fit_u <- fit_lca(Y, nclass = 2, control = list(n_starts = 5))
  fit_o <- fit_lca(Y, nclass = 2, ordering = "increasing",
                   control = list(n_starts = 5))

  expect_lt(fit_o$logLik, fit_u$logLik - 1)
  expect_true(all(apply(fit_o$item_probs, 1, diff) >= -1e-8))
})

test_that("ordered LCA constrains gaussian indicator means", {
  set.seed(23)
  n <- 700
  cls <- sample(1:2, n, TRUE)
  Y <- cbind(rbinom(n, 1, c(0.2, 0.8)[cls]),
             rbinom(n, 1, c(0.3, 0.7)[cls]),
             rnorm(n, c(0, 1.5)[cls], 1))

  fit <- fit_lca(Y, nclass = 2, ordering = "increasing",
                 control = list(n_starts = 5))
  expect_true(all(diff(fit$gaussian_params$means[1, ]) >= -1e-8))
  expect_true(all(apply(fit$item_probs[1:2, ], 1, diff) >= -1e-8))
  # Means recovered on the right scale
  expect_equal(unname(fit$gaussian_params$means[1, ]), c(0, 1.5),
               tolerance = 0.25)
})

test_that("ordering validates its requirements", {
  set.seed(31)
  Y <- cbind(rbinom(100, 1, 0.5), sample(1:3, 100, TRUE))
  expect_error(fit_lca(Y, nclass = 2, ordering = "increasing"),
               "categorical")
  Yb <- cbind(rbinom(100, 1, 0.5), rbinom(100, 1, 0.5))
  expect_error(fit_lca(Yb, nclass = 2, ordering = "increasing",
                       method = "tmb"),
               "method")
})

test_that("ordered LCA works with case weights", {
  set.seed(47)
  n <- 500
  cls <- sample(1:2, n, TRUE)
  Y <- sapply(1:4, function(j) rbinom(n, 1, c(0.2, 0.8)[cls]))
  w <- runif(n, 0.5, 2)

  fit <- fit_lca(Y, nclass = 2, weights = w, ordering = "increasing",
                 control = list(n_starts = 3))
  expect_true(fit$convergence$converged)
  expect_true(all(apply(fit$item_probs, 1, diff) >= -1e-8))
})
