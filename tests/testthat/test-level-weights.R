# Level-specific survey weights (E6). Integer frequency-weight equivalence:
# level-1 weights are exactly equivalent to duplicating rows; level-2
# integer weights are implemented by exact group replication under Laplace
# (identical to duplicating whole groups), and by weighting each group's
# log marginal likelihood under aghq (exact for arbitrary weights).

simulate_weighted <- function(seed = 81, g = 20, n_per = 15) {
  set.seed(seed)
  n <- g * n_per
  grp <- factor(rep(1:g, each = n_per))
  x <- rnorm(n)
  u <- rnorm(g, 0, 0.8)
  d <- data.frame(y = 1 + 0.5 * x + u[as.integer(grp)] + rnorm(n),
                  x = x, grp = grp)
  d$yb <- rbinom(n, 1, plogis(0.5 * x + u[as.integer(grp)]))
  d
}

test_that("level-1 integer weights equal row duplication exactly", {
  d <- simulate_weighted()
  n <- nrow(d)
  w1 <- ifelse(seq_len(n) %% 4 == 0, 2, 1)

  fit_w <- gllamm(yb ~ x + (1 | grp), data = d, family = stats::binomial(),
                  weights = list(level1 = w1))
  fit_dup <- gllamm(yb ~ x + (1 | grp),
                    data = rbind(d, d[seq_len(n) %% 4 == 0, ]),
                    family = stats::binomial())

  expect_equal(unname(coef(fit_w)$fixed), unname(coef(fit_dup)$fixed),
               tolerance = 1e-4)
  expect_equal(fit_w$logLik, fit_dup$logLik, tolerance = 1e-3)
})

test_that("level-2 integer weights match group duplication", {
  d <- simulate_weighted()
  g <- nlevels(d$grp)
  w2 <- ifelse(seq_len(g) <= 5, 2, 1)

  fit_w <- gllamm(y ~ x + (1 | grp), data = d, weights = list(level2 = w2))

  dup_rows <- d[as.integer(d$grp) <= 5, ]
  dup_rows$grp <- factor(paste0("dup", dup_rows$grp))
  d_dup <- rbind(d, dup_rows)
  d_dup$grp <- factor(as.character(d_dup$grp))
  fit_dup <- gllamm(y ~ x + (1 | grp), data = d_dup)

  # Replication-based weighting makes this identity exact
  expect_equal(unname(coef(fit_w)$fixed), unname(coef(fit_dup)$fixed),
               tolerance = 1e-5)
  expect_equal(sqrt(fit_w$coefficients$random_var[[1]][1, 1]),
               sqrt(fit_dup$coefficients$random_var[[1]][1, 1]),
               tolerance = 1e-4)
  expect_equal(fit_w$logLik, fit_dup$logLik, tolerance = 1e-5)
})

test_that("non-integer level-2 weights are rejected under Laplace", {
  d <- simulate_weighted(g = 10, n_per = 10)
  w2 <- runif(10, 0.5, 2)
  expect_error(gllamm(y ~ x + (1 | grp), data = d,
                      weights = list(level2 = w2)),
               "aghq")
})

test_that("aghq level-2 weights equal group duplication exactly", {
  d <- simulate_weighted(g = 12, n_per = 10)
  g <- nlevels(d$grp)
  w2 <- ifelse(seq_len(g) <= 4, 2, 1)
  fit_w <- gllamm(yb ~ x + (1 | grp), data = d, family = stats::binomial(),
                  weights = list(level2 = w2), integration = aghq(15))
  dup_rows <- d[as.integer(d$grp) <= 4, ]
  dup_rows$grp <- factor(paste0("dup", dup_rows$grp))
  d_dup <- rbind(d, dup_rows)
  d_dup$grp <- factor(as.character(d_dup$grp))
  fit_dup <- gllamm(yb ~ x + (1 | grp), data = d_dup,
                    family = stats::binomial(), integration = aghq(15))
  expect_equal(fit_w$logLik, fit_dup$logLik, tolerance = 1e-4)
  expect_equal(unname(coef(fit_w)$fixed), unname(coef(fit_dup)$fixed),
               tolerance = 1e-4)
})

test_that("binomial fits route list weights to the general engine", {
  d <- simulate_weighted(g = 12, n_per = 10)
  g <- nlevels(d$grp)
  w2 <- ifelse(seq_len(g) <= 4, 2, 1)
  fit_w <- gllamm(yb ~ x + (1 | grp), data = d, family = stats::binomial(),
                  weights = list(level2 = w2))
  dup_rows <- d[as.integer(d$grp) <= 4, ]
  dup_rows$grp <- factor(paste0("dup", dup_rows$grp))
  d_dup <- rbind(d, dup_rows)
  d_dup$grp <- factor(as.character(d_dup$grp))
  fit_dup <- gllamm(yb ~ x + (1 | grp), data = d_dup,
                    family = stats::binomial())
  expect_equal(fit_w$logLik, fit_dup$logLik, tolerance = 1e-5)
})

test_that("level-2 weights accept per-observation constant-within-group form", {
  d <- simulate_weighted(g = 10, n_per = 10)
  w2_obs <- ifelse(as.integer(d$grp) <= 3, 2, 1)   # per observation
  fit <- gllamm(y ~ x + (1 | grp), data = d, weights = list(level2 = w2_obs))
  expect_true(fit$convergence$converged)
})

test_that("weights validation errors clearly", {
  d <- simulate_weighted(g = 5, n_per = 8)
  expect_error(gllamm(y ~ x + (1 | grp), data = d,
                      weights = list(wrong = rep(1, 40))), "level1")
  w_bad <- ifelse(seq_len(40) %% 2 == 0, 2, 1)   # varies within group
  expect_error(gllamm(y ~ x + (1 | grp), data = d,
                      weights = list(level2 = w_bad)), "constant within")
})
