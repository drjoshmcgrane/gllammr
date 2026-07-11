# Equivalence proof for the vectorized Monte Carlo marginal-prediction
# integrator. The reference implementations below are verbatim copies of the
# pre-vectorization per-sample Welford loops; every assertion checks that the
# new vectorized code reproduces them bit-for-bit (tolerance 1e-12), and that
# an end-to-end fixed-seed predict(type = "marginal") is unchanged (1e-10).

# ---- Reference: old per-sample Welford integrator (single term) ------------
old_mc_fixed_samples <- function(X, Z, beta, u_samples, inv_link_fn) {
  n_obs <- nrow(X)
  n_sim <- nrow(u_samples)
  eta_fixed <- as.vector(X %*% beta)
  mean_pred <- numeric(n_obs)
  m2_pred <- numeric(n_obs)
  for (s in 1:n_sim) {
    u <- u_samples[s, ]
    eta_random <- if (is.matrix(Z)) as.vector(Z %*% u) else Z * u
    pred_s <- inv_link_fn(eta_fixed + eta_random)
    delta <- pred_s - mean_pred
    mean_pred <- mean_pred + delta / s
    m2_pred <- m2_pred + delta * (pred_s - mean_pred)
  }
  list(fit = mean_pred, se = sqrt(m2_pred / (n_sim - 1)))
}


test_that("vectorized mc_integrate_fixed_samples matches old Welford loop", {
  inv <- plogis

  # (a) Z matrix, q = 2, small n_obs
  set.seed(1)
  n_obs <- 30; p <- 2; q <- 2; n_sim <- 500
  X <- matrix(rnorm(n_obs * p), n_obs, p)
  Z <- matrix(rnorm(n_obs * q), n_obs, q)
  beta <- c(0.4, -0.3)
  U <- matrix(rnorm(n_sim * q), n_sim, q)
  expect_equal(mc_integrate_fixed_samples(X, Z, beta, U, inv),
               old_mc_fixed_samples(X, Z, beta, U, inv), tolerance = 1e-12)

  # (b) Z as a vector (single random effect)
  set.seed(2)
  n_obs <- 50; n_sim <- 400
  X <- matrix(rnorm(n_obs * 2), n_obs, 2)
  Zv <- rnorm(n_obs)
  U1 <- matrix(rnorm(n_sim), n_sim, 1)
  expect_equal(mc_integrate_fixed_samples(X, Zv, beta, U1, inv),
               old_mc_fixed_samples(X, Zv, beta, U1, inv), tolerance = 1e-12)

  # (c) Larger n_obs, q = 2, Poisson-style exp link
  set.seed(3)
  n_obs <- 2000; q <- 2; n_sim <- 200
  X <- matrix(rnorm(n_obs * 2), n_obs, 2)
  Z <- matrix(rnorm(n_obs * q), n_obs, q)
  U <- matrix(rnorm(n_sim * q) * 0.2, n_sim, q)
  expect_equal(mc_integrate_fixed_samples(X, Z, beta, U, exp),
               old_mc_fixed_samples(X, Z, beta, U, exp), tolerance = 1e-12)
})


test_that("chunked reduction matches non-chunked and the old Welford loop", {
  inv <- plogis
  set.seed(4)
  n_obs <- 500; q <- 2; n_sim <- 300
  X <- matrix(rnorm(n_obs * 2), n_obs, 2)
  Z <- matrix(rnorm(n_obs * q), n_obs, q)
  beta <- c(0.2, -0.1)
  U <- matrix(rnorm(n_sim * q), n_sim, q)
  eta_fixed <- as.vector(X %*% beta)
  terms <- list(list(Z = Z, U = U))

  full  <- .mc_integrate_columns(eta_fixed, terms, inv, n_sim)
  # Force the memory-guard chunk path with a tiny threshold (150000 > 1e4)
  chunk <- .mc_integrate_columns(eta_fixed, terms, inv, n_sim,
                                 chunk_threshold = 1e4, chunk_cols = 50L)

  expect_equal(chunk$mean, full$mean, tolerance = 1e-12)
  expect_equal(chunk$m2,   full$m2,   tolerance = 1e-12)

  # Chunked moments must reproduce the old Welford se as well
  old <- old_mc_fixed_samples(X, Z, beta, U, inv)
  expect_equal(chunk$mean, old$fit, tolerance = 1e-12)
  expect_equal(sqrt(chunk$m2 / (n_sim - 1)), old$se, tolerance = 1e-12)
})


# ---- Reference: old per-replicate loop inside predict_marginal_gllamm -------
old_predict_marginal <- function(fit, n_sim, seed) {
  parts <- .gllamm_re_parts(fit, fit$data)
  X <- if (!is.null(fit$X)) fit$X else parts$md$X
  beta <- fit$coefficients$fixed
  inv_link <- get_inverse_link(fit$family)
  eta_fixed <- as.numeric(as.matrix(X) %*% beta)
  n_obs <- length(eta_fixed)
  set.seed(seed)
  mean_acc <- numeric(n_obs)
  m2_acc <- numeric(n_obs)
  for (sims in seq_len(n_sim)) {
    eta <- eta_fixed
    for (t in seq_len(parts$n_terms)) {
      u_t <- rmvnorm_chol(1, parts$Sigmas[[t]])
      eta <- eta + as.vector(parts$md$Z[[t]] %*% as.numeric(u_t))
    }
    mu <- inv_link(eta)
    delta <- mu - mean_acc
    mean_acc <- mean_acc + delta / sims
    m2_acc <- m2_acc + delta * (mu - mean_acc)
  }
  list(fit = mean_acc, se = sqrt(m2_acc / (n_sim * pmax(n_sim - 1, 1))))
}


test_that("end-to-end predict(type='marginal') reproduces the old loop", {
  skip_if_not_installed("Matrix")

  # Random-intercept binomial fit (single term, q = 1)
  set.seed(42)
  g <- 40; n <- 400
  grp <- factor(rep(1:g, length.out = n))
  x <- rnorm(n); u <- rnorm(g)
  yb <- rbinom(n, 1, plogis(0.3 + 0.5 * x + u[as.integer(grp)]))
  d <- data.frame(yb = yb, x = x, grp = grp)
  fit <- gllamm(yb ~ x + (1 | grp), data = d, family = stats::binomial())

  seed <- 2024L
  set.seed(seed)
  new <- predict(fit, type = "marginal", n_sim = 1000, se.fit = TRUE)
  old <- old_predict_marginal(fit, 1000, seed)
  expect_equal(unname(new$fit),    unname(old$fit), tolerance = 1e-10)
  expect_equal(unname(new$se.fit), unname(old$se),  tolerance = 1e-10)

  # Random-slope binomial fit (single term, q = 2) exercises the q>1 draw order
  set.seed(43)
  g2 <- 30; n2 <- 600
  grp2 <- factor(rep(1:g2, length.out = n2)); x2 <- rnorm(n2)
  U <- matrix(rnorm(g2 * 2), g2) %*% chol(matrix(c(1, .3, .3, .5), 2))
  eta2 <- 0.2 + U[as.integer(grp2), 1] + (0.5 + U[as.integer(grp2), 2]) * x2
  yb2 <- rbinom(n2, 1, plogis(eta2))
  d2 <- data.frame(y = yb2, x = x2, grp = grp2)
  fit2 <- gllamm(y ~ x + (x | grp), data = d2, family = stats::binomial())

  set.seed(seed)
  new2 <- predict(fit2, type = "marginal", n_sim = 800, se.fit = TRUE)
  old2 <- old_predict_marginal(fit2, 800, seed)
  expect_equal(unname(new2$fit),    unname(old2$fit), tolerance = 1e-10)
  expect_equal(unname(new2$se.fit), unname(old2$se),  tolerance = 1e-10)
})
