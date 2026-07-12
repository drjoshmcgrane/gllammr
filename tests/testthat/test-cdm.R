# Cognitive diagnosis models (DINA / DINO / G-DINA)

sim_dina <- function(n = 800, seed = 1) {
  set.seed(seed)
  Q <- rbind(diag(3), diag(3), c(1, 1, 0), c(0, 1, 1), c(1, 0, 1),
             c(1, 1, 1))
  J <- nrow(Q)
  alpha <- matrix(rbinom(n * 3, 1, 0.55), n, 3)
  g <- runif(J, 0.05, 0.2)
  s <- runif(J, 0.05, 0.2)
  eta <- sapply(seq_len(J), function(j) {
    as.integer(rowSums(alpha[, Q[j, ] == 1, drop = FALSE]) == sum(Q[j, ]))
  })
  Y <- sapply(seq_len(J), function(j) {
    rbinom(n, 1, ifelse(eta[, j] == 1, 1 - s[j], g[j]))
  })
  list(Y = Y, Q = Q, alpha = alpha, g = g, s = s)
}

test_that("DINA recovers guess/slip and attribute profiles", {
  d <- sim_dina(n = 1000, seed = 101)
  fit <- fit_cdm(d$Y, d$Q, model = "dina")

  expect_true(fit$convergence$converged)
  ghat <- sapply(fit$item_params, function(e) e$guess)
  shat <- sapply(fit$item_params, function(e) e$slip)
  expect_lt(max(abs(ghat - d$g)), 0.07)
  expect_lt(max(abs(shat - d$s)), 0.07)
  expect_gt(mean(fit$modal_attributes == d$alpha), 0.9)

  # Posterior structure
  expect_equal(dim(fit$posterior), c(1000, 8))
  expect_equal(unname(rowSums(fit$posterior)), rep(1, 1000), tolerance = 1e-8)
  expect_true(all(fit$attribute_posteriors >= 0 &
                    fit$attribute_posteriors <= 1))
  # Monotonicity: guess <= 1 - slip
  expect_true(all(ghat <= 1 - shat + 1e-8))
})

test_that("DINA matches CDM::din", {
  skip_if_not_installed("CDM")
  skip_on_cran()  # cross-package agreement is CI-only
  d <- sim_dina(n = 1000, seed = 101)
  fit <- fit_cdm(d$Y, d$Q, model = "dina")
  ref <- CDM::din(d$Y, q.matrix = d$Q, rule = "DINA", progress = FALSE)

  expect_equal(fit$logLik, ref$loglike, tolerance = 0.05)
  ghat <- sapply(fit$item_params, function(e) e$guess)
  shat <- sapply(fit$item_params, function(e) e$slip)
  expect_lt(mean(abs(ghat - ref$guess$est)), 0.01)
  expect_lt(mean(abs(shat - ref$slip$est)), 0.01)
})

test_that("saturated G-DINA matches CDM::gdina and nests DINA", {
  skip_if_not_installed("CDM")
  skip_on_cran()  # cross-package agreement is CI-only
  set.seed(202)
  Q <- rbind(diag(3), diag(3), c(1, 1, 0), c(0, 1, 1), c(1, 0, 1),
             c(1, 1, 1))
  J <- nrow(Q); n <- 1200
  alpha <- matrix(rbinom(n * 3, 1, 0.5), n, 3)
  Y <- sapply(seq_len(J), function(j) {
    meas <- which(Q[j, ] == 1)
    lev <- rowSums(alpha[, meas, drop = FALSE])
    base <- runif(1, 0.08, 0.2); top <- runif(1, 0.8, 0.92)
    rbinom(n, 1, base + (top - base) * lev / length(meas))
  })

  fit <- fit_cdm(Y, Q, model = "gdina")
  ref <- CDM::gdina(Y, q.matrix = Q, progress = FALSE)
  expect_equal(fit$logLik, ref$loglike, tolerance = 0.05)

  # G-DINA nests DINA: saturated logLik >= DINA logLik on the same data
  fit_d <- fit_cdm(Y, Q, model = "dina")
  expect_gte(fit$logLik, fit_d$logLik - 1e-6)
})

test_that("DINO handles disjunctive items", {
  set.seed(33)
  Q <- rbind(diag(2), diag(2), c(1, 1), c(1, 1))
  J <- nrow(Q); n <- 1000
  alpha <- matrix(rbinom(n * 2, 1, 0.5), n, 2)
  eta <- sapply(seq_len(J), function(j) {
    as.integer(rowSums(alpha[, Q[j, ] == 1, drop = FALSE]) > 0)
  })
  Y <- sapply(seq_len(J), function(j) {
    rbinom(n, 1, ifelse(eta[, j] == 1, 0.85, 0.15))
  })
  fit <- fit_cdm(Y, Q, model = "dino")
  ghat <- sapply(fit$item_params, function(e) e$guess)
  shat <- sapply(fit$item_params, function(e) e$slip)
  expect_lt(max(abs(ghat - 0.15)), 0.08)
  expect_lt(max(abs(shat - 0.15)), 0.08)
})

test_that("monotonicity constrains and unconstrained nests it", {
  set.seed(44)
  d <- sim_dina(n = 400, seed = 44)
  fit_m <- fit_cdm(d$Y, d$Q, model = "gdina", monotone = TRUE)
  fit_u <- fit_cdm(d$Y, d$Q, model = "gdina", monotone = FALSE)
  expect_gte(fit_u$logLik, fit_m$logLik - 1e-6)
  # Monotone fit satisfies the lattice order for every item
  for (j in seq_along(fit_m$item_params)) {
    pr <- fit_m$item_params[[j]]$prob
    pats <- do.call(rbind, lapply(strsplit(names(pr), ""), as.integer))
    for (a in seq_along(pr)) {
      for (b in seq_along(pr)) {
        if (a != b && all(pats[a, ] <= pats[b, ])) {
          expect_lte(pr[a], pr[b] + 1e-8)
        }
      }
    }
  }
})

test_that("attribute hierarchy prunes the profile space", {
  set.seed(55)
  # Linear hierarchy A1 -> A2 -> A3: only 4 admissible profiles
  Q <- rbind(diag(3), diag(3), c(1, 1, 0), c(0, 1, 1))
  J <- nrow(Q); n <- 900
  # Generate respecting the hierarchy: mastery level 0-3
  lev <- sample(0:3, n, TRUE)
  alpha <- cbind(lev >= 1, lev >= 2, lev >= 3) * 1L
  eta <- sapply(seq_len(J), function(j) {
    as.integer(rowSums(alpha[, Q[j, ] == 1, drop = FALSE]) == sum(Q[j, ]))
  })
  Y <- sapply(seq_len(J), function(j) {
    rbinom(n, 1, ifelse(eta[, j] == 1, 0.9, 0.12))
  })
  fit <- fit_cdm(Y, Q, model = "dina",
                 hierarchy = list(c(1, 2), c(2, 3)))
  expect_equal(fit$nclass, 4)
  expect_setequal(fit$profile_labels, c("000", "100", "110", "111"))
  expect_gt(mean(fit$modal_attributes == alpha), 0.9)
})

test_that("fit_cdm validates inputs", {
  set.seed(66)
  Y <- matrix(rbinom(200, 1, 0.5), 50, 4)
  Q <- rbind(c(1, 0), c(0, 1), c(1, 1), c(0, 0))
  expect_error(fit_cdm(Y, Q), "all zero")
  expect_error(fit_cdm(Y, Q[1:3, ]), "one row per item")
  expect_error(fit_cdm(Y + 1, rbind(Q[1:3, ], c(1, 0))), "binary")
  Qok <- rbind(Q[1:3, ], c(1, 0))
  expect_error(fit_cdm(Y, Qok, hierarchy = list(c(1, 2), c(2, 1))), "cycle")
  expect_error(fit_cdm(Y, Qok, weights = rep(1, 10)), "weights length")
})

test_that("fit_cdm handles missing data and weights", {
  d <- sim_dina(n = 600, seed = 77)
  Yna <- d$Y
  Yna[sample(length(Yna), 300)] <- NA
  fit <- fit_cdm(Yna, d$Q, model = "dina")
  expect_true(fit$convergence$converged)
  ghat <- sapply(fit$item_params, function(e) e$guess)
  expect_lt(mean(abs(ghat - d$g)), 0.06)

  w <- runif(600, 0.5, 2)
  fit_w <- fit_cdm(d$Y, d$Q, model = "dina", weights = w)
  expect_true(fit_w$convergence$converged)
})
