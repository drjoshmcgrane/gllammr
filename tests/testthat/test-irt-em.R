# MML-EM estimation path for IRT (method = "em").

simulate_2pl <- function(seed = 7, np = 800, ni = 20) {
  set.seed(seed)
  theta <- rnorm(np)
  a <- runif(ni, 0.6, 2.0); b <- rnorm(ni)
  resp <- matrix(rbinom(np * ni, 1,
                        plogis(outer(theta, b, "-") *
                                 matrix(a, np, ni, byrow = TRUE))), np, ni)
  list(resp = resp, a = a, b = b, theta = theta)
}

test_that("EM 2PL matches mirt to numerical agreement", {
  skip_if_not_installed("mirt")
  skip_on_cran()  # cross-package agreement is CI-only
  s <- simulate_2pl()
  fit <- fit_irt(s$resp, model = "2PL", method = "em", se = FALSE)
  ref <- mirt::mirt(as.data.frame(s$resp), 1, itemtype = "2PL", verbose = FALSE)
  co <- mirt::coef(ref, simplify = TRUE)$items
  b_ref <- -co[, "d"] / co[, "a1"]

  expect_gt(cor(fit$item_parameters$discrimination, co[, "a1"]), 0.9999)
  expect_gt(cor(fit$item_parameters$difficulty, b_ref), 0.9999)
  expect_equal(fit$logLik, mirt::extract.mirt(ref, "logLik"), tolerance = 1e-4)
})

test_that("EM Rasch agrees with the Laplace path and improves the logLik", {
  set.seed(123)
  np <- 600; ni <- 25
  resp <- matrix(rbinom(np * ni, 1,
                        plogis(outer(rnorm(np), rnorm(ni), "-"))), np, ni)
  fit_em <- fit_irt(resp, model = "Rasch", method = "em", se = FALSE)
  fit_lap <- fit_irt(resp, model = "Rasch", se = FALSE)

  expect_gt(cor(fit_em$item_parameters$difficulty,
                fit_lap$item_parameters$difficulty), 0.9999)
  expect_equal(fit_em$ability_sd, fit_lap$ability_sd, tolerance = 0.02)
  # Quadrature evaluates the marginal likelihood more exactly than Laplace
  expect_gte(fit_em$logLik, fit_lap$logLik - 1e-6)
})

test_that("EM GRM matches mirt; EM PCM agrees with Laplace", {
  skip_if_not_installed("mirt")
  skip_on_cran()  # cross-package agreement is CI-only
  set.seed(7)
  np <- 800
  theta <- rnorm(np)
  taus <- t(sapply(rnorm(12), function(b0) b0 + c(-1, 0, 1)))
  resp <- sapply(1:12, function(j) {
    cum <- sapply(1:3, function(k) plogis(theta - taus[j, k]))
    1L + rowSums(matrix(runif(np), np, 3) < cum)
  })

  fit_g <- fit_irt(resp, model = "GRM", method = "em", se = FALSE)
  ref <- mirt::mirt(as.data.frame(resp), 1, itemtype = "graded", verbose = FALSE)
  co <- mirt::coef(ref, simplify = TRUE)$items
  expect_gt(cor(fit_g$item_parameters$discrimination, co[, "a1"]), 0.999)
  expect_equal(fit_g$logLik, mirt::extract.mirt(ref, "logLik"), tolerance = 1e-3)

  fit_p <- fit_irt(resp, model = "PCM", method = "em", se = FALSE)
  fit_pl <- fit_irt(resp, model = "PCM", se = FALSE)
  expect_equal(fit_p$ability_sd, fit_pl$ability_sd, tolerance = 0.05)
  expect_gte(fit_p$logLik, fit_pl$logLik - 1e-6)
})

test_that("EM handles the short-test 2PL case where Laplace diverges", {
  skip_if_not_installed("ltm")
  skip_on_cran()  # cross-package agreement is CI-only
  data("LSAT", package = "ltm", envir = environment())
  fit <- fit_irt(as.matrix(LSAT), model = "2PL", method = "em", se = FALSE)
  ref_coef <- coef(ltm::ltm(LSAT ~ z1))

  expect_equal(unname(fit$item_parameters$difficulty),
               unname(ref_coef[, "Dffclt"]), tolerance = 0.05)
  expect_equal(unname(fit$item_parameters$discrimination),
               unname(ref_coef[, "Dscrmn"]), tolerance = 0.05)
  expect_true(all(fit$item_parameters$discrimination < 2))  # no divergence
})

test_that("EM respects person weights and missing data", {
  s <- simulate_2pl(np = 300, ni = 10)
  resp <- s$resp
  resp[sample(length(resp), 200)] <- NA
  fit <- fit_irt(resp, model = "2PL", method = "em", se = FALSE)
  expect_true(fit$convergence$converged)

  # Integer weights == duplication
  w <- rep(c(1, 2), length.out = 300)
  fit_w <- fit_irt(s$resp, model = "Rasch", method = "em", weights = w, se = FALSE)
  fit_d <- fit_irt(rbind(s$resp, s$resp[w == 2, ]), model = "Rasch",
                   method = "em", se = FALSE)
  expect_equal(fit_w$item_parameters$difficulty,
               fit_d$item_parameters$difficulty, tolerance = 1e-3)
})

test_that("EM rejects multi-level specifications clearly", {
  s <- simulate_2pl(np = 100, ni = 5)
  pd <- data.frame(g = factor(rep(1:5, each = 20)))
  expect_error(fit_irt(s$resp, model = "Rasch", method = "em",
                       person_data = pd, random = ~ (1 | g), se = FALSE),
               "single-level")
})


test_that("C++ EM matches mirt exactly on a large polytomous battery", {
  skip_if_not_installed("mirt")
  skip_on_cran()   # ~15s with the mirt reference fit
  set.seed(9)
  np <- 2000; ni <- 60
  theta <- rnorm(np)
  a <- runif(ni, 0.7, 1.8)
  taus <- t(sapply(rnorm(ni), function(b0) b0 + c(-1.5, -0.5, 0.5, 1.5)))
  resp <- sapply(1:ni, function(j) {
    cum <- sapply(1:4, function(k) plogis(a[j] * (theta - taus[j, k])))
    1L + rowSums(matrix(runif(np), np, 4) < cum)
  })

  fit <- fit_irt(resp, model = "GRM", se = FALSE)  # se = FALSE: C++ EM
  ref <- mirt::mirt(as.data.frame(resp), 1, itemtype = "graded",
                    verbose = FALSE)
  co <- mirt::coef(ref, simplify = TRUE)$items

  expect_gt(cor(fit$item_parameters$discrimination, co[, "a1"]), 0.99999)
  expect_equal(fit$logLik, mirt::extract.mirt(ref, "logLik"), tolerance = 1e-6)
})
