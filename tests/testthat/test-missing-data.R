# Package-wide missing-data policy: formula-based fitters listwise-delete
# with a warning and stay internally aligned (an NA fit must equal the
# explicit complete-case fit); matrix-response latent variable models use
# all observed responses.

test_that("GLMM with NA matches the explicit complete-case fit", {
  set.seed(1)
  d <- data.frame(y = rnorm(300), x = rnorm(300),
                  g = factor(rep(1:30, 10)))
  d$x[5] <- NA; d$y[10] <- NA
  dc <- d[stats::complete.cases(d), ]

  expect_warning(f_na <- gllamm(y ~ x + (1 | g), data = d), "listwise")
  f_cc <- gllamm(y ~ x + (1 | g), data = dc)
  expect_equal(f_na$logLik, f_cc$logLik, tolerance = 1e-6)
  expect_equal(unname(coef(f_na)$fixed), unname(coef(f_cc)$fixed),
               tolerance = 1e-6)
})

test_that("ordinal with NA matches the complete-case fit", {
  set.seed(2)
  d <- data.frame(x = rnorm(300), g = factor(rep(1:30, 10)))
  d$r <- sample(1:4, 300, TRUE)
  d$r[3] <- NA; d$x[8] <- NA
  dc <- d[stats::complete.cases(d), ]

  expect_warning(f_na <- fit_ordinal(r ~ x + (1 | g), data = d),
                 "listwise")
  f_cc <- fit_ordinal(r ~ x + (1 | g), data = dc)
  expect_equal(f_na$logLik, f_cc$logLik, tolerance = 1e-6)
})

test_that("survival with NA covariate or event matches complete-case", {
  set.seed(3)
  ds <- data.frame(time = rexp(200), status = rbinom(200, 1, .7),
                   x = rnorm(200), g = factor(rep(1:20, 10)))
  ds$x[2] <- NA; ds$status[5] <- NA
  dsc <- ds[stats::complete.cases(ds), ]

  f_na <- suppressWarnings(
    fit_survival(Surv(time, status) ~ x + (1 | g), data = ds,
                 distribution = "exponential"))
  f_cc <- fit_survival(Surv(time, status) ~ x + (1 | g), data = dsc,
                       distribution = "exponential")
  expect_equal(f_na$logLik, f_cc$logLik, tolerance = 1e-6)
})

test_that("npml with NA matches complete-case", {
  set.seed(4)
  dn <- data.frame(yb = rbinom(300, 1, .5), x = rnorm(300),
                   g = factor(rep(1:30, 10)))
  dn$x[7] <- NA
  dnc <- dn[stats::complete.cases(dn), ]

  f_na <- suppressWarnings(
    fit_npml(yb ~ x + (1 | g), data = dn, k = 2, family = binomial()))
  f_cc <- fit_npml(yb ~ x + (1 | g), data = dnc, k = 2,
                   family = binomial())
  expect_equal(f_na$logLik, f_cc$logLik, tolerance = 1e-4)
})

test_that("weights are aligned with listwise-deleted rows", {
  set.seed(5)
  d <- data.frame(yy = rbinom(300, 1, .5), x = rnorm(300),
                  g = factor(rep(1:30, 10)))
  d$x[5] <- NA
  w <- runif(300, 0.5, 2)
  keep <- stats::complete.cases(d)

  f_w <- suppressWarnings(
    gllamm(yy ~ x + (1 | g), data = d, family = binomial(), weights = w))
  f_wc <- gllamm(yy ~ x + (1 | g), data = d[keep, ],
                 family = binomial(), weights = w[keep])
  expect_equal(f_w$logLik, f_wc$logLik, tolerance = 1e-6)
})

test_that("crossed-RE models survive NA rows", {
  set.seed(6)
  d <- data.frame(yy = rbinom(400, 1, .5), x = rnorm(400),
                  g1 = factor(rep(1:40, 10)),
                  g2 = factor(rep(1:10, 40)))
  d$x[c(3, 17)] <- NA
  dc <- d[stats::complete.cases(d), ]

  f_na <- suppressWarnings(
    gllamm(yy ~ x + (1 | g1) + (1 | g2), data = d, family = binomial()))
  f_cc <- gllamm(yy ~ x + (1 | g1) + (1 | g2), data = dc,
                 family = binomial())
  expect_equal(f_na$logLik, f_cc$logLik, tolerance = 1e-5)
})

test_that("mixed responses listwise-delete with a warning", {
  set.seed(7)
  dm <- data.frame(y1 = rnorm(200), y2 = rbinom(200, 1, .5),
                   x = rnorm(200), g = factor(rep(1:20, 10)))
  dm$y1[4] <- NA
  expect_warning(
    f <- fit_mixed(list(gaussian = y1 ~ x, binomial = y2 ~ x),
                   random = ~ 1 | g, data = dm),
    "listwise")
  expect_true(is.finite(f$logLik))
})

test_that("matrix-response models use observed responses under MAR", {
  set.seed(8)
  # IRT: parameters from NA-laced data close to complete-data truth
  np <- 600; ni <- 12
  theta <- rnorm(np)
  b <- seq(-1.5, 1.5, length.out = ni)
  p <- plogis(outer(theta, b, "-"))
  resp <- matrix(rbinom(np * ni, 1, p), np, ni)
  resp_na <- resp
  resp_na[sample(length(resp_na), 0.15 * length(resp_na))] <- NA

  f_full <- fit_irt(resp, model = "Rasch")
  f_na <- fit_irt(resp_na, model = "Rasch")
  expect_lt(max(abs(f_full$item_parameters$difficulty -
                      f_na$item_parameters$difficulty)), 0.25)

  # LCA with NA already covered in EM; sanity-check it runs and classifies
  Y <- sapply(1:5, function(j) rbinom(400, 1, c(0.2, 0.8)[1 + (theta[1:400] > 0)]))
  Y[sample(length(Y), 200)] <- NA
  f_lca <- fit_lca(Y, nclass = 2)
  expect_true(f_lca$convergence$converged)
})

test_that("level-specific weights with incomplete rows error clearly", {
  set.seed(9)
  d <- data.frame(y = rnorm(200), x = rnorm(200),
                  g = factor(rep(1:20, 10)),
                  w2 = rep(runif(20, 0.5, 2), each = 10))
  d$x[5] <- NA
  expect_error(
    suppressWarnings(
      gllamm(y ~ x + (1 | g), data = d, weights = list(level2 = d$w2))),
    "complete data")
})
