# EIRT deep-audit regressions: saturated identities against descriptive
# IRT, multilevel cross-walks against glmer, exact frequency-weight
# semantics, and method correctness on multilevel fits.

make_dichot <- function(np = 300, ni = 8, seed = 42, school = FALSE) {
  set.seed(seed)
  theta <- rnorm(np)
  b <- seq(-1.2, 1.2, length.out = ni)
  eta_extra <- 0
  sch <- NULL
  if (school) {
    sch <- factor(rep(seq_len(np / 20), each = 20))
    eta_extra <- rnorm(nlevels(sch), 0, 0.7)[sch]
  }
  resp <- sapply(seq_len(ni), function(j)
    rbinom(np, 1, plogis(theta + eta_extra - b[j])))
  list(resp = resp, sch = sch)
}

test_that("saturated EIRT difficulty model equals descriptive IRT exactly", {
  d <- make_dichot()
  ni <- ncol(d$resp)
  fe <- fit_eirt(d$resp, data.frame(item = factor(seq_len(ni))),
                 difficulty_formula = ~ item, model = "Rasch",
                 item_residuals = FALSE)
  fi <- fit_irt(d$resp, model = "Rasch", method = "laplace", se = FALSE)
  expect_equal(as.numeric(logLik(fe)), as.numeric(logLik(fi)),
               tolerance = 1e-5)
})

test_that("EIRT person fweights reproduce duplicated-data fits exactly", {
  d <- make_dichot(np = 120, ni = 8, seed = 8)
  np <- nrow(d$resp)
  z <- rnorm(ncol(d$resp))
  w <- sample(1:3, np, replace = TRUE)
  idx <- rep(seq_len(np), w)

  fa <- fit_eirt(d$resp, data.frame(z = z), difficulty_formula = ~ z,
                 model = "Rasch", item_residuals = FALSE, weights = w)
  fb <- fit_eirt(d$resp[idx, ], data.frame(z = z),
                 difficulty_formula = ~ z, model = "Rasch",
                 item_residuals = FALSE)
  expect_equal(as.numeric(logLik(fa)), as.numeric(logLik(fb)),
               tolerance = 1e-6)
  expect_equal(fa$regression_coefficients$difficulty,
               fb$regression_coefficients$difficulty, tolerance = 1e-5)
  expect_equal(fa$ability_sd, fb$ability_sd, tolerance = 1e-5)

  # EM route weights the log marginal likelihood directly: also exact
  ea <- fit_irt(d$resp, model = "Rasch", weights = w, se = FALSE)
  eb <- fit_irt(d$resp[idx, ], model = "Rasch", se = FALSE)
  expect_equal(as.numeric(logLik(ea)), as.numeric(logLik(eb)),
               tolerance = 1e-4)

  # Laplace fits refuse non-integer person weights (scaling the joint
  # contribution would make the objective unbounded in sigma_theta)
  expect_error(
    fit_eirt(d$resp, data.frame(z = z), difficulty_formula = ~ z,
             model = "Rasch", weights = w + 0.5),
    "Non-integer")
  expect_error(
    fit_irt(d$resp, model = "Rasch", method = "laplace", weights = w + 0.5, se = FALSE),
    "Non-integer")
})

test_that("multilevel EIRT matches glmer on the crossed Rasch cross-walk", {
  skip_if_not_installed("lme4")
  d <- make_dichot(np = 200, ni = 8, seed = 99, school = TRUE)
  np <- nrow(d$resp); ni <- ncol(d$resp)
  fe <- fit_eirt(d$resp, data.frame(item = factor(seq_len(ni))),
                 difficulty_formula = ~ item, model = "Rasch",
                 item_residuals = FALSE,
                 person_data = data.frame(sch = d$sch),
                 random = ~ (1 | sch))
  long <- data.frame(y = as.vector(t(d$resp)),
                     item = factor(rep(seq_len(ni), np)),
                     id = factor(rep(seq_len(np), each = ni)),
                     sch = factor(rep(d$sch, each = ni)))
  fg <- suppressWarnings(
    lme4::glmer(y ~ 0 + item + (1 | id) + (1 | sch), data = long,
                family = binomial()))
  expect_equal(as.numeric(logLik(fe)), as.numeric(logLik(fg)),
               tolerance = 0.02)
  vc <- as.data.frame(lme4::VarCorr(fg))
  expect_equal(unname(fe$ability_sd), vc$sdcor[vc$grp == "id"],
               tolerance = 0.02)
  expect_equal(unname(fe$random_effects$sigma_random),
               vc$sdcor[vc$grp == "sch"], tolerance = 0.02)
})

test_that("EIRT difficulty standard errors match glmer", {
  skip_if_not_installed("lme4")
  d <- make_dichot(np = 300, ni = 8, seed = 7)
  np <- nrow(d$resp); ni <- ncol(d$resp)
  fe <- fit_eirt(d$resp, data.frame(item = factor(seq_len(ni))),
                 difficulty_formula = ~ item, model = "Rasch",
                 item_residuals = FALSE)
  sdr <- summary(fe$tmb_sdr)
  d_est <- sdr[rownames(sdr) == "difficulty", "Estimate"]
  d_se <- sdr[rownames(sdr) == "difficulty", "Std. Error"]
  long <- data.frame(y = as.vector(t(d$resp)),
                     item = factor(rep(seq_len(ni), np)),
                     id = factor(rep(seq_len(np), each = ni)))
  fg <- suppressWarnings(
    lme4::glmer(y ~ 0 + item + (1 | id), data = long,
                family = binomial()))
  expect_lt(max(abs(d_est - (-lme4::fixef(fg)))), 0.02)
  expect_lt(max(abs(d_se - sqrt(diag(as.matrix(vcov(fg)))))), 0.005)
})

test_that("multilevel EIRT simulate and marginal use the group REs", {
  d <- make_dichot(np = 200, ni = 10, seed = 99, school = TRUE)
  ni <- ncol(d$resp)
  z <- rnorm(ni)
  fm <- fit_eirt(d$resp, data.frame(z = z), difficulty_formula = ~ z,
                 model = "Rasch", item_residuals = TRUE,
                 person_data = data.frame(sch = d$sch),
                 random = ~ (1 | sch))

  # Simulation must reproduce the school-level variance, not just theta
  sm <- simulate(fm, nsim = 20, seed = 2)
  smean <- mean(sapply(sm, function(s)
    var(tapply(rowMeans(s), d$sch, mean))))
  emp <- var(tapply(rowMeans(d$resp), d$sch, mean))
  expect_gt(smean, 0.5 * emp)

  # Marginal predictions integrate over the TOTAL latent distribution
  mg <- predict(fm, type = "marginal")
  expect_lt(abs(mean(mg) - mean(d$resp)), 0.03)
})

test_that("fit_eirt rejects misuse instead of silently proceeding", {
  d <- make_dichot(np = 150, ni = 8, seed = 5)
  z <- rnorm(8)

  # Collinear item predictors -> singular Hessian; named, early error
  expect_error(
    fit_eirt(d$resp, data.frame(z = z, z2 = 2 * z),
             difficulty_formula = ~ z + z2, model = "Rasch"),
    "rank deficient")

  # threshold_formula has no GRM interpretation; must not be ignored
  set.seed(5)
  resp3 <- sapply(1:6, function(j)
    1L + rowSums(outer(rnorm(150) - (j - 3) / 2, c(-0.8, 0.8), ">")))
  expect_error(
    fit_eirt(resp3, data.frame(z = rnorm(6)), difficulty_formula = ~ 1,
             threshold_formula = ~ z, model = "GRM"),
    "only supported for PCM and GPCM")
})

test_that("multilevel fit_irt simulate and marginal use the group REs", {
  d <- make_dichot(np = 200, ni = 10, seed = 21, school = TRUE)
  fm <- fit_irt(d$resp, model = "Rasch",
                person_data = data.frame(sch = d$sch),
                random = ~ (1 | sch), se = FALSE)
  sm <- simulate(fm, nsim = 20, seed = 2)
  smean <- mean(sapply(sm, function(s)
    var(tapply(rowMeans(s), d$sch, mean))))
  emp <- var(tapply(rowMeans(d$resp), d$sch, mean))
  expect_gt(smean, 0.5 * emp)

  mg <- predict(fm, type = "marginal")
  expect_lt(abs(mean(mg) - mean(d$resp)), 0.03)
})

test_that("step-level predictors: recovery, identity, and identification", {
  set.seed(9)
  np <- 400; ni <- 10; K <- 4
  theta <- rnorm(np)
  z_item <- rnorm(ni)
  x_step <- matrix(rnorm(ni * (K - 1)), ni, K - 1)
  b <- 0.3 + 0.6 * z_item
  resp <- sapply(seq_len(ni), function(j) {
    delta <- b[j] + 0.7 * x_step[j, ] + c(-0.9, 0, 0.9)
    sapply(seq_len(np), function(p) {
      cs <- cumsum(c(0, theta[p] - delta))
      sample.int(K, 1, prob = exp(cs) / sum(exp(cs)))
    })
  })
  step_data <- data.frame(x = as.vector(t(x_step)))

  # Item-level and step-level predictors estimated together, with SEs
  f <- fit_eirt(resp, data.frame(z_item = z_item),
                difficulty_formula = ~ z_item,
                step_formula = ~ x, step_data = step_data,
                model = "PCM", item_residuals = TRUE)
  expect_true(f$convergence$converged)
  expect_equal(unname(f$regression_coefficients$difficulty["z_item"]), 0.6,
               tolerance = 0.15)
  eta <- f$regression_coefficients$step
  expect_equal(unname(eta["x", "estimate"]), 0.7, tolerance = 0.1)
  expect_true(is.finite(eta["x", "se"]) && eta["x", "se"] > 0)

  # A step covariate constant within items is EXACTLY an item covariate
  zc <- rnorm(ni)
  resp2 <- sapply(seq_len(ni), function(j) {
    delta <- 0.5 * zc[j] + c(-0.8, 0, 0.8)
    sapply(seq_len(np), function(p) {
      cs <- cumsum(c(0, theta[p] - delta))
      sample.int(K, 1, prob = exp(cs) / sum(exp(cs)))
    })
  })
  sd2 <- data.frame(zc = rep(zc, each = K - 1))
  fa <- fit_eirt(resp2, data.frame(one = rep(1, ni)),
                 difficulty_formula = ~ 0 + one,
                 step_formula = ~ zc, step_data = sd2,
                 model = "PCM", item_residuals = FALSE)
  fb <- fit_eirt(resp2, data.frame(zc = zc), difficulty_formula = ~ zc,
                 threshold_formula = ~ 1, model = "PCM",
                 item_residuals = FALSE)
  expect_equal(as.numeric(logLik(fa)), as.numeric(logLik(fb)),
               tolerance = 1e-5)

  # Cross-level collinearity is rejected, not silently ridged
  expect_error(
    fit_eirt(resp2, data.frame(zc = zc), difficulty_formula = ~ zc,
             step_formula = ~ zc, step_data = sd2, model = "PCM"),
    "rank deficient")
})

test_that("threshold regression is identified with shared covariates", {
  set.seed(9)
  np <- 400; ni <- 10; K <- 4
  theta <- rnorm(np); z <- rnorm(ni)
  resp <- sapply(seq_len(ni), function(j) {
    delta <- 0.3 + 0.5 * z[j] + c(-0.8, 0, 0.8) + c(-0.3, 0.1, 0.2) * z[j]
    sapply(seq_len(np), function(p) {
      cs <- cumsum(c(0, theta[p] - delta))
      sample.int(K, 1, prob = exp(cs) / sum(exp(cs)))
    })
  })
  # Same covariate at the item level AND with step-specific effects:
  # previously a flat ridge (all SEs NaN); xi rows now sum to zero so
  # the location effect and the deviations separate
  f <- fit_eirt(resp, data.frame(z = z), difficulty_formula = ~ z,
                threshold_formula = ~ z, model = "PCM",
                item_residuals = FALSE)
  sdr <- summary(f$tmb_sdr, select = "fixed")
  g <- sdr[rownames(sdr) == "gamma", ]
  expect_true(all(is.finite(g[, "Std. Error"])))
  expect_equal(unname(g[2, "Estimate"]), 0.5, tolerance = 0.1)
  xi <- f$regression_coefficients$threshold
  expect_equal(unname(rowSums(xi)), c(0, 0), tolerance = 1e-8)
})

test_that("step_formula guards reject unsupported models and missing step_data", {
  d <- make_dichot(np = 150, ni = 8, seed = 11)
  ni <- ncol(d$resp)
  sd_dummy <- data.frame(x = rnorm(ni))

  # Guard 1: step_formula has no interpretation outside PCM/GPCM (the
  # adjacent-categories / LPCM framework) - must not be silently ignored
  expect_error(
    fit_eirt(d$resp, data.frame(z = rnorm(ni)), difficulty_formula = ~ 1,
             step_formula = ~ x, step_data = sd_dummy, model = "Rasch"),
    "step_formula is only supported for PCM and GPCM")

  # Guard 2: step_formula requires step_data (one row per item-step
  # combination); omitting it must error rather than proceed with NULL
  expect_error(
    fit_eirt(d$resp, data.frame(z = rnorm(ni)), difficulty_formula = ~ 1,
             step_formula = ~ x, model = "PCM"),
    "step_data must be provided with step_formula")

  # Guard 3: step-level regression needs adjacent-categories steps, so
  # items with fewer than 3 categories (e.g. dichotomous) are rejected
  sd_full <- data.frame(x = rnorm(ni))
  expect_error(
    fit_eirt(d$resp, data.frame(z = rnorm(ni)), difficulty_formula = ~ 1,
             step_formula = ~ x, step_data = sd_full, model = "PCM"),
    "require items with at least")
})
