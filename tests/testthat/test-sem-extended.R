# Extended SEM: correlated factors, MIMIC, SEs, fit indices, FIML

sim_two_factor <- function(n = 400, seed = 42, rho = 0.6) {
  set.seed(seed)
  f <- cbind(rnorm(n), rnorm(n))
  f[, 2] <- rho * f[, 1] + sqrt(1 - rho^2) * f[, 2]
  data.frame(
    x1 = f[, 1] + rnorm(n, 0, .6), x2 = 0.8 * f[, 1] + rnorm(n, 0, .6),
    x3 = 1.2 * f[, 1] + rnorm(n, 0, .6),
    y1 = f[, 2] + rnorm(n, 0, .5), y2 = 0.9 * f[, 2] + rnorm(n, 0, .5),
    y3 = 1.1 * f[, 2] + rnorm(n, 0, .5))
}

test_that("correlated-factor CFA matches lavaan exactly", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # cross-package agreement is CI-only
  d <- sim_two_factor()
  fit <- fit_sem(measurement = list(f1 = ~ x1 + x2 + x3,
                                    f2 = ~ y1 + y2 + y3), data = d)
  lav <- lavaan::cfa("f1 =~ x1+x2+x3\nf2 =~ y1+y2+y3", data = d)
  pe <- lavaan::parameterEstimates(lav)

  expect_equal(fit$logLik, as.numeric(lavaan::fitMeasures(lav, "logl")), tolerance = 1e-5)
  i <- fit$param_table$label == "f1~~f2"
  j <- pe$lhs == "f1" & pe$op == "~~" & pe$rhs == "f2"
  expect_equal(fit$param_table$est[i], pe$est[j], tolerance = 1e-4)
  expect_equal(fit$param_table$se[i], pe$se[j], tolerance = 5e-3)
  # Off-diagonal latent covariance present in the output matrix
  expect_gt(fit$latent_covariance["f1", "f2"], 0.3)
})

test_that("fit measures match lavaan", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # cross-package agreement is CI-only
  d <- sim_two_factor(seed = 7)
  fit <- fit_sem(measurement = list(f1 = ~ x1 + x2 + x3,
                                    f2 = ~ y1 + y2 + y3), data = d)
  lav <- lavaan::cfa("f1 =~ x1+x2+x3\nf2 =~ y1+y2+y3", data = d)
  fm_l <- lavaan::fitMeasures(lav, c("chisq", "df", "cfi", "tli", "rmsea",
                                     "rmsea.ci.lower", "rmsea.ci.upper",
                                     "srmr"))
  fm <- fit$fit_measures
  expect_equal(unname(fm["chisq"]), unname(fm_l["chisq"]), tolerance = 1e-3)
  expect_equal(unname(fm["df"]), unname(fm_l["df"]))
  expect_equal(unname(fm["cfi"]), unname(fm_l["cfi"]), tolerance = 1e-4)
  expect_equal(unname(fm["tli"]), unname(fm_l["tli"]), tolerance = 1e-3)
  expect_equal(unname(fm["rmsea"]), unname(fm_l["rmsea"]), tolerance = 1e-4)
  expect_equal(unname(fm["rmsea_ci_upper"]),
               unname(fm_l["rmsea.ci.upper"]), tolerance = 1e-3)
  expect_equal(unname(fm["srmr"]), unname(fm_l["srmr"]), tolerance = 1e-4)
})

test_that("MIMIC models match lavaan with fixed.x = FALSE", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # cross-package agreement is CI-only
  set.seed(42)
  n <- 500
  w <- rnorm(n)
  f1 <- 0.5 * w + rnorm(n)
  f2 <- 0.6 * f1 + 0.3 * w + rnorm(n, 0, .8)
  d <- data.frame(
    x1 = f1 + rnorm(n, 0, .6), x2 = 0.8 * f1 + rnorm(n, 0, .6),
    x3 = 1.2 * f1 + rnorm(n, 0, .6),
    y1 = f2 + rnorm(n, 0, .5), y2 = 0.9 * f2 + rnorm(n, 0, .5),
    y3 = 1.1 * f2 + rnorm(n, 0, .5), w = w)

  fit <- fit_sem(measurement = list(f1 = ~ x1 + x2 + x3,
                                    f2 = ~ y1 + y2 + y3),
                 structural = list(f1 ~ w, f2 ~ f1 + w), data = d)
  lav <- lavaan::sem("f1 =~ x1+x2+x3\nf2 =~ y1+y2+y3\nf1 ~ w\nf2 ~ f1 + w",
                     data = d, fixed.x = FALSE)
  pe <- lavaan::parameterEstimates(lav)

  expect_equal(fit$logLik, as.numeric(lavaan::fitMeasures(lav, "logl")), tolerance = 1e-4)
  for (lab in c("f1~w", "f2~f1", "f2~w")) {
    parts <- strsplit(lab, "~")[[1]]
    i <- fit$param_table$label == lab
    j <- pe$lhs == parts[1] & pe$op == "~" & pe$rhs == parts[2]
    expect_equal(fit$param_table$est[i], pe$est[j], tolerance = 1e-3)
    expect_equal(fit$param_table$se[i], pe$se[j], tolerance = 5e-3)
  }
  expect_identical(fit$covariates, "w")
})

test_that("FIML matches lavaan missing = 'fiml'", {
  skip_if_not_installed("lavaan")
  skip_on_cran()  # cross-package agreement is CI-only
  d <- sim_two_factor(n = 500, seed = 9)
  set.seed(10)
  for (v in c("x1", "x2", "y1", "y3")) d[[v]][sample(500, 60)] <- NA

  fit <- fit_sem(measurement = list(f1 = ~ x1 + x2 + x3,
                                    f2 = ~ y1 + y2 + y3),
                 structural = list(f2 ~ f1), data = d, missing = "fiml")
  lav <- lavaan::sem("f1 =~ x1+x2+x3\nf2 =~ y1+y2+y3\nf2 ~ f1", data = d,
                     missing = "fiml", fixed.x = FALSE)
  pe <- lavaan::parameterEstimates(lav)

  expect_equal(fit$logLik, as.numeric(lavaan::fitMeasures(lav, "logl")), tolerance = 1e-3)
  i <- fit$param_table$label == "f2~f1"
  j <- pe$lhs == "f2" & pe$op == "~" & pe$rhs == "f1"
  expect_equal(fit$param_table$est[i], pe$est[j], tolerance = 1e-3)
  expect_equal(fit$param_table$se[i], pe$se[j], tolerance = 5e-3)
  expect_equal(unname(fit$fit_measures["chisq"]),
               as.numeric(lavaan::fitMeasures(lav, "chisq")), tolerance = 1e-2)
  # All rows retained (no listwise deletion)
  expect_equal(fit$n_obs, 500)
  expect_identical(fit$missing, "fiml")
  # Factor scores computed for every case despite missingness
  expect_false(anyNA(fit$factor_scores))
})

test_that("listwise default warns and drops incomplete rows", {
  d <- sim_two_factor(n = 300, seed = 11)
  d$x1[1:30] <- NA
  expect_warning(
    fit <- fit_sem(measurement = list(f1 = ~ x1 + x2 + x3,
                                      f2 = ~ y1 + y2 + y3), data = d),
    "listwise")
  expect_equal(fit$n_obs, 270)
})

test_that("standardized solution bounds and summary run", {
  d <- sim_two_factor(seed = 13)
  fit <- fit_sem(measurement = list(f1 = ~ x1 + x2 + x3,
                                    f2 = ~ y1 + y2 + y3), data = d)
  std <- fit$standardized
  # Standardized loadings and the factor correlation are within [-1, 1]
  expect_true(all(abs(std$est_std[grepl("=~|~~", std$label) &
                                    !grepl("(x|y)\\d~~", std$label)]) <= 1.001))
  out <- capture.output(summary(fit))
  expect_true(any(grepl("RMSEA", out)))
  expect_true(any(grepl("Parameter estimates", out)))
})

test_that("laplace path warns about orthogonal exogenous factors", {
  d <- sim_two_factor(seed = 17, n = 200)
  expect_warning(
    fit_sem(measurement = list(f1 = ~ x1 + x2 + x3, f2 = ~ y1 + y2 + y3),
            data = d, method = "laplace"),
    "uncorrelated")
  expect_error(
    fit_sem(measurement = list(f1 = ~ x1 + x2 + x3),
            structural = list(f1 ~ x1), data = d, method = "laplace"),
    "latent")
})
