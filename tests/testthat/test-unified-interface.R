test_that("gllamm() provides unified interface for all models", {
  set.seed(12345)
  n <- 60

  # Test 1: Gaussian GLMM (original functionality)
  data_gauss <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    group = rep(1:6, each = 10)
  )

  fit_gauss <- gllamm(y ~ x + (1 | group), data = data_gauss)
  expect_s3_class(fit_gauss, "gllamm")
  expect_equal(fit_gauss$family$family, "gaussian")


  # Test 2: Ordinal model via gllamm()
  data_ord <- data.frame(
    y = ordered(sample(1:4, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:6, each = 10)
  )

  fit_ord <- gllamm(y ~ x + (1 | group),
                   data = data_ord,
                   family = ordinal(link = "logit"))

  expect_s3_class(fit_ord, "gllamm_ordinal")
  expect_s3_class(fit_ord, "gllamm")  # Should inherit from gllamm
  expect_equal(fit_ord$link, "logit")
  expect_equal(fit_ord$n_categories, 4)


  # Test 3: Ordinal with different link
  fit_ord_acl <- gllamm(y ~ x + (1 | group),
                       data = data_ord,
                       family = ordinal(link = "acl"))

  expect_s3_class(fit_ord_acl, "gllamm_ordinal")
  expect_equal(fit_ord_acl$link, "acl")


  # Test 4: fit() works with unified interface
  fit_stats_gauss <- fit(fit_gauss, quiet = TRUE)
  expect_s3_class(fit_stats_gauss, "fit_gllamm")

  fit_stats_ord <- fit(fit_ord, test_po = FALSE)
  expect_s3_class(fit_stats_ord, "fit_ordinal")


  # Test 5: plot() dispatches correctly
  expect_silent({
    pdf(tempfile())
    plot(fit_gauss, which = 1)
    dev.off()
  })

  expect_silent({
    pdf(tempfile())
    plot(fit_ord, which = 1, covariate = "x")
    dev.off()
  })
})


test_that("gllamm() and fit_ordinal() produce equivalent results", {
  set.seed(999)
  n <- 50

  data <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:5, each = 10)
  )

  # Fit via gllamm()
  fit1 <- gllamm(y ~ x + (1 | group),
                data = data,
                family = ordinal(link = "logit"))

  # Fit via fit_ordinal()
  fit2 <- fit_ordinal(y ~ x + (1 | group),
                     data = data,
                     link = "logit")

  # Should produce same results
  expect_equal(fit1$logLik, fit2$logLik, tolerance = 1e-6)
  expect_equal(fit1$coefficients$fixed, fit2$coefficients$fixed, tolerance = 1e-6)
  expect_equal(fit1$coefficients$thresholds, fit2$coefficients$thresholds, tolerance = 1e-6)
  expect_equal(fit1$n_categories, fit2$n_categories)
})


test_that("all ordinal links work via gllamm() interface", {
  set.seed(777)
  n <- 40

  data <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:4, each = 10)
  )

  # Test each link function
  links <- c("logit", "probit", "acl")

  for (link in links) {
    fit <- gllamm(y ~ x + (1 | group),
                 data = data,
                 family = ordinal(link = link))

    expect_s3_class(fit, "gllamm_ordinal")
    expect_equal(fit$link, link)
    expect_true(fit$convergence$converged)
  }
})


test_that("unified interface works with all methods", {
  set.seed(333)
  n <- 40

  data <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:4, each = 10)
  )

  fit <- gllamm(y ~ x + (1 | group),
               data = data,
               family = ordinal(link = "logit"))

  # Test S3 methods work
  expect_output(print(fit), "Ordinal Regression Model")
  expect_output(summary(fit), "Ordinal Regression Model")

  # Test fit() method
  fit_stats <- fit(fit, test_po = FALSE)
  expect_s3_class(fit_stats, "fit_ordinal")
  expect_output(print(fit_stats), "Ordinal Model Fit")

  # Test plot() method
  expect_silent({
    pdf(tempfile())
    plot(fit, which = 1, covariate = "x")
    dev.off()
  })

  # Test coef() if available
  coefs <- coef(fit)
  expect_true(!is.null(coefs))
})
