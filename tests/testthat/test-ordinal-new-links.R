test_that("ordinal() family constructor works", {
  # Test default (logit)
  fam1 <- ordinal()
  expect_s3_class(fam1, "ordinal_family")
  expect_equal(fam1$family, "ordinal")
  expect_equal(fam1$link, "logit")
  expect_equal(fam1$link_code, 1L)

  # Test probit
  fam2 <- ordinal("probit")
  expect_equal(fam2$link, "probit")
  expect_equal(fam2$link_code, 2L)

  # Test ACL
  fam3 <- ordinal("acl")
  expect_equal(fam3$link, "acl")
  expect_equal(fam3$link_code, 3L)

  # Test CRL forward
  fam4 <- ordinal("crl_forward")
  expect_equal(fam4$link, "crl_forward")
  expect_equal(fam4$link_code, 4L)

  # Test CRL backward
  fam5 <- ordinal("crl_backward")
  expect_equal(fam5$link, "crl_backward")
  expect_equal(fam5$link_code, 5L)

  # Test PPO
  fam6 <- ordinal("ppo")
  expect_equal(fam6$link, "ppo")
  expect_equal(fam6$link_code, 6L)

  # Test invalid link
  expect_error(ordinal("invalid"), "should be one of")
})


test_that("ordinal() family prints correctly", {
  fam <- ordinal("acl")
  output <- capture.output(print(fam))
  expect_true(any(grepl("Family: ordinal", output)))
  expect_true(any(grepl("Link function: acl", output)))
})


test_that("gllamm() works with ordinal() family", {
  skip_if_not_installed("MASS")

  # Simulate simple ordinal data
  set.seed(123)
  n <- 100
  n_groups <- 10

  data <- data.frame(
    y = sample(1:4, n, replace = TRUE),
    x = rnorm(n),
    group = rep(1:n_groups, each = n / n_groups)
  )
  data$y <- ordered(data$y)

  # Test gllamm() with ordinal(link = "logit")
  expect_silent({
    fit_via_gllamm <- gllamm(y ~ x + (1 | group),
                             data = data,
                             family = ordinal(link = "logit"))
  })
  expect_s3_class(fit_via_gllamm, "gllamm_ordinal")
  expect_equal(fit_via_gllamm$link, "logit")
  expect_equal(fit_via_gllamm$n_categories, 4)
  expect_true(length(fit_via_gllamm$coefficients$thresholds) == 3)

  # Test gllamm() with ordinal(link = "probit")
  expect_silent({
    fit_probit <- gllamm(y ~ x + (1 | group),
                        data = data,
                        family = ordinal(link = "probit"))
  })
  expect_s3_class(fit_probit, "gllamm_ordinal")
  expect_equal(fit_probit$link, "probit")

  # Test gllamm() with ordinal(link = "acl")
  expect_silent({
    fit_acl <- gllamm(y ~ x + (1 | group),
                     data = data,
                     family = ordinal(link = "acl"))
  })
  expect_s3_class(fit_acl, "gllamm_ordinal")
  expect_equal(fit_acl$link, "acl")
})


test_that("fit_ordinal works with standard links", {
  skip_if_not_installed("MASS")

  # Simulate simple ordinal data
  set.seed(123)
  n <- 100
  n_groups <- 10

  data <- data.frame(
    y = sample(1:4, n, replace = TRUE),
    x = rnorm(n),
    group = rep(1:n_groups, each = n / n_groups)
  )
  data$y <- ordered(data$y)

  # Test logit link via fit_ordinal() directly
  expect_silent({
    fit_logit <- fit_ordinal(y ~ x + (1 | group), data = data, link = "logit")
  })
  expect_s3_class(fit_logit, "gllamm_ordinal")
  expect_equal(fit_logit$link, "logit")
  expect_equal(fit_logit$n_categories, 4)
  expect_true(length(fit_logit$coefficients$thresholds) == 3)

  # Test probit link via fit_ordinal() directly
  expect_silent({
    fit_probit <- fit_ordinal(y ~ x + (1 | group), data = data, link = "probit")
  })
  expect_s3_class(fit_probit, "gllamm_ordinal")
  expect_equal(fit_probit$link, "probit")
})


test_that("test_proportional_odds returns correct structure", {
  # Create a simple ordinal model
  set.seed(456)
  n <- 50
  data <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:5, each = 10)
  )

  fit <- fit_ordinal(y ~ x + (1 | group), data = data, link = "logit")

  # PO test refits the model as PPO and runs a likelihood ratio test
  # (requires the data and prints a progress message)
  expect_output({
    po_test <- test_proportional_odds(fit, data = data)
  }, "partial proportional odds")

  expect_s3_class(po_test, "po_test")
  expect_true("statistic" %in% names(po_test))
  expect_true("df" %in% names(po_test))
  expect_true("p_value" %in% names(po_test))
  expect_true("conclusion" %in% names(po_test))

  # Test print method
  expect_output(print(po_test), "Proportional Odds Assumption Test")
})


test_that("test_proportional_odds errors for non-ordinal models", {
  # Create a regular GLMM
  data <- data.frame(
    y = rnorm(50),
    x = rnorm(50),
    group = rep(1:5, each = 10)
  )

  fit <- gllamm(y ~ x + (1 | group), data = data)

  expect_error(
    test_proportional_odds(fit),
    "must be of class 'gllamm_ordinal'"
  )
})


test_that("test_proportional_odds errors for ACL/CRL models", {
  skip_if_not_installed("MASS")

  set.seed(171819)
  n <- 200
  n_groups <- 5
  data <- data.frame(
    y = ordered(sample(1:4, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:n_groups, each = n / n_groups)
  )

  fit_acl <- fit_ordinal(y ~ x + (1 | group), data = data, link = "acl")
  expect_error(test_proportional_odds(fit_acl),
               "only applies to logit or probit")

  fit_crl <- fit_ordinal(y ~ x + (1 | group), data = data,
                         link = "crl_forward")
  expect_error(test_proportional_odds(fit_crl),
               "only applies to logit or probit")
})


test_that("ordinal model coefficient structure is correct", {
  set.seed(789)
  n <- 60
  data <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x1 = rnorm(n),
    x2 = rnorm(n),
    group = rep(1:6, each = 10)
  )

  fit <- fit_ordinal(y ~ x1 + x2 + (1 | group), data = data, link = "logit")

  # Check coefficient structure
  expect_true("fixed" %in% names(fit$coefficients))
  expect_true("thresholds" %in% names(fit$coefficients))
  expect_true("random_var" %in% names(fit$coefficients))

  # Fixed effects: the two covariates (intercept is absorbed into thresholds)
  expect_equal(length(fit$coefficients$fixed), 2)

  # Thresholds should be K-1 (3 categories = 2 thresholds)
  expect_equal(length(fit$coefficients$thresholds), 2)

  # Thresholds should be ordered
  expect_true(fit$coefficients$thresholds[1] < fit$coefficients$thresholds[2])
})


test_that("ordinal model convergence is tracked", {
  set.seed(101112)
  n <- 40
  data <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:4, each = 10)
  )

  fit <- fit_ordinal(y ~ x + (1 | group), data = data, link = "logit")

  expect_true("convergence" %in% names(fit))
  expect_true("converged" %in% names(fit$convergence))
  expect_true("message" %in% names(fit$convergence))
})


test_that("ordinal model handles different numbers of categories", {
  set.seed(131415)
  n <- 50

  # Test with 3 categories
  data3 <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:5, each = 10)
  )
  fit3 <- fit_ordinal(y ~ x + (1 | group), data = data3)
  expect_equal(fit3$n_categories, 3)
  expect_equal(length(fit3$coefficients$thresholds), 2)

  # Test with 5 categories
  data5 <- data.frame(
    y = ordered(sample(1:5, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:5, each = 10)
  )
  fit5 <- fit_ordinal(y ~ x + (1 | group), data = data5)
  expect_equal(fit5$n_categories, 5)
  expect_equal(length(fit5$coefficients$thresholds), 4)
})


test_that("ordinal model handles factor vs numeric response", {
  set.seed(161718)
  n <- 40

  # Test with ordered factor
  data_factor <- data.frame(
    y = ordered(sample(c("Low", "Med", "High"), n, replace = TRUE),
                levels = c("Low", "Med", "High")),
    x = rnorm(n),
    group = rep(1:4, each = 10)
  )
  fit_factor <- fit_ordinal(y ~ x + (1 | group), data = data_factor)
  expect_equal(fit_factor$n_categories, 3)
  expect_equal(fit_factor$category_labels, c("Low", "Med", "High"))

  # Test with numeric
  data_numeric <- data.frame(
    y = sample(1:3, n, replace = TRUE),
    x = rnorm(n),
    group = rep(1:4, each = 10)
  )
  fit_numeric <- fit_ordinal(y ~ x + (1 | group), data = data_numeric)
  expect_equal(fit_numeric$n_categories, 3)
})
