test_that("plot.gllamm dispatches to model-specific methods", {
  # This test checks that the dispatcher works correctly
  set.seed(123)
  n <- 50

  # Standard GLMM - should use default diagnostics
  data_glmm <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    group = rep(1:5, each = 10)
  )
  fit_glmm <- gllamm(y ~ x + (1 | group), data = data_glmm)

  # Should not error
  expect_silent({
    pdf(tempfile())
    plot(fit_glmm, which = 1:3)
    dev.off()
  })

  # Ordinal - should dispatch to plot.gllamm_ordinal
  data_ord <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:5, each = 10)
  )
  # Use gllamm() interface
  fit_ord <- gllamm(y ~ x + (1 | group), data = data_ord,
                   family = ordinal(link = "logit"))

  # Should not error
  expect_silent({
    pdf(tempfile())
    plot(fit_ord, which = 1:3, covariate = "x")
    dev.off()
  })
})


test_that("plot.gllamm_ordinal produces all plot types", {
  set.seed(456)
  n <- 60
  data <- data.frame(
    y = ordered(sample(1:4, n, replace = TRUE)),
    x1 = rnorm(n),
    x2 = rnorm(n),
    group = rep(1:6, each = 10)
  )

  fit <- fit_ordinal(y ~ x1 + x2 + (1 | group), data = data)

  # Test each plot type individually
  expect_silent({
    pdf(tempfile())
    plot(fit, which = 1, covariate = "x1")  # Cumulative probs
    dev.off()
  })

  expect_silent({
    pdf(tempfile())
    plot(fit, which = 2, covariate = "x1")  # Category probs
    dev.off()
  })

  expect_silent({
    pdf(tempfile())
    plot(fit, which = 3)  # Thresholds
    dev.off()
  })

  expect_silent({
    pdf(tempfile())
    plot(fit, which = 4, covariate = "x1")  # Covariate effects
    dev.off()
  })

  # Test all together
  expect_silent({
    pdf(tempfile())
    plot(fit, which = 1:4, covariate = "x2")
    dev.off()
  })
})


test_that("plot.gllamm_ordinal errors with invalid covariate", {
  set.seed(789)
  n <- 40
  data <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:4, each = 10)
  )

  fit <- fit_ordinal(y ~ x + (1 | group), data = data)

  expect_error({
    pdf(tempfile())
    plot(fit, which = 1, covariate = "nonexistent")
    dev.off()
  }, "Covariate.*not found")
})


test_that("plot_ordinal_effects works", {
  set.seed(101112)
  n <- 50
  data <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n),
    group = rep(1:5, each = 10)
  )

  fit <- fit_ordinal(y ~ x1 + x2 + x3 + (1 | group), data = data)

  # Should not error
  expect_silent({
    pdf(tempfile())
    plot_ordinal_effects(fit)
    dev.off()
  })

  # Test sorting
  expect_silent({
    pdf(tempfile())
    plot_ordinal_effects(fit, sort_by = "magnitude")
    dev.off()
  })

  expect_silent({
    pdf(tempfile())
    plot_ordinal_effects(fit, sort_by = "name")
    dev.off()
  })
})


test_that("IRT plotting functions work", {
  skip_if_not_installed("TMB")

  set.seed(2024)
  n_persons <- 100
  n_items <- 6
  theta <- rnorm(n_persons)
  b <- seq(-1, 1, length.out = n_items)
  responses <- sapply(seq_len(n_items), function(j)
    rbinom(n_persons, 1, plogis(theta - b[j])))

  fit_irt_2pl <- fit_irt(responses, model = "2PL")

  expect_silent({
    pdf(tempfile())
    plot(fit_irt_2pl, which = 1:4, items = 1:3)
    dev.off()
  })

  expect_silent({
    pdf(tempfile())
    plot(fit_irt_2pl, which = 1, items = 1:3)
    dev.off()
  })
})


test_that("LCA plotting functions work", {
  skip_if_not_installed("TMB")

  set.seed(2025)
  n <- 150
  n_items <- 5
  class1_probs <- c(0.8, 0.7, 0.9, 0.75, 0.85)
  class2_probs <- c(0.2, 0.3, 0.1, 0.25, 0.15)
  true_class <- sample(1:2, n, replace = TRUE, prob = c(0.6, 0.4))
  data_lca <- matrix(NA, n, n_items)
  for (i in seq_len(n)) {
    probs <- if (true_class[i] == 1) class1_probs else class2_probs
    data_lca[i, ] <- rbinom(n_items, 1, probs)
  }

  fit_lca_2 <- fit_lca(data_lca, nclass = 2)

  expect_silent({
    pdf(tempfile())
    plot(fit_lca_2, which = 1:3)
    dev.off()
  })

  expect_silent({
    pdf(tempfile())
    plot(fit_lca_2, which = 1)
    dev.off()
  })
})


test_that("plotting respects par settings", {
  set.seed(131415)
  n <- 40
  data <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:4, each = 10)
  )

  fit <- fit_ordinal(y ~ x + (1 | group), data = data)

  # Save par settings
  old_par <- par(no.readonly = TRUE)

  # Create plot
  pdf(tempfile())
  plot(fit, which = 1:3, covariate = "x")
  new_par <- par(no.readonly = TRUE)
  dev.off()

  # Par should be restored (mfrow in particular)
  expect_equal(old_par$mfrow, new_par$mfrow)
})


test_that("plot handles single plot request", {
  set.seed(161718)
  n <- 40
  data <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:4, each = 10)
  )

  fit <- fit_ordinal(y ~ x + (1 | group), data = data)

  # Single plot should not change par settings
  expect_silent({
    pdf(tempfile())
    plot(fit, which = 1, covariate = "x")
    dev.off()
  })
})


test_that("ordinal plots handle different numbers of categories", {
  set.seed(192021)
  n <- 50

  # 3 categories
  data3 <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:5, each = 10)
  )
  fit3 <- fit_ordinal(y ~ x + (1 | group), data = data3)

  expect_silent({
    pdf(tempfile())
    plot(fit3, which = 1:3, covariate = "x")
    dev.off()
  })

  # 5 categories
  data5 <- data.frame(
    y = ordered(sample(1:5, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:5, each = 10)
  )
  fit5 <- fit_ordinal(y ~ x + (1 | group), data = data5)

  expect_silent({
    pdf(tempfile())
    plot(fit5, which = 1:3, covariate = "x")
    dev.off()
  })
})


test_that("ordinal plots work with custom covariate values", {
  set.seed(222324)
  n <- 40
  data <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:4, each = 10)
  )

  fit <- fit_ordinal(y ~ x + (1 | group), data = data)

  # Custom covariate range
  expect_silent({
    pdf(tempfile())
    plot(fit, which = 1, covariate = "x", covariate_values = seq(-3, 3, length.out = 50))
    dev.off()
  })
})


test_that("default GLMM diagnostics still work", {
  set.seed(252627)
  n <- 50
  data <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    group = rep(1:5, each = 10)
  )

  fit <- gllamm(y ~ x + (1 | group), data = data)

  # Should produce standard diagnostic plots
  expect_silent({
    pdf(tempfile())
    plot(fit, which = c(1, 2, 3, 5))
    dev.off()
  })
})


test_that("plot automatically selects covariate when not specified", {
  set.seed(282930)
  n <- 40
  data <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x1 = rnorm(n),
    x2 = rnorm(n),
    group = rep(1:4, each = 10)
  )

  fit <- fit_ordinal(y ~ x1 + x2 + (1 | group), data = data)

  # Should automatically select first non-intercept covariate (x1)
  expect_silent({
    pdf(tempfile())
    plot(fit, which = 1)  # Don't specify covariate
    dev.off()
  })
})


test_that("plot errors for ordinal model with only intercept", {
  set.seed(313233)
  n <- 40
  data <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    group = rep(1:4, each = 10)
  )

  fit <- fit_ordinal(y ~ 1 + (1 | group), data = data)

  # Should error when trying to plot without covariate
  expect_error({
    pdf(tempfile())
    plot(fit, which = 1)
    dev.off()
  }, "No covariates found|Specify covariate")
})
