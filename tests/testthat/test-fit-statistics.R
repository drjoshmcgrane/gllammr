test_that("fit() generic dispatches correctly", {
  # Create simple models of different types
  set.seed(123)
  n <- 50

  # GLMM
  data_glmm <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    group = rep(1:5, each = 10)
  )
  fit_glmm <- gllamm(y ~ x + (1 | group), data = data_glmm)

  # Ordinal (using gllamm() interface)
  data_ord <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:5, each = 10)
  )
  fit_ord <- gllamm(y ~ x + (1 | group), data = data_ord,
                   family = ordinal(link = "logit"))

  # Test that fit() works
  expect_silent(fit_stats_glmm <- fit(fit_glmm, quiet = TRUE))
  expect_silent(fit_stats_ord <- fit(fit_ord, test_po = FALSE))

  # Check classes
  expect_s3_class(fit_stats_glmm, "fit_statistics")
  expect_s3_class(fit_stats_glmm, "fit_gllamm")
  expect_s3_class(fit_stats_ord, "fit_statistics")
  expect_s3_class(fit_stats_ord, "fit_ordinal")
})


test_that("fit() for GLLAMM returns correct structure", {
  set.seed(456)
  n <- 60
  data <- data.frame(
    y = rnorm(n, mean = 10, sd = 2),
    x = rnorm(n),
    group = rep(1:6, each = 10)
  )

  fit_model <- gllamm(y ~ x + (1 | group), data = data)
  fit_stats <- fit(fit_model, quiet = TRUE)

  # Check components
  expect_true("model_type" %in% names(fit_stats))
  expect_equal(fit_stats$model_type, "GLMM")

  expect_true("logLik" %in% names(fit_stats))
  expect_true("AIC" %in% names(fit_stats))
  expect_true("BIC" %in% names(fit_stats))
  expect_true("n_obs" %in% names(fit_stats))
  expect_true("n_params" %in% names(fit_stats))

  expect_equal(fit_stats$n_obs, n)

  # For Gaussian models, should have R-squared
  expect_true("R2_marginal" %in% names(fit_stats))
  expect_true("R2_conditional" %in% names(fit_stats))
  expect_true(fit_stats$R2_marginal >= 0 && fit_stats$R2_marginal <= 1)
  expect_true(fit_stats$R2_conditional >= fit_stats$R2_marginal)

  # Should have ICC
  expect_true("ICC" %in% names(fit_stats))
})


test_that("fit() for ordinal models returns correct structure", {
  set.seed(789)
  n <- 50
  data <- data.frame(
    y = ordered(sample(1:4, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:5, each = 10)
  )

  fit_model <- fit_ordinal(y ~ x + (1 | group), data = data, link = "logit")
  fit_stats <- fit(fit_model, test_po = FALSE)  # Skip PO test for speed

  # Check components
  expect_equal(fit_stats$model_type, "Ordinal")
  expect_equal(fit_stats$link, "logit")

  expect_true("logLik" %in% names(fit_stats))
  expect_true("AIC" %in% names(fit_stats))
  expect_true("BIC" %in% names(fit_stats))
  expect_true("n_categories" %in% names(fit_stats))

  expect_equal(fit_stats$n_categories, 4)
  expect_equal(fit_stats$n_obs, n)

  # Should have pseudo-R²
  expect_true("pseudo_R2" %in% names(fit_stats))
  expect_true(fit_stats$pseudo_R2 >= 0 && fit_stats$pseudo_R2 <= 1)
})


test_that("fit() for ordinal models includes PO test when requested", {
  set.seed(101112)
  n <- 40
  data <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:4, each = 10)
  )

  fit_model <- fit_ordinal(y ~ x + (1 | group), data = data, link = "logit")

  # With PO test (default for logit/probit): refits the model as PPO and
  # reports a likelihood ratio test (prints a progress message)
  expect_output({
    fit_stats <- fit(fit_model, test_po = TRUE)
  }, "partial proportional odds")

  expect_true("proportional_odds_test" %in% names(fit_stats))
  expect_s3_class(fit_stats$proportional_odds_test, "po_test")
  expect_true(is.finite(fit_stats$proportional_odds_test$p_value))

  # Without PO test
  fit_stats_no_po <- fit(fit_model, test_po = FALSE)
  expect_false("proportional_odds_test" %in% names(fit_stats_no_po))
})


test_that("fit() pseudo-R² is computed correctly for ordinal models", {
  set.seed(131415)
  n <- 60

  # Model with strong effect (should have higher pseudo-R²)
  data_strong <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE, prob = c(0.2, 0.3, 0.5))),
    x = rnorm(n, mean = 0, sd = 2),  # Strong predictor
    group = rep(1:6, each = 10)
  )
  # Make y depend on x
  data_strong$y <- ordered(pmin(3, pmax(1, round(2 + 0.5 * data_strong$x))))

  fit_strong <- fit_ordinal(y ~ x + (1 | group), data = data_strong)
  fit_stats_strong <- fit(fit_strong, test_po = FALSE)

  # Pseudo-R\^2 should be positive for model with predictor
  expect_true(fit_stats_strong$pseudo_R2 >= 0)
  expect_true(fit_stats_strong$pseudo_R2 <= 1)
})


test_that("print.fit_statistics works for different model types", {
  set.seed(161718)
  n <- 40

  # GLMM
  data_glmm <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    group = rep(1:4, each = 10)
  )
  fit_glmm <- gllamm(y ~ x + (1 | group), data = data_glmm)
  fit_stats_glmm <- fit(fit_glmm, quiet = TRUE)

  expect_output(print(fit_stats_glmm), "Model Fit Statistics")
  expect_output(print(fit_stats_glmm), "Model type: GLMM")
  expect_output(print(fit_stats_glmm), "R\\^2")

  # Ordinal
  data_ord <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n),
    group = rep(1:4, each = 10)
  )
  fit_ord <- fit_ordinal(y ~ x + (1 | group), data = data_ord)
  fit_stats_ord <- fit(fit_ord, test_po = FALSE)

  expect_output(print(fit_stats_ord), "Model type: Ordinal")
  expect_output(print(fit_stats_ord), "Pseudo-R.2")
})


test_that("fit() errors for unsupported classes", {
  # Create a non-GLLAMM object
  lm_fit <- lm(mpg ~ wt, data = mtcars)

  expect_error(
    fit(lm_fit),
    "No fit method for class"
  )
})


test_that("fit() handles models without random effects", {
  # gllamm() has no fixed-effects-only path (by design: use glm()/lm()
  # for that); this documents the actual behavior rather than a
  # hypothetical no-RE fit() output
  set.seed(181920)
  n <- 50
  data <- data.frame(y = rnorm(n), x = rnorm(n))
  expect_error(suppressWarnings(gllamm(y ~ x, data = data)),
               "No random effects specified")
})


test_that("fit() computes ICC correctly", {
  set.seed(192021)
  n <- 60
  data <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    group = rep(1:6, each = 10)
  )

  fit_model <- gllamm(y ~ x + (1 | group), data = data)
  fit_stats <- fit(fit_model, quiet = TRUE)

  # ICC should be between 0 and 1
  expect_true(all(fit_stats$ICC >= 0))
  expect_true(all(fit_stats$ICC <= 1))
})


test_that("fit() for LCA models returns correct structure", {
  skip_if_not_installed("TMB")

  set.seed(3031)
  n <- 150; n_items <- 5
  class1_probs <- c(0.8, 0.7, 0.9, 0.75, 0.85)
  class2_probs <- c(0.2, 0.3, 0.1, 0.25, 0.15)
  true_class <- sample(1:2, n, replace = TRUE, prob = c(0.6, 0.4))
  data_lca <- matrix(NA, n, n_items)
  for (i in seq_len(n)) {
    probs <- if (true_class[i] == 1) class1_probs else class2_probs
    data_lca[i, ] <- rbinom(n_items, 1, probs)
  }
  fit_lca_model <- fit_lca(data_lca, nclass = 2)
  fit_stats <- fit(fit_lca_model)

  expect_true("entropy" %in% names(fit_stats))
  expect_true("class_proportions" %in% names(fit_stats))
  expect_true(fit_stats$entropy >= 0 && fit_stats$entropy <= 1)
  expect_equal(length(fit_stats$class_proportions), 2)
  expect_equal(sum(fit_stats$class_proportions), 1, tolerance = 1e-8)
})


test_that("fit() for IRT models returns correct structure", {
  skip_if_not_installed("TMB")

  set.seed(3233)
  n_persons <- 150; n_items <- 8
  theta <- rnorm(n_persons)
  b <- seq(-1.2, 1.2, length.out = n_items)
  responses <- sapply(seq_len(n_items), function(j)
    rbinom(n_persons, 1, plogis(theta - b[j])))
  fit_irt_model <- fit_irt(responses, model = "2PL")
  fit_stats <- fit(fit_irt_model)

  expect_true("item_fit" %in% names(fit_stats))
  expect_true("person_fit" %in% names(fit_stats))
  expect_true("reliability" %in% names(fit_stats))
  expect_equal(nrow(fit_stats$item_fit), n_items)
  expect_equal(nrow(fit_stats$person_fit), n_persons)
})


test_that("fit statistics are consistent with model object", {
  set.seed(222324)
  n <- 50
  data <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    group = rep(1:5, each = 10)
  )

  fit_model <- gllamm(y ~ x + (1 | group), data = data)
  fit_stats <- fit(fit_model, quiet = TRUE)

  # Check that fit statistics match model object
  expect_equal(fit_stats$logLik, fit_model$logLik)
  expect_equal(fit_stats$AIC, fit_model$AIC)
  expect_equal(fit_stats$BIC, fit_model$BIC)
  expect_equal(fit_stats$n_obs, fit_model$n_obs)
})
