test_that("LCA data simulation works", {
  set.seed(123)
  n <- 200
  n_items <- 5

  # Two classes with different response patterns
  class1_probs <- c(0.8, 0.7, 0.9, 0.75, 0.85)
  class2_probs <- c(0.2, 0.3, 0.1, 0.25, 0.15)

  true_class <- sample(1:2, n, replace = TRUE, prob = c(0.6, 0.4))
  data <- matrix(NA, n, n_items)

  for (i in 1:n) {
    probs <- if (true_class[i] == 1) class1_probs else class2_probs
    data[i, ] <- rbinom(n_items, 1, probs)
  }

  expect_equal(dim(data), c(200, 5))
  expect_true(all(data %in% c(0, 1)))
})


test_that("fit_lca accepts valid input", {
  skip_if_not_installed("TMB")

  set.seed(456)
  n <- 100
  data <- matrix(rbinom(100 * 4, 1, 0.5), 100, 4)

  fit <- fit_lca(data, nclass = 2)

  expect_s3_class(fit, "gllamm_lca")
  expect_equal(fit$nclass, 2)
  expect_equal(length(fit$class_probs), 2)
  expect_equal(dim(fit$item_probs), c(4, 2))
})


test_that("LCA class probabilities sum to 1", {
  skip_if_not_installed("TMB")

  set.seed(789)
  data <- matrix(rbinom(150 * 6, 1, 0.5), 150, 6)

  fit <- fit_lca(data, nclass = 3)

  expect_equal(sum(fit$class_probs), 1, tolerance = 1e-6)
})


test_that("LCA recovers known classes", {
  skip_if_not_installed("TMB")

  set.seed(111)
  n <- 300
  n_items <- 6

  # Clear class separation
  class1_probs <- rep(0.9, n_items)
  class2_probs <- rep(0.1, n_items)

  true_class <- sample(1:2, n, replace = TRUE, prob = c(0.5, 0.5))
  data <- matrix(NA, n, n_items)

  for (i in 1:n) {
    probs <- if (true_class[i] == 1) class1_probs else class2_probs
    data[i, ] <- rbinom(n_items, 1, probs)
  }

  fit <- fit_lca(data, nclass = 2)

  # Check that modal class matches true class reasonably well
  # (may need to relabel classes)
  accuracy1 <- mean(fit$modal_class == true_class)
  accuracy2 <- mean(fit$modal_class == (3 - true_class))  # Flipped labels

  expect_true(max(accuracy1, accuracy2) > 0.7)
})


test_that("LCA posterior probabilities are valid", {
  skip_if_not_installed("TMB")

  set.seed(222)
  data <- matrix(rbinom(100 * 4, 1, 0.5), 100, 4)

  fit <- fit_lca(data, nclass = 2)

  # Posterior probabilities should sum to 1 for each person
  row_sums <- rowSums(fit$posterior)
  expect_true(all(abs(row_sums - 1) < 1e-6))

  # All probabilities should be in [0, 1]
  expect_true(all(fit$posterior >= 0))
  expect_true(all(fit$posterior <= 1))
})


test_that("LCA print and summary methods work", {
  skip_if_not_installed("TMB")

  set.seed(333)
  data <- matrix(rbinom(80 * 5, 1, 0.5), 80, 5)

  fit <- fit_lca(data, nclass = 2)

  expect_output(print(fit), "Latent Class Analysis")
  expect_output(print(fit), "Number of classes: 2")
  expect_output(summary(fit), "Class sizes")
})


test_that("LCA handles different numbers of classes", {
  skip_if_not_installed("TMB")

  set.seed(444)
  data <- matrix(rbinom(100 * 4, 1, 0.5), 100, 4)

  fit2 <- fit_lca(data, nclass = 2)
  fit3 <- fit_lca(data, nclass = 3)
  fit4 <- fit_lca(data, nclass = 4)

  expect_equal(fit2$nclass, 2)
  expect_equal(fit3$nclass, 3)
  expect_equal(fit4$nclass, 4)

  # BIC should penalize more complex models
  # (though with random data, not guaranteed)
  expect_true(fit4$BIC > fit2$BIC)
})


test_that("LCA item probabilities are in valid range", {
  skip_if_not_installed("TMB")

  set.seed(555)
  data <- matrix(rbinom(120 * 5, 1, 0.5), 120, 5)

  fit <- fit_lca(data, nclass = 2)

  expect_true(all(fit$item_probs >= 0))
  expect_true(all(fit$item_probs <= 1))
})
