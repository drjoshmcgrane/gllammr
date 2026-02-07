test_that("DIF test accepts valid input", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(123)
  n_persons <- 200
  n_items <- 10

  # Simulate data without DIF
  theta <- rnorm(n_persons, 0, 1)
  difficulty <- rnorm(n_items, 0, 1)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      p <- plogis(theta[i] - difficulty[j])
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  # Create group variable
  group <- rep(c("A", "B"), each = n_persons/2)

  # Test for DIF
  dif_result <- dif_test_with_data(responses, group = group, model = "Rasch")

  expect_s3_class(dif_result, "dif_analysis")
  expect_equal(nrow(dif_result$dif_results), n_items)
})


test_that("DIF detects uniform DIF in simulated data", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(456)
  n_persons <- 300
  n_items <- 10

  # Simulate data WITH uniform DIF on item 5
  theta <- rnorm(n_persons, 0, 1)
  difficulty <- rnorm(n_items, 0, 1)

  group <- rep(c(1, 2), each = n_persons/2)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      # Add DIF to item 5 for group 2
      diff_effective <- difficulty[j]
      if (j == 5 && group[i] == 2) {
        diff_effective <- diff_effective + 1.5  # Large uniform DIF
      }

      p <- plogis(theta[i] - diff_effective)
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  # Test for DIF
  dif_result <- dif_test_with_data(responses, group = group, model = "Rasch", alpha = 0.05)

  # Item 5 should be flagged
  expect_true(5 %in% dif_result$flagged_items)

  # Item 5 should have low p-value
  item5_pval <- dif_result$dif_results$p_value[dif_result$dif_results$item == 5]
  expect_lt(item5_pval, 0.05)
})


test_that("DIF test with 2PL model", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(789)
  n_persons <- 200
  n_items <- 8

  theta <- rnorm(n_persons, 0, 1)
  difficulty <- rnorm(n_items, 0, 1)
  discrimination <- runif(n_items, 0.5, 2)

  group <- rep(c("Male", "Female"), each = n_persons/2)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      p <- plogis(discrimination[j] * (theta[i] - difficulty[j]))
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  dif_result <- dif_test_with_data(responses, group = group, model = "2PL")

  expect_s3_class(dif_result, "dif_analysis")
  expect_equal(dif_result$model, "2PL")
  expect_equal(dif_result$group_labels, c("Male", "Female"))
})


test_that("DIF test validates group variable", {
  skip_if_not_installed("TMB")

  set.seed(111)
  n_persons <- 100
  n_items <- 5

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)

  # Invalid: 3 groups
  group_invalid <- rep(c("A", "B", "C"), length.out = n_persons)

  expect_error(dif_test_with_data(responses, group = group_invalid, model = "Rasch"),
               "must have exactly 2 unique values")

  # Invalid: wrong length
  group_wrong_length <- rep(c("A", "B"), each = 25)

  expect_error(dif_test_with_data(responses, group = group_wrong_length, model = "Rasch"),
               "must match number of persons")
})


test_that("DIF test with polytomous model (GRM)", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(222)
  n_persons <- 200
  n_items <- 8
  n_categories <- 4

  responses <- matrix(sample(1:n_categories, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  group <- rep(c(1, 2), each = n_persons/2)

  dif_result <- dif_test_with_data(responses, group = group, model = "GRM")

  expect_s3_class(dif_result, "dif_analysis")
  expect_equal(dif_result$model, "GRM")
})


test_that("DIF print method works", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(333)
  n_persons <- 150
  n_items <- 8

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)
  group <- rep(c("Control", "Treatment"), each = n_persons/2)

  dif_result <- dif_test_with_data(responses, group = group, model = "Rasch")

  # Print should not error
  expect_output(print(dif_result), "Differential Item Functioning")
  expect_output(print(dif_result), "Rasch")
  expect_output(print(dif_result), "Control vs Treatment")
})


test_that("DIF summary method works", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(444)
  n_persons <- 100
  n_items <- 6

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)
  group <- rep(c(0, 1), each = n_persons/2)

  dif_result <- dif_test_with_data(responses, group = group, model = "Rasch")

  expect_output(summary(dif_result), "All Items")
})


test_that("DIF plot for dichotomous items", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(555)
  n_persons <- 200
  n_items <- 8

  theta <- rnorm(n_persons, 0, 1)
  difficulty <- rnorm(n_items, 0, 1)
  discrimination <- rep(1.5, n_items)

  group <- rep(c(1, 2), each = n_persons/2)

  # Add DIF to item 3
  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      diff_eff <- difficulty[j]
      if (j == 3 && group[i] == 2) {
        diff_eff <- diff_eff + 1
      }
      p <- plogis(discrimination[j] * (theta[i] - diff_eff))
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  dif_result <- dif_test_with_data(responses, group = group, model = "2PL")

  # Plot should not error
  expect_silent(dif_plot(dif_result, item = 3))
})


test_that("DIF effect size computation", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(666)
  n_persons <- 200
  n_items <- 10

  theta <- rnorm(n_persons, 0, 1)
  difficulty <- rnorm(n_items, 0, 1)

  group <- rep(c(1, 2), each = n_persons/2)

  # Large DIF on item 7
  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      diff_eff <- difficulty[j]
      if (j == 7 && group[i] == 2) {
        diff_eff <- diff_eff + 2  # Very large DIF
      }
      p <- plogis(theta[i] - diff_eff)
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  dif_result <- dif_test_with_data(responses, group = group, model = "Rasch")

  # Item 7 should have large effect size
  item7_effect <- dif_result$dif_results$effect_size[dif_result$dif_results$item == 7]
  expect_gt(abs(item7_effect), 0.5)
})


test_that("DIF test with no flagged items", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(777)
  n_persons <- 150
  n_items <- 10

  # Simulate data WITHOUT DIF
  theta <- rnorm(n_persons, 0, 1)
  difficulty <- rnorm(n_items, 0, 1)

  group <- rep(c(1, 2), each = n_persons/2)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      p <- plogis(theta[i] - difficulty[j])
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  dif_result <- dif_test_with_data(responses, group = group, model = "Rasch", alpha = 0.01)

  # With strict alpha, should have few/no flagged items
  expect_lte(length(dif_result$flagged_items), 2)
})


test_that("DIF test with custom alpha level", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(888)
  n_persons <- 100
  n_items <- 8

  responses <- matrix(rbinom(n_persons * n_items, 1, 0.5), n_persons, n_items)
  group <- rep(c(1, 2), each = n_persons/2)

  # Test with different alpha levels
  dif_result_05 <- dif_test_with_data(responses, group = group, model = "Rasch", alpha = 0.05)
  dif_result_01 <- dif_test_with_data(responses, group = group, model = "Rasch", alpha = 0.01)

  # Stricter alpha should flag fewer items
  expect_lte(length(dif_result_01$flagged_items),
             length(dif_result_05$flagged_items))
})
