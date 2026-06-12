test_that("Formula parsing works for simple random intercept", {
  data <- data.frame(
    y = rnorm(100),
    x = rnorm(100),
    group = rep(1:10, each = 10)
  )

  formula <- y ~ x + (1 | group)
  parsed <- parse_formula(formula, data)

  expect_equal(length(parsed$random_terms), 1)
  expect_equal(parsed$response_name, "y")
  expect_s3_class(parsed$fixed_formula, "formula")
})


test_that("Formula parsing handles multiple fixed effects", {
  data <- data.frame(
    y = rnorm(100),
    x1 = rnorm(100),
    x2 = rnorm(100),
    group = rep(1:10, each = 10)
  )

  formula <- y ~ x1 + x2 + (1 | group)
  parsed <- parse_formula(formula, data)

  expect_equal(length(all.vars(parsed$fixed_formula)), 3)  # y, x1, x2
})


test_that("Random term parsing extracts grouping variable", {
  data <- data.frame(y = 1:10, x = 1:10, g = rep(1:2, 5))
  rt <- parse_random_term("(1|g)", data)

  expect_equal(rt$grouping, "g")
  expect_s3_class(rt$formula, "formula")
})


test_that("Random term parsing detects nested structure", {
  data <- data.frame(y = 1:10, x = 1:10, school = rep(1:2, 5), class = 1:10)
  rt <- parse_random_term("(1|school/class)", data)

  expect_equal(rt$nested, TRUE)
  expect_equal(length(rt$grouping), 2)
})


test_that("Formula validation catches missing response", {
  data <- data.frame(x = 1:10, g = rep(1:2, 5))
  formula <- ~ x + (1 | g)

  expect_error(validate_formula(formula, data), "must have a response")
})


test_that("Formula validation catches missing variables", {
  data <- data.frame(y = 1:10, x = 1:10)
  formula <- y ~ x + (1 | group)

  expect_error(validate_formula(formula, data), "not found in data")
})


test_that("Model matrices are created correctly", {
  set.seed(123)
  data <- data.frame(
    y = rnorm(20),
    x = rnorm(20),
    group = rep(1:5, each = 4)
  )

  parsed <- parse_formula(y ~ x + (1 | group), data)
  mats <- make_model_matrices(parsed, data)

  expect_equal(mats$n_obs, 20)
  expect_equal(mats$n_fixed, 2)  # Intercept + x
  expect_equal(mats$n_groups[1], 5)
  expect_equal(ncol(mats$X), 2)
  expect_equal(nrow(mats$X), 20)
})


test_that("Grouping indices are 0-indexed", {
  data <- data.frame(
    y = 1:10,
    x = 1:10,
    group = rep(c("A", "B"), each = 5)
  )

  parsed <- parse_formula(y ~ x + (1 | group), data)
  mats <- make_model_matrices(parsed, data)

  # Should be 0 and 1, not 1 and 2
  expect_true(all(mats$groups[[1]] %in% c(0, 1)))
  expect_equal(min(mats$groups[[1]]), 0)
})
