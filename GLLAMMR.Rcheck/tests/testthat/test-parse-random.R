# Test random effects formula parsing

test_that("parse_random_formula handles simple grouping", {
  data <- data.frame(
    student_id = 1:100,
    class_id = rep(1:10, each = 10)
  )

  formula <- ~ (1 | class_id)
  result <- parse_random_formula(formula, data)

  expect_equal(length(result), 1)
  expect_equal(result[[1]]$group_var, "class_id")
  expect_equal(result[[1]]$type, "intercept")
  expect_false(result[[1]]$is_nested)
})

test_that("parse_random_formula handles nested notation", {
  data <- data.frame(
    student_id = 1:100,
    school_id = rep(1:5, each = 20),
    class_id = rep(1:20, each = 5)
  )

  formula <- ~ (1 | school_id/class_id)
  result <- parse_random_formula(formula, data)

  # Should expand to school + school:class
  expect_equal(length(result), 2)
  expect_equal(result[[1]]$group_var, "school_id")
  expect_equal(result[[2]]$group_var, "school_id:class_id")
  expect_true(result[[2]]$is_interaction %||% FALSE)
})

test_that("parse_random_formula handles multiple terms", {
  data <- data.frame(
    student_id = 1:100,
    class_id = rep(1:10, each = 10),
    school_id = rep(1:5, each = 20)
  )

  formula <- ~ (1 | school_id) + (1 | class_id)
  result <- parse_random_formula(formula, data)

  expect_equal(length(result), 2)
  expect_equal(result[[1]]$group_var, "school_id")
  expect_equal(result[[2]]$group_var, "class_id")
})

test_that("parse_random_formula validates grouping variables", {
  data <- data.frame(
    student_id = 1:100,
    class_id = rep(1:10, each = 10)
  )

  formula <- ~ (1 | nonexistent)
  expect_error(parse_random_formula(formula, data), "not found in person_data")
})

test_that("create_grouping_matrix works correctly", {
  data <- data.frame(
    student_id = 1:100,
    class_id = rep(1:10, each = 10)
  )

  formula <- ~ (1 | class_id)
  terms <- parse_random_formula(formula, data)
  result <- create_grouping_matrix(terms, data)

  expect_equal(result$n_re, 1)
  expect_equal(result$n_groups, 10)
  expect_equal(result$group_names, "class_id")
  expect_equal(dim(result$group_ids), c(100, 1))

  # Check 0-indexed
  expect_true(all(result$group_ids >= 0 & result$group_ids <= 9))

  # Check grouping structure
  expect_equal(result$group_ids[1:10, 1], rep(0, 10))
  expect_equal(result$group_ids[11:20, 1], rep(1, 10))
})

test_that("create_grouping_matrix handles NA values", {
  data <- data.frame(
    student_id = 1:100,
    class_id = c(rep(1:8, each = 10), rep(NA, 20))
  )

  formula <- ~ (1 | class_id)
  terms <- parse_random_formula(formula, data)
  result <- create_grouping_matrix(terms, data)

  # NA should be coded as -1
  expect_equal(result$group_ids[81:100, 1], rep(-1L, 20))
  expect_equal(result$n_groups, 8)
})

test_that("extract_random_terms handles parentheses", {
  # This was a bug that was fixed
  formula <- ~ (1 | class_id)
  result <- extract_random_terms(formula)

  expect_equal(length(result), 1)
  expect_equal(as.character(result[[1]][[3]]), "class_id")
})

test_that("expand_nested_terms creates interaction variables", {
  data <- data.frame(
    student_id = 1:40,
    school_id = rep(1:2, each = 20),
    class_id = rep(1:4, each = 10)
  )

  # Create nested term manually
  nested_term <- list(
    type = "intercept",
    group_var = "school_id",
    nested_in = "class_id",
    is_nested = TRUE
  )

  result <- expand_nested_terms(list(nested_term), data)

  expect_equal(length(result), 2)
  expect_equal(result[[1]]$group_var, "school_id")
  expect_equal(result[[2]]$group_var, "school_id:class_id")
  expect_true(result[[2]]$is_interaction %||% FALSE)
})
