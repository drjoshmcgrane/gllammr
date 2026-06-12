# Test S3 methods for multi-level IRT

# Helper to create a simple multilevel fit
create_multilevel_fit <- function() {
  set.seed(123)
  n_persons <- 100
  n_items <- 10
  n_groups <- 10

  person_data <- data.frame(
    person_id = 1:n_persons,
    group_id = rep(1:n_groups, each = n_persons / n_groups)
  )

  theta_0 <- rnorm(n_persons, 0, 1)
  u_group <- rnorm(n_groups, 0, 0.5)
  theta_total <- theta_0 + u_group[person_data$group_id]
  difficulty <- rnorm(n_items, 0, 1)

  responses <- matrix(rbinom(n_persons * n_items, 1,
                             plogis(outer(theta_total, difficulty, '-'))),
                      n_persons, n_items)

  fit_irt(responses, model = "Rasch",
          person_data = person_data,
          random = ~ (1 | group_id))
}

test_that("VarCorr returns correct structure", {
  fit <- create_multilevel_fit()
  vc <- VarCorr(fit)

  expect_s3_class(vc, "VarCorr.gllamm")
  expect_s3_class(vc, "data.frame")
  expect_equal(ncol(vc), 3)  # Groups, Variance, Std.Dev
  expect_equal(nrow(vc), 3)  # group_id, Person, Residual
  expect_true(all(c("Groups", "Variance", "Std.Dev") %in% names(vc)))
})

test_that("VarCorr values are sensible", {
  fit <- create_multilevel_fit()
  vc <- VarCorr(fit)

  # All variances should be positive
  expect_true(all(vc$Variance > 0))

  # Std.Dev should be sqrt of Variance
  expect_equal(vc$Std.Dev, sqrt(vc$Variance))

  # Residual variance should be pi^2/3 for logistic
  expect_equal(vc$Variance[3], pi^2/3, tolerance = 1e-6)
})

test_that("print.VarCorr works", {
  fit <- create_multilevel_fit()
  vc <- VarCorr(fit)

  expect_output(print(vc), "Variance Components")
  expect_output(print(vc), "group_id")
  expect_output(print(vc), "Person")
  expect_output(print(vc), "Residual")
})

test_that("icc returns correct structure", {
  fit <- create_multilevel_fit()
  icc_all <- icc(fit)

  expect_type(icc_all, "double")
  expect_length(icc_all, 2)  # group_id + Person
  expect_equal(names(icc_all), c("group_id", "Person"))
})

test_that("icc values are in valid range", {
  fit <- create_multilevel_fit()
  icc_all <- icc(fit)

  # All ICCs should be between 0 and 1
  expect_true(all(icc_all >= 0))
  expect_true(all(icc_all <= 1))

  # Sum of ICCs should be less than 1 (residual variance exists)
  expect_true(sum(icc_all) < 1)
})

test_that("icc for specific level works", {
  fit <- create_multilevel_fit()

  icc_group <- icc(fit, level = "group_id")

  expect_length(icc_group, 1)
  expect_named(icc_group, "group_id")
  expect_true(icc_group >= 0 & icc_group <= 1)
})

test_that("icc errors for invalid level", {
  fit <- create_multilevel_fit()

  expect_error(icc(fit, level = "nonexistent"), "not found")
})

test_that("ranef returns correct structure", {
  fit <- create_multilevel_fit()
  re <- ranef(fit)

  expect_type(re, "list")
  expect_equal(names(re), "group_id")
  expect_length(re$group_id, 10)
  expect_type(re$group_id, "double")
})

test_that("ranef for specific level works", {
  fit <- create_multilevel_fit()

  re_group <- ranef(fit, level = "group_id")

  expect_length(re_group, 10)
  expect_type(re_group, "double")
  expect_true(all(!is.na(re_group)))
})

test_that("ranef errors for invalid level", {
  fit <- create_multilevel_fit()

  expect_error(ranef(fit, level = "nonexistent"), "not found")
})

test_that("abilities returns correct structure", {
  fit <- create_multilevel_fit()

  theta_0 <- abilities(fit, composite = FALSE)
  theta_comp <- abilities(fit, composite = TRUE)

  expect_length(theta_0, 100)
  expect_length(theta_comp, 100)
  expect_type(theta_0, "double")
  expect_type(theta_comp, "double")
})

test_that("composite abilities differ from person deviations", {
  fit <- create_multilevel_fit()

  theta_0 <- abilities(fit, composite = FALSE)
  theta_comp <- abilities(fit, composite = TRUE)

  # They should be different (group effects added to composite)
  expect_false(all(theta_0 == theta_comp))

  # But should be correlated (theoretical cor with sigma_theta = 1 and
  # sigma_u = 0.5 is 1/sqrt(1.25) ~ 0.89, lower still with estimation noise)
  expect_true(cor(theta_0, theta_comp) > 0.7)
})

test_that("coef.gllamm_irt_multilevel works for item type", {
  fit <- create_multilevel_fit()

  item_params <- coef(fit, type = "item")

  expect_s3_class(item_params, "data.frame")
  expect_equal(nrow(item_params), 10)
  expect_true("difficulty" %in% names(item_params))
})

test_that("coef.gllamm_irt_multilevel works for person type", {
  fit <- create_multilevel_fit()

  person_abilities <- coef(fit, type = "person")

  expect_length(person_abilities, 100)
  expect_type(person_abilities, "double")
})

test_that("coef.gllamm_irt_multilevel works for random type", {
  fit <- create_multilevel_fit()

  random_effects <- coef(fit, type = "random")

  expect_type(random_effects, "list")
  expect_equal(names(random_effects), "group_id")
  expect_length(random_effects$group_id, 10)
})

test_that("VarCorr fails for standard IRT model", {
  set.seed(123)
  responses <- matrix(rbinom(500, 1, 0.5), 50, 10)
  fit_std <- fit_irt(responses, model = "Rasch")

  expect_error(VarCorr(fit_std), "random effects")
})

test_that("icc fails for standard IRT model", {
  set.seed(123)
  responses <- matrix(rbinom(500, 1, 0.5), 50, 10)
  fit_std <- fit_irt(responses, model = "Rasch")

  expect_error(icc(fit_std), "random effects")
})

test_that("ranef fails for standard IRT model", {
  set.seed(123)
  responses <- matrix(rbinom(500, 1, 0.5), 50, 10)
  fit_std <- fit_irt(responses, model = "Rasch")

  expect_error(ranef(fit_std), "random effects")
})

test_that("abilities works for standard IRT model", {
  set.seed(123)
  responses <- matrix(rbinom(500, 1, 0.5), 50, 10)
  fit_std <- fit_irt(responses, model = "Rasch")

  theta <- abilities(fit_std)

  expect_length(theta, 50)
  expect_type(theta, "double")
})

test_that("Methods work with multiple random effects", {
  set.seed(123)
  n_persons <- 100
  n_items <- 10

  person_data <- data.frame(
    person_id = 1:n_persons,
    level1 = rep(1:5, each = 20),
    level2 = rep(1:20, each = 5)
  )

  theta <- rnorm(n_persons)
  difficulty <- rnorm(n_items)
  responses <- matrix(rbinom(n_persons * n_items, 1,
                             plogis(outer(theta, difficulty, '-'))),
                      n_persons, n_items)

  fit <- fit_irt(responses, model = "Rasch",
                 person_data = person_data,
                 random = ~ (1 | level1) + (1 | level2))

  # VarCorr should have 4 rows (2 REs + Person + Residual)
  vc <- VarCorr(fit)
  expect_equal(nrow(vc), 4)

  # icc should have 3 levels (2 REs + Person)
  icc_all <- icc(fit)
  expect_length(icc_all, 3)

  # ranef should return list with 2 elements
  re <- ranef(fit)
  expect_length(re, 2)
  expect_equal(names(re), c("level1", "level2"))

  # ranef for specific level
  re1 <- ranef(fit, level = "level1")
  expect_length(re1, 5)

  re2 <- ranef(fit, level = "level2")
  expect_length(re2, 20)
})
