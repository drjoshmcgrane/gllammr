# Test multi-level IRT models

# Helper function to simulate multi-level data
simulate_multilevel_data <- function(n_persons = 100, n_items = 10, n_groups = 10,
                                      sigma_person = 1.0, sigma_group = 0.5, seed = 123) {
  set.seed(seed)

  person_data <- data.frame(
    person_id = 1:n_persons,
    group_id = rep(1:n_groups, each = n_persons / n_groups)
  )

  theta_0 <- rnorm(n_persons, 0, sigma_person)
  u_group <- rnorm(n_groups, 0, sigma_group)
  theta_total <- theta_0 + u_group[person_data$group_id]
  difficulty <- rnorm(n_items, 0, 1)

  responses <- matrix(NA, n_persons, n_items)
  for (i in 1:n_persons) {
    for (j in 1:n_items) {
      p <- plogis(theta_total[i] - difficulty[j])
      responses[i, j] <- rbinom(1, 1, p)
    }
  }

  list(responses = responses, person_data = person_data,
       true_sigma_person = sigma_person, true_sigma_group = sigma_group)
}

test_that("fit_irt works without random effects (standard model)", {
  data <- simulate_multilevel_data(n_persons = 50, n_items = 8)

  fit <- fit_irt(data$responses, model = "Rasch", se = FALSE)

  expect_s3_class(fit, "gllamm_irt")
  expect_false(inherits(fit, "gllamm_irt_multilevel"))
  expect_null(fit$random_effects)
  expect_true(fit$convergence$converged)
})

test_that("fit_irt with random effects creates multilevel model", {
  data <- simulate_multilevel_data(n_persons = 100, n_items = 10)

  fit <- fit_irt(data$responses, model = "Rasch",
                 person_data = data$person_data,
                 random = ~ (1 | group_id), se = FALSE)

  expect_s3_class(fit, "gllamm_irt_multilevel")
  expect_s3_class(fit, "gllamm_irt")
  expect_false(is.null(fit$random_effects))
  expect_true(fit$convergence$converged)
})

test_that("multilevel model recovers variance components", {
  data <- simulate_multilevel_data(n_persons = 200, n_items = 15,
                                    sigma_person = 1.0, sigma_group = 0.5)

  fit <- fit_irt(data$responses, model = "Rasch",
                 person_data = data$person_data,
                 random = ~ (1 | group_id), se = FALSE)

  # Recovered SDs should be reasonably close to true values
  # (allowing for sampling variability)
  expect_true(abs(fit$ability_sd - data$true_sigma_person) < 0.3)
  expect_true(abs(fit$random_effects$sigma_random - data$true_sigma_group) < 0.3)
})

test_that("multilevel model improves fit over standard model", {
  data <- simulate_multilevel_data(n_persons = 100, n_items = 10,
                                    sigma_group = 0.6)  # Substantial group effect

  fit_std <- fit_irt(data$responses, model = "Rasch", se = FALSE)
  fit_ml <- fit_irt(data$responses, model = "Rasch",
                    person_data = data$person_data,
                    random = ~ (1 | group_id), se = FALSE)

  # Multilevel should have better (higher) log-likelihood
  expect_true(fit_ml$logLik >= fit_std$logLik)
})

test_that("multilevel model returns correct random effects structure", {
  data <- simulate_multilevel_data(n_persons = 100, n_items = 10, n_groups = 10)

  fit <- fit_irt(data$responses, model = "Rasch",
                 person_data = data$person_data,
                 random = ~ (1 | group_id), se = FALSE)

  expect_false(is.null(fit$random_effects))
  expect_equal(fit$random_effects$group_names, "group_id")
  expect_equal(fit$random_effects$n_groups, 10)
  expect_equal(dim(fit$random_effects$u_random)[2], 1)  # 1 RE level
  expect_length(fit$random_effects$sigma_random, 1)
  expect_length(fit$random_effects$icc, 2)  # group_id + Person
  expect_length(fit$random_effects$composite_theta, 100)
})

test_that("VarCorr extracts variance components correctly", {
  data <- simulate_multilevel_data(n_persons = 100, n_items = 10)

  fit <- fit_irt(data$responses, model = "Rasch",
                 person_data = data$person_data,
                 random = ~ (1 | group_id), se = FALSE)

  vc <- VarCorr(fit)

  expect_s3_class(vc, "VarCorr.gllamm")
  expect_equal(nrow(vc), 3)  # group_id, Person, Residual
  expect_equal(vc$Groups, c("group_id", "Person", "Residual"))
  expect_true(all(c("Variance", "Std.Dev") %in% names(vc)))
  expect_true(all(vc$Variance > 0))
})

test_that("icc computes intraclass correlations correctly", {
  data <- simulate_multilevel_data(n_persons = 100, n_items = 10)

  fit <- fit_irt(data$responses, model = "Rasch",
                 person_data = data$person_data,
                 random = ~ (1 | group_id), se = FALSE)

  # All ICCs
  icc_all <- icc(fit)
  expect_length(icc_all, 2)  # group_id + Person
  expect_true(all(icc_all >= 0 & icc_all <= 1))
  expect_true(sum(icc_all) <= 1.0)  # Total variance includes residual

  # Specific level
  icc_group <- icc(fit, level = "group_id")
  expect_length(icc_group, 1)
  expect_equal(icc_group, icc_all["group_id"])
})

test_that("ranef extracts random effects correctly", {
  data <- simulate_multilevel_data(n_persons = 100, n_items = 10, n_groups = 10)

  fit <- fit_irt(data$responses, model = "Rasch",
                 person_data = data$person_data,
                 random = ~ (1 | group_id), se = FALSE)

  # All random effects
  re_all <- ranef(fit)
  expect_type(re_all, "list")
  expect_equal(names(re_all), "group_id")
  expect_length(re_all$group_id, 10)

  # Specific level
  re_group <- ranef(fit, level = "group_id")
  expect_length(re_group, 10)
  expect_equal(re_group, re_all$group_id)
})

test_that("abilities extracts person abilities correctly", {
  data <- simulate_multilevel_data(n_persons = 100, n_items = 10)

  fit <- fit_irt(data$responses, model = "Rasch",
                 person_data = data$person_data,
                 random = ~ (1 | group_id), se = FALSE)

  # Person deviations only
  theta_0 <- abilities(fit, composite = FALSE)
  expect_length(theta_0, 100)

  # Composite abilities (person + group effects)
  theta_comp <- abilities(fit, composite = TRUE)
  expect_length(theta_comp, 100)

  # Composite should differ from deviations
  expect_false(all(theta_0 == theta_comp))
})

test_that("multilevel 2PL model works", {
  data <- simulate_multilevel_data(n_persons = 100, n_items = 10)

  fit <- fit_irt(data$responses, model = "2PL",
                 person_data = data$person_data,
                 random = ~ (1 | group_id), se = FALSE)

  expect_s3_class(fit, "gllamm_irt_multilevel")
  expect_true(fit$convergence$converged)
  expect_false(is.null(fit$random_effects))
})

test_that("multilevel 3PL model works", {
  data <- simulate_multilevel_data(n_persons = 100, n_items = 10)

  fit <- fit_irt(data$responses, model = "3PL",
                 person_data = data$person_data,
                 random = ~ (1 | group_id), se = FALSE)

  expect_s3_class(fit, "gllamm_irt_multilevel")
  expect_true(fit$convergence$converged)
  expect_false(is.null(fit$random_effects))
})

test_that("multilevel polytomous model works", {
  set.seed(123)
  n_persons <- 100
  n_items <- 8
  responses <- matrix(sample(1:4, n_persons * n_items, replace = TRUE),
                      n_persons, n_items)

  person_data <- data.frame(
    person_id = 1:n_persons,
    group_id = rep(1:10, each = 10)
  )

  fit <- fit_irt(responses, model = "GRM",
                 person_data = person_data,
                 random = ~ (1 | group_id), se = FALSE)

  expect_s3_class(fit, "gllamm_irt_multilevel")
  expect_true(fit$convergence$converged)
  expect_false(is.null(fit$random_effects))
})

test_that("partial nesting (NA in grouping variable) works", {
  data <- simulate_multilevel_data(n_persons = 100, n_items = 10)

  # Set some students to NA (not in any group)
  data$person_data$group_id[91:100] <- NA

  fit <- fit_irt(data$responses, model = "Rasch",
                 person_data = data$person_data,
                 random = ~ (1 | group_id), se = FALSE)

  expect_s3_class(fit, "gllamm_irt_multilevel")
  expect_true(fit$convergence$converged)

  # Should have 9 groups (original 10, but one fully NA from simulation + setting last 10 to NA)
  # The simulation creates 10 groups of 10 people each, setting last 10 to NA leaves 9 groups
  expect_equal(fit$random_effects$n_groups, 9)
})

test_that("multiple random effects work", {
  set.seed(123)
  n_persons <- 100
  n_items <- 10

  person_data <- data.frame(
    person_id = 1:n_persons,
    school_id = rep(1:5, each = 20),
    class_id = rep(1:20, each = 5)
  )

  theta <- rnorm(n_persons)
  difficulty <- rnorm(n_items)
  responses <- matrix(rbinom(n_persons * n_items, 1,
                             plogis(outer(theta, difficulty, '-'))),
                      n_persons, n_items)

  fit <- fit_irt(responses, model = "Rasch",
                 person_data = person_data,
                 random = ~ (1 | school_id) + (1 | class_id), se = FALSE)

  expect_s3_class(fit, "gllamm_irt_multilevel")
  expect_equal(length(fit$random_effects$group_names), 2)
  expect_equal(fit$random_effects$group_names, c("school_id", "class_id"))
  expect_equal(fit$random_effects$n_groups, c(5, 20))
})

test_that("print method works for multilevel models", {
  data <- simulate_multilevel_data(n_persons = 100, n_items = 10)

  fit <- fit_irt(data$responses, model = "Rasch",
                 person_data = data$person_data,
                 random = ~ (1 | group_id), se = FALSE)

  # Should not error
  expect_output(print(fit), "Multi-Level IRT Model")
  expect_output(print(fit), "Random Effects")
  expect_output(print(fit), "Variance Components")
  expect_output(print(fit), "Intraclass Correlations")
})

test_that("validation catches errors", {
  data <- simulate_multilevel_data(n_persons = 100, n_items = 10)

  # random without person_data
  expect_error(
    fit_irt(data$responses, model = "Rasch", random = ~ (1 | group_id), se = FALSE),
    "person_data must be provided"
  )

  # person_data wrong size
  wrong_data <- data$person_data[1:50, ]
  expect_error(
    fit_irt(data$responses, model = "Rasch",
            person_data = wrong_data, random = ~ (1 | group_id), se = FALSE),
    "same number of rows"
  )

  # nonexistent grouping variable
  expect_error(
    fit_irt(data$responses, model = "Rasch",
            person_data = data$person_data, random = ~ (1 | nonexistent), se = FALSE),
    "not found in person_data"
  )
})
