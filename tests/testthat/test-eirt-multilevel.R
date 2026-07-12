# Multi-level explanatory IRT: person-level group random effects combined
# with item predictors, for dichotomous (Rasch/2PL) and polytomous (PCM) models.

simulate_ml_eirt <- function(seed = 99, n_class = 20, per_class = 25, ni = 12,
                             sigma_class = 0.7) {
  set.seed(seed)
  np <- n_class * per_class
  class_id <- rep(1:n_class, each = per_class)
  u_class <- rnorm(n_class, 0, sigma_class)
  theta <- rnorm(np, 0, 1) + u_class[class_id]
  x_item <- runif(ni, -1, 1)
  b <- 0.4 + 1.2 * x_item + rnorm(ni, 0, 0.3)
  resp <- matrix(rbinom(np * ni, 1, plogis(outer(theta, b, "-"))), np, ni)
  list(resp = resp, theta = theta, b = b,
       item_data = data.frame(x = x_item),
       person_data = data.frame(class = factor(class_id)))
}

test_that("multilevel Rasch-EIRT recovers item regression and variance components", {
  s <- simulate_ml_eirt()
  fit <- fit_eirt(s$resp, item_data = s$item_data, difficulty_formula = ~ x,
                  person_data = s$person_data, random = ~ (1 | class),
                  model = "Rasch")

  expect_true(fit$convergence$converged)
  expect_s3_class(fit, "gllamm_eirt_multilevel")

  gamma <- fit$regression_coefficients$difficulty
  expect_equal(unname(gamma[2]), 1.2, tolerance = 0.15)   # item predictor slope

  re <- fit$random_effects
  expect_equal(unname(re$sigma_random), 0.7, tolerance = 0.25)
  expect_equal(unname(fit$ability_sd), 1.0, tolerance = 0.15)
  expect_length(re$composite_theta, nrow(s$resp))
  expect_named(re$icc, c("class", "Person"))
})

test_that("multilevel fit beats single-level on clustered data", {
  skip_on_cran()  # two multilevel EIRT fits; multilevel smoke kept above
  s <- simulate_ml_eirt()
  fit_ml <- fit_eirt(s$resp, item_data = s$item_data, difficulty_formula = ~ x,
                     person_data = s$person_data, random = ~ (1 | class),
                     model = "Rasch")
  fit_sl <- fit_eirt(s$resp, item_data = s$item_data, difficulty_formula = ~ x,
                     model = "Rasch")

  expect_gt(fit_ml$logLik, fit_sl$logLik)
  # Single-level total ability SD should absorb the class variance
  expect_equal(unname(fit_sl$ability_sd), sqrt(1 + 0.7^2), tolerance = 0.15)
  expect_s3_class(fit_sl, "gllamm_eirt")
  expect_false(inherits(fit_sl, "gllamm_eirt_multilevel"))
})

test_that("single-level results are unchanged when random = NULL", {
  s <- simulate_ml_eirt(seed = 5, n_class = 10, per_class = 20, ni = 8)
  f1 <- fit_eirt(s$resp, item_data = s$item_data, difficulty_formula = ~ x,
                 model = "Rasch")
  f2 <- fit_eirt(s$resp, item_data = s$item_data, difficulty_formula = ~ x,
                 person_data = s$person_data, random = NULL, model = "Rasch")
  expect_equal(f1$logLik, f2$logLik, tolerance = 1e-6)
  expect_equal(f1$regression_coefficients$difficulty,
               f2$regression_coefficients$difficulty, tolerance = 1e-6)
})

test_that("multilevel polytomous (PCM) EIRT fits and recovers structure", {
  skip_on_cran()  # multilevel polytomous EIRT fit; CI-only
  s <- simulate_ml_eirt()
  np <- nrow(s$resp); ni <- ncol(s$resp)
  set.seed(100)
  delta1 <- s$b - 0.5; delta2 <- s$b + 0.5
  respP <- matrix(0L, np, ni)
  for (i in 1:np) for (j in 1:ni) {
    num <- c(0, s$theta[i] - delta1[j],
             (s$theta[i] - delta1[j]) + (s$theta[i] - delta2[j]))
    p <- exp(num - max(num)); p <- p / sum(p)
    respP[i, j] <- sample(1:3, 1, prob = p)
  }

  fit <- fit_eirt(respP, item_data = s$item_data, difficulty_formula = ~ x,
                  person_data = s$person_data, random = ~ (1 | class),
                  model = "PCM")

  expect_true(fit$convergence$converged)
  expect_s3_class(fit, "gllamm_eirt_multilevel")
  gamma <- fit$regression_coefficients$difficulty
  expect_equal(unname(gamma[2]), 1.2, tolerance = 0.2)
  expect_equal(unname(fit$random_effects$sigma_random), 0.7, tolerance = 0.25)
})

test_that("random without person_data errors clearly", {
  s <- simulate_ml_eirt(seed = 5, n_class = 6, per_class = 10, ni = 6)
  expect_error(
    fit_eirt(s$resp, item_data = s$item_data, random = ~ (1 | class),
             model = "Rasch"),
    "person_data"
  )
})
