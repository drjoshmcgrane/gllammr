# E5: rank-ordered (exploded) logit.

simulate_ranks <- function(seed = 111, n_case = 400, n_alt = 4, n_grp = 25,
                           sigma_u = 0.6) {
  set.seed(seed)
  grp <- sample(1:n_grp, n_case, TRUE)
  u <- rnorm(n_grp, 0, sigma_u)
  rows <- n_case * n_alt
  d <- data.frame(
    id = rep(1:n_case, each = n_alt),
    region = factor(rep(grp, each = n_alt)),
    price = runif(rows, 0, 2),
    quality = rnorm(rows)
  )
  gumbel <- -log(-log(runif(rows)))
  util <- -1.0 * d$price + 0.8 * d$quality + u[grp[d$id]] * d$price + gumbel
  d$rank <- ave(-util, d$id, FUN = rank)
  d
}

test_that("rank-ordered logit recovers preferences and taste heterogeneity", {
  d <- simulate_ranks()
  fit <- fit_rank(rank ~ price + quality, case = ~ id,
                  random = ~ (0 + price | region), data = d)

  expect_true(fit$convergence$converged)
  expect_equal(unname(coef(fit)$fixed["price"]), -1.0, tolerance = 0.2)
  expect_equal(unname(coef(fit)$fixed["quality"]), 0.8, tolerance = 0.15)
  expect_equal(unname(fit$coefficients$random_sd), 0.6, tolerance = 0.3)
})

test_that("fixed-effects rank logit and partial rankings work", {
  d <- simulate_ranks(seed = 112)
  fit0 <- fit_rank(rank ~ price + quality, case = ~ id, data = d)
  expect_true(fit0$convergence$converged)
  expect_true(is.na(fit0$coefficients$random_sd))

  d$rank[d$rank > 2] <- NA   # top-2 rankings only
  fitp <- fit_rank(rank ~ price + quality, case = ~ id,
                   random = ~ (0 + price | region), data = d)
  expect_true(fitp$convergence$converged)
  expect_equal(unname(coef(fitp)$fixed["price"]), -1.0, tolerance = 0.25)
})

test_that("rank input validation errors clearly", {
  d <- simulate_ranks(n_case = 20)
  d$rank[d$id == 1] <- c(1, 1, 2, 3)   # tie
  expect_error(fit_rank(rank ~ price, case = ~ id, data = d), "ties")
})
