test_that("Simulation from known parameters can be recovered", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(42)

  # True parameters
  beta_0 <- 2.0
  beta_1 <- -0.5
  sigma_u <- 0.8
  sigma <- 1.2

  # Generate data
  n_groups <- 50
  n_per_group <- 20
  n_obs <- n_groups * n_per_group

  group <- rep(1:n_groups, each = n_per_group)
  x <- rnorm(n_obs, mean = 0, sd = 1)

  # Random effects
  u <- rnorm(n_groups, mean = 0, sd = sigma_u)
  u_expanded <- u[group]

  # Response
  eta <- beta_0 + beta_1 * x + u_expanded
  y <- rnorm(n_obs, mean = eta, sd = sigma)

  data <- data.frame(y = y, x = x, group = group)

  # Fit model
  fit <- gllamm(y ~ x + (1 | group), data = data)

  # Check recovery (with generous tolerances for finite sample)
  expect_equal(fixef(fit)[["(Intercept)"]], beta_0, tolerance = 0.15)
  expect_equal(fixef(fit)[["x"]], beta_1, tolerance = 0.15)
  expect_equal(sqrt(VarCorr(fit)[[1]]), sigma_u, tolerance = 0.20)
})


test_that("Small sample simulation recovery", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(99)

  # True parameters
  beta_0 <- 1.0
  sigma_u <- 0.5

  # Small sample
  n_groups <- 10
  n_per_group <- 5

  group <- rep(1:n_groups, each = n_per_group)
  u <- rnorm(n_groups, mean = 0, sd = sigma_u)
  y <- rnorm(n_groups * n_per_group, mean = beta_0 + u[group], sd = 1)

  data <- data.frame(y = y, group = group)

  # Fit model
  fit <- gllamm(y ~ 1 + (1 | group), data = data)

  # Check that estimates are in reasonable range
  expect_true(abs(fixef(fit)[1] - beta_0) < 0.5)
  expect_true(sqrt(VarCorr(fit)[[1]]) > 0)
})


test_that("simulate.gllamm produces correct dimensions", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(123)
  data <- data.frame(
    y = rnorm(50),
    x = rnorm(50),
    group = rep(1:10, each = 5)
  )

  fit <- gllamm(y ~ x + (1 | group), data = data)

  # Single simulation
  sim1 <- simulate(fit, nsim = 1, seed = 456)
  expect_equal(length(sim1), 50)
  expect_type(sim1, "double")

  # Multiple simulations
  sim5 <- simulate(fit, nsim = 5, seed = 456)
  expect_equal(dim(sim5), c(50, 5))
})


test_that("Simulated data has correct properties", {
  skip_if_not_installed("TMB")
  skip("TMB compilation required")

  set.seed(789)

  # Simple model with known structure
  n_groups <- 20
  n_per_group <- 10
  data <- data.frame(
    y = rnorm(n_groups * n_per_group),
    group = rep(1:n_groups, each = n_per_group)
  )

  fit <- gllamm(y ~ 1 + (1 | group), data = data)

  # Simulate many datasets
  sims <- simulate(fit, nsim = 100, seed = 321)

  # Check mean is approximately correct (should be near intercept)
  sim_means <- apply(sims, 2, mean)
  expect_equal(mean(sim_means), fixef(fit)[1], tolerance = 0.5)

  # Check variance structure is reasonable
  sim_vars <- apply(sims, 2, var)
  expect_true(all(sim_vars > 0))
})
