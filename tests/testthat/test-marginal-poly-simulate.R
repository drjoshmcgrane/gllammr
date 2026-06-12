# Stage B3/B4/B5 features: polytomous marginal predictions, the shared
# category-probability helper, simulate() across families, Cook's distance.

simulate_poly_data <- function(seed = 42, np = 200, ni = 6) {
  set.seed(seed)
  theta <- rnorm(np)
  resp <- sapply(seq(-1, 1, length.out = ni), function(b) {
    p1 <- plogis(theta - b - 0.8)
    p2 <- plogis(theta - b + 0.8)
    1 + (runif(np) < p2) + (runif(np) < p1)
  })
  list(resp = resp, theta = theta)
}

simulate_glmm_data <- function(seed = 42, g = 25, n_per = 20) {
  set.seed(seed)
  n <- g * n_per
  grp <- factor(rep(1:g, each = n_per))
  x <- rnorm(n)
  u <- rnorm(g, 0, 1)
  data.frame(
    x = x, grp = grp,
    y  = 1 + 0.5 * x + u[as.integer(grp)] + rnorm(n),
    yb = rbinom(n, 1, plogis(0.5 * x + u[as.integer(grp)])),
    yp = rpois(n, exp(0.2 + 0.3 * x + 0.5 * u[as.integer(grp)])),
    yo = cut(0.8 * x + u[as.integer(grp)] + rlogis(n),
             c(-Inf, -1, 0.7, Inf), labels = FALSE)
  )
}

test_that("irt_category_probs matches the template orientation for GRM", {
  # P(lowest category) must DECREASE with theta; P(highest) must INCREASE
  p <- gllammr:::irt_category_probs("GRM", c(-3, 0, 3),
                                    thresholds = c(-1, 1), discrimination = 1)
  expect_equal(rowSums(p), rep(1, 3), tolerance = 1e-10)
  expect_true(p[1, 1] > p[3, 1])   # low ability -> low category
  expect_true(p[3, 3] > p[1, 3])   # high ability -> high category
})

test_that("irt_category_probs rows sum to one for all models", {
  theta <- seq(-3, 3, length.out = 11)
  for (m in c("GRM", "PCM", "GPCM", "NRM")) {
    p <- gllammr:::irt_category_probs(m, theta,
                                      thresholds = c(-0.5, 0.5),
                                      discrimination = 1.2)
    expect_equal(rowSums(p), rep(1, 11), tolerance = 1e-10)
    expect_true(all(p >= 0))
  }
})

test_that("polytomous IRT marginal predictions return category matrix", {
  s <- simulate_poly_data()
  fit <- fit_irt(s$resp, model = "GRM")
  m <- predict(fit, type = "marginal")

  expect_true(is.matrix(m))
  expect_equal(nrow(m), ncol(s$resp))
  expect_equal(unname(rowSums(m, na.rm = TRUE)), rep(1, nrow(m)), tolerance = 0.01)
})

test_that("polytomous EIRT marginal and conditional predictions work", {
  s <- simulate_poly_data()
  idata <- data.frame(x = seq(-1, 1, length.out = ncol(s$resp)))
  fit <- fit_eirt(s$resp, item_data = idata, difficulty_formula = ~ x,
                  model = "PCM")

  m <- predict(fit, type = "marginal")
  expect_equal(unname(rowSums(m, na.rm = TRUE)), rep(1, nrow(m)), tolerance = 0.01)

  pr <- predict(fit, type = "probability", ability = c(-1, 0, 1))
  expect_type(pr, "list")
  expect_length(pr, ncol(s$resp))
  expect_equal(rowSums(pr[[1]]), rep(1, 3), tolerance = 1e-8)

  # New items: PCM uses predicted location with steps at population zero
  mn <- predict(fit, newdata = data.frame(x = c(-0.5, 0.5)), type = "marginal")
  expect_equal(nrow(mn), 2)
})

test_that("GRM marginal predictions for new items error informatively", {
  s <- simulate_poly_data()
  idata <- data.frame(x = seq(-1, 1, length.out = ncol(s$resp)))
  fit <- fit_eirt(s$resp, item_data = idata, difficulty_formula = ~ x,
                  model = "GRM")
  expect_error(predict(fit, newdata = data.frame(x = 0), type = "marginal"),
               "item-specific")
})

test_that("simulate() follows the stats contract across families", {
  d <- simulate_glmm_data()

  fits <- list(
    gaussian = gllamm(y ~ x + (1 | grp), data = d),
    binomial = gllamm(yb ~ x + (1 | grp), data = d, family = stats::binomial()),
    poisson  = gllamm(yp ~ x + (1 | grp), data = d, family = stats::poisson())
  )

  for (nm in names(fits)) {
    s <- simulate(fits[[nm]], nsim = 3, seed = 7)
    expect_s3_class(s, "data.frame")
    expect_equal(dim(s), c(nrow(d), 3L))
    expect_false(is.null(attr(s, "seed")))
  }

  expect_true(all(unlist(simulate(fits$binomial, nsim = 2, seed = 1)) %in% 0:1))
  expect_true(all(unlist(simulate(fits$poisson, nsim = 2, seed = 1)) >= 0))

  # newdata simulation
  sn <- simulate(fits$gaussian, nsim = 2, newdata = d[1:100, ], seed = 2)
  expect_equal(nrow(sn), 100)
})

test_that("simulate() reproduces with the same seed", {
  d <- simulate_glmm_data()
  fit <- gllamm(y ~ x + (1 | grp), data = d)
  s1 <- simulate(fit, nsim = 2, seed = 99)
  s2 <- simulate(fit, nsim = 2, seed = 99)
  expect_equal(s1, s2)
})

test_that("simulate.gllamm_ordinal produces valid categories", {
  d <- simulate_glmm_data()
  fit <- fit_ordinal(yo ~ x + (1 | grp), data = d, link = "logit")
  s <- simulate(fit, nsim = 2, seed = 5)
  expect_s3_class(s, "data.frame")
  expect_true(all(unlist(s) %in% 1:3))
})

test_that("cooks.distance.gllamm flags an outlying cluster", {
  d <- simulate_glmm_data(g = 15, n_per = 15)
  # Contaminate one cluster heavily
  idx <- d$grp == "3"
  d$y[idx] <- d$y[idx] + 8

  fit <- gllamm(y ~ x + (1 | grp), data = d)
  D <- cooks.distance(fit)

  expect_length(D, 15)
  expect_true(all(is.finite(D)))
  expect_equal(unname(which.max(D)), 3)
})

test_that("cooks.distance respects max_groups guard", {
  d <- simulate_glmm_data(g = 25, n_per = 8)
  fit <- gllamm(y ~ x + (1 | grp), data = d)
  expect_error(cooks.distance(fit, max_groups = 10), "max_groups")
})

test_that("polytomous ICC and DIF plots draw without error", {
  s <- simulate_poly_data()
  fit_grm <- fit_irt(s$resp, model = "GRM")
  fit_pcm <- fit_irt(s$resp, model = "PCM")

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(fit_grm, which = 1))
  expect_no_error(plot(fit_pcm, which = 1))
})
