# Every model class is reachable through gllamm(): each unified-interface
# route must reproduce the direct fit_* call exactly.

test_that("eirt() family routes through gllamm()", {
  set.seed(1)
  np <- 300; ni <- 8
  theta <- rnorm(np)
  hard <- rep(c(0, 1), length.out = ni)
  b <- -0.5 + 0.8 * hard + rnorm(ni, 0, 0.2)
  resp <- sapply(1:ni, function(j) rbinom(np, 1, plogis(theta - b[j])))
  idata <- data.frame(hard = hard)

  set.seed(2)
  f1 <- gllamm(resp, family = eirt(idata, difficulty_formula = ~ hard,
                                   model = "Rasch"))
  set.seed(2)
  f2 <- fit_eirt(resp, idata, difficulty_formula = ~ hard,
                 model = "Rasch")
  expect_s3_class(f1, "gllamm_eirt")
  expect_equal(f1$logLik, f2$logLik, tolerance = 1e-6)
})

test_that("eirt() family routes step_formula/step_data through gllamm()", {
  set.seed(9)
  np <- 300; ni <- 8; K <- 3
  theta <- rnorm(np)
  z_item <- rnorm(ni)
  x_step <- rnorm(ni * (K - 1))
  b <- 0.2 + 0.5 * z_item
  resp <- sapply(seq_len(ni), function(j) {
    delta <- b[j] + 0.6 * x_step[((j - 1) * (K - 1) + 1):(j * (K - 1))] +
      c(-0.7, 0.7)
    sapply(seq_len(np), function(p) {
      cs <- cumsum(c(0, theta[p] - delta))
      sample.int(K, 1, prob = exp(cs) / sum(exp(cs)))
    })
  })
  idata <- data.frame(z_item = z_item)
  step_data <- data.frame(x = x_step)

  set.seed(2)
  f1 <- gllamm(resp, family = eirt(idata, difficulty_formula = ~ z_item,
                                   step_formula = ~ x, step_data = step_data,
                                   model = "PCM"))
  set.seed(2)
  f2 <- fit_eirt(resp, idata, difficulty_formula = ~ z_item,
                 step_formula = ~ x, step_data = step_data, model = "PCM")
  expect_s3_class(f1, "gllamm_eirt")
  expect_equal(f1$logLik, f2$logLik, tolerance = 1e-6)
  expect_equal(f1$regression_coefficients$step, f2$regression_coefficients$step)
})

test_that("sem() family routes through gllamm()", {
  set.seed(3)
  n <- 300
  fl <- rnorm(n)
  d <- data.frame(x1 = fl + rnorm(n, 0, .6), x2 = 0.8 * fl + rnorm(n, 0, .6),
                  x3 = 1.2 * fl + rnorm(n, 0, .6))
  f1 <- gllamm(d, family = sem(measurement = list(f1 = ~ x1 + x2 + x3)))
  f2 <- fit_sem(measurement = list(f1 = ~ x1 + x2 + x3), data = d)
  expect_s3_class(f1, "gllamm_sem")
  expect_equal(f1$logLik, f2$logLik, tolerance = 1e-6)
})

test_that("mixed_response() family routes through gllamm()", {
  set.seed(4)
  n <- 300
  g <- factor(rep(1:30, 10))
  u <- rnorm(30, 0, 0.7)
  d <- data.frame(g = g,
                  y1 = 1 + u[g] + rnorm(n),
                  y2 = rbinom(n, 1, plogis(u[g])))
  f1 <- gllamm(~ 1 | g, data = d,
               family = mixed_response(gaussian = y1 ~ 1,
                                       binomial = y2 ~ 1))
  f2 <- fit_mixed(list(gaussian = y1 ~ 1, binomial = y2 ~ 1),
                  random = ~ 1 | g, data = d)
  expect_equal(f1$logLik, f2$logLik, tolerance = 1e-6)
})

test_that("survival_family() routes through gllamm()", {
  set.seed(5)
  n <- 400
  g <- factor(rep(1:40, 10))
  u <- rnorm(40, 0, 0.5)
  x <- rnorm(n)
  t_true <- rexp(n, rate = exp(-1 + 0.5 * x + u[g]))
  cens <- rexp(n, 0.2)
  d <- data.frame(time = pmin(t_true, cens),
                  status = as.integer(t_true <= cens), x = x, g = g)
  f1 <- gllamm(Surv(time, status) ~ x + (1 | g), data = d,
               family = survival_family("exponential"))
  f2 <- fit_survival(Surv(time, status) ~ x + (1 | g), data = d,
                     distribution = "exponential")
  expect_equal(f1$logLik, f2$logLik, tolerance = 1e-6)
})

test_that("ranking() family routes through gllamm()", {
  set.seed(6)
  n_cases <- 80; n_alt <- 4
  d <- expand.grid(alt = 1:n_alt, chooser = 1:n_cases)
  d$price <- rnorm(nrow(d))
  d$rank <- as.vector(replicate(n_cases, sample(1:n_alt)))
  f1 <- gllamm(rank ~ price, data = d, family = ranking(case = ~ chooser))
  f2 <- fit_rank(rank ~ price, case = ~ chooser, data = d)
  expect_equal(f1$logLik, f2$logLik, tolerance = 1e-6)
})

test_that("integration = npml(k) routes through gllamm()", {
  set.seed(7)
  g <- 60; n_per <- 8
  grp <- factor(rep(1:g, each = n_per))
  cls <- sample(1:2, g, TRUE)
  locs <- c(-1, 1)
  d <- data.frame(grp = grp,
                  yb = rbinom(g * n_per, 1, plogis(locs[cls[grp]])))
  f1 <- gllamm(yb ~ 1 + (1 | grp), data = d, family = binomial(),
               integration = npml(2))
  f2 <- fit_npml(yb ~ 1 + (1 | grp), data = d, k = 2,
                 family = stats::binomial())
  expect_equal(f1$logLik, f2$logLik, tolerance = 1e-4)
})

test_that("binomial() with integration = aghq(k) is honored", {
  set.seed(8)
  g <- 50; n_per <- 6
  grp <- factor(rep(1:g, each = n_per))
  u <- rnorm(g, 0, 1.5)
  d <- data.frame(grp = grp,
                  yb = rbinom(g * n_per, 1, plogis(u[grp])))
  f_aghq <- gllamm(yb ~ 1 + (1 | grp), data = d, family = binomial(),
                   integration = aghq(15))
  f_lap <- gllamm(yb ~ 1 + (1 | grp), data = d, family = binomial())
  # AGHQ refines the Laplace approximation: likelihoods must differ
  # (previously the integration argument was silently ignored)
  expect_false(isTRUE(all.equal(f_aghq$logLik, f_lap$logLik,
                                tolerance = 1e-8)))
  if (requireNamespace("lme4", quietly = TRUE)) {
    ref <- ref_fit(lme4::glmer(yb ~ 1 + (1 | grp), data = d,
                               family = stats::binomial(), nAGQ = 15,
                               control = lme4::glmerControl(optimizer = "bobyqa")))
    expect_equal(f_aghq$logLik, as.numeric(logLik(ref)), tolerance = 0.05)
  }
})
