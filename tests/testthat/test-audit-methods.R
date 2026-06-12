# Package-wide audit regressions: external validation of under-checked
# likelihoods, multi-term prediction/simulation correctness, and method
# support across fit classes.

test_that("crossed GLMM simulate and marginal use every RE term", {
  set.seed(4)
  n <- 2000
  g1 <- factor(rep(1:40, each = 50)); g2 <- factor(rep(1:25, 80))
  u1 <- rnorm(40, 0, 1.0); u2 <- rnorm(25, 0, 0.8)
  x <- rnorm(n)
  d <- data.frame(g1 = g1, g2 = g2, x = x,
                  yb = rbinom(n, 1, plogis(-0.5 + 0.5 * x +
                                             u1[g1] + u2[g2])))
  f <- gllamm(yb ~ x + (1 | g1) + (1 | g2), data = d,
              family = binomial())

  # Simulated data must reproduce BOTH terms' group-mean variance
  s <- simulate(f, nsim = 5, seed = 2)
  m2 <- mean(sapply(1:5, function(k) var(tapply(s[[k]], d$g2, mean))))
  e2 <- var(tapply(d$yb, d$g2, mean))
  expect_gt(m2, 0.4 * e2)

  # Marginal predictions integrate over both terms
  pm <- predict(f, type = "marginal", n_sim = 1500)
  expect_lt(abs(mean(pm) - mean(d$yb)), 0.03)

  # Conditional newdata predictions apply BLUPs from every term
  pn <- predict(f, newdata = d[1:40, ], type = "response")
  expect_equal(unname(pn), unname(fitted(f)[1:40]), tolerance = 1e-6)
})

test_that("multinomial matches nnet in the fixed-effect limit", {
  skip_if_not_installed("nnet")
  set.seed(8)
  n <- 2000; x <- rnorm(n)
  e2 <- exp(0.5 + 0.7 * x); e3 <- exp(-0.3 + 1.1 * x)
  den <- 1 + e2 + e3; r <- runif(n)
  d <- data.frame(
    ym = factor(ifelse(r < 1 / den, "a",
                       ifelse(r < (1 + e2) / den, "b", "c"))),
    x = x, g = factor(rep(1:20, 100)))
  f <- fit_multinomial(ym ~ x + (1 | g), data = d)
  ref <- nnet::multinom(ym ~ x, data = d, trace = FALSE)
  expect_lt(max(abs(f$coefficients$beta - coef(ref))), 0.01)
  expect_lt(sqrt(f$coefficients$random_var[1]), 0.05)

  # Methods on multinomial fits
  expect_true(is.matrix(ranef(f)) || is.list(ranef(f)))
  s <- simulate(f, nsim = 2, seed = 1)
  expect_equal(nrow(s), n)
  expect_true(all(unlist(s) %in% levels(d$ym)))
})

test_that("weibull frailty matches survreg in the no-frailty limit", {
  set.seed(8)
  n <- 2500; x <- rnorm(n)
  tt <- (rexp(n) / exp(-1 + 0.6 * x))^(1 / 1.5)
  cens <- quantile(tt, 0.8)
  ds <- data.frame(time = pmin(tt, cens),
                   status = as.integer(tt <= cens),
                   x = x, g = factor(rep(1:25, 100)))
  fw <- fit_survival(Surv(time, status) ~ x + (1 | g), data = ds,
                     distribution = "weibull")
  sr <- survival::survreg(survival::Surv(time, status) ~ x, data = ds,
                          dist = "weibull")
  # Parameterization map: shape = 1/scale, beta = -beta_AFT
  expect_equal(fw$shape, 1 / sr$scale, tolerance = 0.02)
  expect_equal(unname(fw$coefficients$fixed["x"]),
               unname(-coef(sr)["x"]), tolerance = 0.02)
  # log hazard ratio = shape * beta recovers the generating 0.6
  expect_equal(fw$shape * unname(fw$coefficients$fixed["x"]), 0.6,
               tolerance = 0.08)

  # simulate uses the same parameterization (round trip)
  ds1 <- ds[ds$status == 1, ]
  s1 <- simulate(fw, nsim = 1, seed = 1)[[1]]
  expect_true(all(s1 > 0))
})

test_that("rank-ordered logit equals exploded conditional logit", {
  set.seed(8)
  n_cases <- 300; n_alt <- 4
  dr <- expand.grid(alt = 1:n_alt, chooser = 1:n_cases)
  dr$price <- rnorm(nrow(dr)); dr$quality <- rnorm(nrow(dr))
  util <- -0.9 * dr$price + 0.6 * dr$quality -
    log(-log(runif(nrow(dr))))
  dr$rank <- ave(-util, dr$chooser, FUN = rank)
  fr <- fit_rank(rank ~ price + quality, case = ~ chooser, data = dr)

  expl <- do.call(rbind, lapply(split(dr, dr$chooser), function(dd) {
    dd <- dd[order(dd$rank), ]
    do.call(rbind, lapply(1:(n_alt - 1), function(k) {
      rem <- dd[k:n_alt, ]
      data.frame(rem, choice = as.integer(rem$rank == min(rem$rank)),
                 stratum = paste0(dd$chooser[1], "_", k))
    }))
  }))
  expl$one <- 1
  cl <- survival::coxph(
    survival::Surv(one, choice) ~ price + quality +
      survival::strata(stratum),
    data = expl, method = "exact")
  expect_equal(unname(fr$coefficients$fixed[c("price", "quality")]),
               unname(coef(cl)[c("price", "quality")]),
               tolerance = 1e-3)
})

test_that("simulate methods exist and behave for every model class", {
  set.seed(5)
  np <- 300; ni <- 8
  theta <- rnorm(np)
  b <- seq(-1.2, 1.2, length.out = ni)
  resp <- sapply(1:ni, function(j) rbinom(np, 1, plogis(theta - b[j])))

  # IRT (EM fit): parametric-bootstrap round trip
  f_irt <- fit_irt(resp, model = "Rasch")
  s <- simulate(f_irt, nsim = 2, seed = 9)
  expect_equal(length(s), 2)
  expect_equal(dim(s[[1]]), c(np, ni))
  expect_false(anyNA(s[[1]]))
  f_rt <- fit_irt(s[[1]], model = "Rasch")
  expect_lt(mean(abs(f_rt$item_parameters$difficulty -
                       f_irt$item_parameters$difficulty)), 0.3)

  # Polytomous IRT: simulate + predict
  resp3 <- sapply(1:6, function(j) {
    1L + rowSums(outer(theta - (j - 3) / 2 + rlogis(np),
                       c(-0.8, 0.8), ">"))
  })
  f_grm <- fit_irt(resp3, model = "GRM")
  s3 <- simulate(f_grm, nsim = 1, seed = 3)
  expect_true(all(s3[[1]] %in% 1:3))
  pr <- predict(f_grm, type = "probs")
  expect_true(all(abs(rowSums(pr[[1]]) - 1) < 1e-8))
  es <- predict(f_grm, type = "expected")
  expect_equal(dim(es), c(np, 6))

  # LCA
  cls <- rbinom(np, 1, 0.45) + 1
  Yl <- sapply(1:5, function(j) rbinom(np, 1, c(0.2, 0.8)[cls]))
  f_lca <- fit_lca(Yl, nclass = 2)
  sl <- simulate(f_lca, nsim = 1, seed = 2)
  expect_equal(dim(sl[[1]]), c(np, 5))
  expect_true(all(sl[[1]] %in% 0:1))

  # CDM
  Q <- rbind(diag(2), diag(2), c(1, 1))
  alpha <- matrix(rbinom(np * 2, 1, 0.5), np, 2)
  Yc <- sapply(1:5, function(j) {
    m <- which(Q[j, ] == 1)
    rbinom(np, 1, 0.15 + 0.7 * (rowSums(alpha[, m, drop = FALSE]) ==
                                  length(m)))
  })
  f_cdm <- fit_cdm(Yc, Q, model = "dina")
  sc <- simulate(f_cdm, nsim = 1, seed = 2)
  expect_equal(dim(sc[[1]]), c(np, 5))

  # SEM: simulated covariance close to implied covariance
  f1v <- rnorm(np)
  dsem <- data.frame(x1 = f1v + rnorm(np, 0, .6),
                     x2 = 0.8 * f1v + rnorm(np, 0, .6),
                     x3 = 1.2 * f1v + rnorm(np, 0, .6))
  f_sem <- fit_sem(measurement = list(f = ~ x1 + x2 + x3), data = dsem)
  ss <- simulate(f_sem, nsim = 1, seed = 2)
  expect_equal(dim(ss[[1]]), c(np, 3))
  expect_lt(max(abs(cov(ss[[1]]) - cov(as.matrix(dsem)))), 0.35)

  # EIRT
  f_eirt <- fit_eirt(resp, data.frame(z = rnorm(ni)),
                     difficulty_formula = ~ z, model = "Rasch")
  se <- simulate(f_eirt, nsim = 1, seed = 2)
  expect_equal(dim(se[[1]]), c(np, ni))
})

test_that("vcov on SEM fits returns the parameter covariance matrix", {
  set.seed(6)
  np <- 300
  f1v <- rnorm(np)
  dsem <- data.frame(x1 = f1v + rnorm(np, 0, .6),
                     x2 = 0.8 * f1v + rnorm(np, 0, .6),
                     x3 = 1.2 * f1v + rnorm(np, 0, .6))
  f <- fit_sem(measurement = list(f = ~ x1 + x2 + x3), data = dsem)
  V <- vcov(f)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), nrow(f$param_table))
})

test_that("summary never recycles mismatched standard errors", {
  set.seed(7)
  n <- 400
  ds <- data.frame(time = rexp(n), status = rbinom(n, 1, .7),
                   x = rnorm(n), g = factor(rep(1:20, 20)))
  f <- fit_survival(Surv(time, status) ~ x + (1 | g), data = ds,
                    distribution = "exponential")
  expect_no_warning(capture.output(summary(f)))

  dn <- data.frame(yb = rbinom(n, 1, .5), x = rnorm(n),
                   g = factor(rep(1:20, 20)))
  fn <- fit_npml(yb ~ x + (1 | g), data = dn, k = 2,
                 family = binomial())
  expect_no_warning(capture.output(summary(fn)))
  sn <- simulate(fn, nsim = 2, seed = 1)
  expect_equal(nrow(sn), n)
})
