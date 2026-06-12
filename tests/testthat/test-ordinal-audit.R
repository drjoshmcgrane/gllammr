# Ordinal GLMM audit: link likelihood correctness (vs VGAM), method
# support on multi-term fits, latent-scale ICC

sim_link_data <- function(link, n = 3000, seed = 3, beta = 0.9,
                          th = c(-1.2, 0, 1.2)) {
  set.seed(seed)
  x <- rnorm(n)
  eta <- beta * x
  pr <- switch(link,
    crl_backward = {
      p <- matrix(0, n, 4); surv <- rep(1, n)
      for (c in 4:2) {
        b <- plogis(th[c - 1] - eta); p[, c] <- b * surv
        surv <- surv * (1 - b)
      }
      p[, 1] <- surv; p
    },
    acl = {
      lp <- matrix(0, n, 4)
      for (k in 2:4) lp[, k] <- lp[, k - 1] + th[k - 1] + eta
      el <- exp(lp - apply(lp, 1, max)); el / rowSums(el)
    },
    crl_forward = {
      p <- matrix(0, n, 4); surv <- rep(1, n)
      for (k in 1:3) {
        h <- plogis(th[k] - eta); p[, k] <- surv * h
        surv <- surv * (1 - h)
      }
      p[, 4] <- surv; p
    })
  y <- vapply(seq_len(n), function(i) {
    sample.int(4, 1, prob = pr[i, ])
  }, 1L)
  data.frame(y = y, x = x, g = factor(rep(1:20, length.out = n)))
}

test_that("crl_backward is a proper distribution and matches VGAM", {
  skip_if_not_installed("VGAM")
  d <- sim_link_data("crl_backward")
  f <- fit_ordinal(y ~ x + (1 | g), data = d, link = "crl_backward")
  v <- VGAM::vglm(ordered(y) ~ x, data = d,
                  family = VGAM::sratio(reverse = TRUE, parallel = TRUE))
  # Fixed effect agrees with the no-RE VGAM fit (our RE variance ~ 0)
  expect_equal(unname(f$coefficients$fixed["x"]),
               unname(-coef(v)[["x"]]), tolerance = 0.03)
  # Category probabilities sum to one (the previous likelihood summed
  # to 2 * plogis(tau_max - eta))
  pr <- predict(f, type = "probs")
  expect_true(all(abs(rowSums(pr) - 1) < 1e-10))
})

test_that("acl and crl_forward match VGAM", {
  skip_if_not_installed("VGAM")
  d1 <- sim_link_data("acl", seed = 5)
  f1 <- fit_ordinal(y ~ x + (1 | g), data = d1, link = "acl")
  v1 <- VGAM::vglm(ordered(y) ~ x, data = d1,
                   family = VGAM::acat(parallel = TRUE))
  expect_equal(unname(f1$coefficients$fixed["x"]),
               unname(coef(v1)[["x"]]), tolerance = 0.03)

  d2 <- sim_link_data("crl_forward", seed = 7)
  f2 <- fit_ordinal(y ~ x + (1 | g), data = d2, link = "crl_forward")
  v2 <- VGAM::vglm(ordered(y) ~ x, data = d2,
                   family = VGAM::sratio(reverse = FALSE, parallel = TRUE))
  expect_equal(unname(f2$coefficients$fixed["x"]),
               unname(-coef(v2)[["x"]]), tolerance = 0.03)
})

test_that("predict supports every link and sums to one", {
  set.seed(11)
  d <- data.frame(x = rnorm(600), g = factor(rep(1:20, 30)))
  d$y <- 1L + rowSums(outer(0.8 * d$x + rlogis(600), c(-1, 0, 1), ">"))
  for (lk in c("logit", "probit", "acl", "crl_forward", "crl_backward",
               "ppo")) {
    f <- fit_ordinal(y ~ x + (1 | g), data = d, link = lk)
    pr <- predict(f, type = "probs")
    expect_true(all(abs(rowSums(pr) - 1) < 1e-8),
                info = paste("link", lk))
    expect_true(all(pr >= 0), info = paste("link", lk))
  }
})

test_that("multi-term ordinal: simulate, icc, marginal all work", {
  set.seed(1)
  n_g <- 50; n <- 2500
  d <- data.frame(g = factor(rep(1:n_g, each = 50)),
                  g2 = factor(rep(1:25, 100)), x = rnorm(n))
  u <- rnorm(n_g, 0, 1); u2 <- rnorm(25, 0, 0.5)
  eta <- 0.8 * d$x + u[d$g] + u2[d$g2]
  d$y <- 1L + rowSums(outer(eta + rlogis(n), c(-1.5, 0, 1.5), ">"))
  fm <- fit_ordinal(y ~ x + (1 | g) + (1 | g2), data = d)

  s1 <- simulate(fm, nsim = 2, seed = 7)
  expect_equal(dim(s1), c(n, 2))
  expect_true(all(unlist(s1) %in% 1:4))

  i1 <- icc(fm, quiet = TRUE)
  expect_equal(length(i1), 2)
  expect_named(i1, c("g", "g2"))
  # Latent-scale: var / (sum var + pi^2/3)
  v <- sapply(fm$coefficients$random_var, function(m) m[1, 1])
  expect_equal(unname(i1), unname(v / (sum(v) + pi^2 / 3)),
               tolerance = 1e-8)

  # Marginal predictions integrate over BOTH terms: category means match
  # the empirical distribution
  pm <- predict(fm, type = "marginal", n_sim = 2000)
  expect_lt(max(abs(colMeans(pm) - as.numeric(table(d$y) / n))), 0.02)
})

test_that("single-term ordinal icc uses the latent logistic scale", {
  set.seed(13)
  n_g <- 80; n <- 4000
  g <- factor(rep(1:n_g, each = 50))
  u <- rnorm(n_g, 0, 1)                      # true ICC = 1/(1+pi^2/3) = .233
  d <- data.frame(g = g, x = rnorm(n))
  d$y <- 1L + rowSums(outer(0.5 * d$x + u[g] + rlogis(n),
                            c(-1, 0, 1), ">"))
  f <- fit_ordinal(y ~ x + (1 | g), data = d)
  i <- icc(f, quiet = TRUE)
  expect_equal(unname(i[1]), 1 / (1 + pi^2 / 3), tolerance = 0.05)
})
