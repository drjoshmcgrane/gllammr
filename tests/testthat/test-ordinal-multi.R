# Multi-term (crossed/nested) random effects for ordinal and multinomial

sim_crossed_ordinal <- function(seed = 77, n_raters = 50, n_items = 30) {
  set.seed(seed)
  d <- expand.grid(rater = factor(1:n_raters), item = factor(1:n_items))
  u_r <- rnorm(n_raters, 0, 1.0)
  u_i <- rnorm(n_items, 0, 0.7)
  d$x <- rnorm(nrow(d))
  eta <- 0.8 * d$x + u_r[d$rater] + u_i[d$item]
  cuts <- c(-1.5, 0, 1.5)
  d$rating <- factor(1L + rowSums(outer(eta + rlogis(nrow(d)), cuts, ">")),
                     levels = 1:4, ordered = TRUE)
  d$rating_num <- as.integer(d$rating)
  d
}

test_that("crossed-RE ordinal matches ordinal::clmm exactly", {
  skip_if_not_installed("ordinal")
  skip_on_cran()  # cross-package agreement is CI-only
  d <- sim_crossed_ordinal()
  fit <- fit_ordinal(rating_num ~ x + (1 | rater) + (1 | item), data = d)
  ref <- ordinal::clmm(rating ~ x + (1 | rater) + (1 | item), data = d,
                       link = "logit")

  expect_equal(fit$logLik, as.numeric(logLik(ref)), tolerance = 1e-4)
  expect_equal(unname(fit$coefficients$fixed["x"]),
               unname(coef(ref)["x"]), tolerance = 1e-3)
  expect_equal(unname(sqrt(fit$coefficients$random_var$rater[1, 1])),
               unname(attr(ordinal::VarCorr(ref)$rater, "stddev")),
               tolerance = 1e-3)
  expect_equal(unname(sqrt(fit$coefficients$random_var$item[1, 1])),
               unname(attr(ordinal::VarCorr(ref)$item, "stddev")),
               tolerance = 1e-3)
  expect_equal(unname(fit$coefficients$thresholds),
               unname(ref$alpha), tolerance = 1e-3)
})

test_that("crossed ordinal works on wine through gllamm() and probit", {
  skip_if_not_installed("ordinal")
  skip_on_cran()  # cross-package agreement is CI-only
  data("wine", package = "ordinal", envir = environment())
  wine$rating_num <- as.integer(wine$rating)

  fit <- gllamm(rating_num ~ temp + (1 | judge) + (1 | bottle),
                data = wine, family = ordinal())
  ref <- ordinal::clmm(rating ~ temp + (1 | judge) + (1 | bottle),
                       data = wine, link = "logit")
  expect_equal(fit$logLik, as.numeric(logLik(ref)), tolerance = 1e-3)

  fit_p <- fit_ordinal(rating_num ~ temp + (1 | judge) + (1 | bottle),
                       data = wine, link = "probit")
  ref_p <- ordinal::clmm(rating ~ temp + (1 | judge) + (1 | bottle),
                         data = wine, link = "probit")
  expect_equal(fit_p$logLik, as.numeric(logLik(ref_p)), tolerance = 1e-3)
})

test_that("multi-term ordinal nests the single-term model", {
  d <- sim_crossed_ordinal(seed = 5, n_raters = 30, n_items = 20)
  fit2 <- fit_ordinal(rating_num ~ x + (1 | rater) + (1 | item), data = d)
  fit1 <- fit_ordinal(rating_num ~ x + (1 | rater), data = d)
  expect_gte(fit2$logLik, fit1$logLik - 1e-6)
  expect_equal(fit2$n_random_terms, 2)
  expect_named(fit2$coefficients$random_var, c("rater", "item"))
})

test_that("PPO link still requires a single term, with a clear error", {
  d <- sim_crossed_ordinal(seed = 6, n_raters = 20, n_items = 10)
  expect_error(
    fit_ordinal(rating_num ~ x + (1 | rater) + (1 | item), data = d,
                link = "ppo"),
    "single random")
})

test_that("crossed-RE multinomial recovers variance components", {
  set.seed(88)
  n_g1 <- 50; n_per <- 40
  d <- data.frame(g1 = factor(rep(1:n_g1, each = n_per)),
                  g2 = factor(rep(rep(1:20, length.out = n_per), n_g1)))
  u1 <- rnorm(n_g1, 0, 0.9); u2 <- rnorm(20, 0, 0.5)
  d$x <- rnorm(nrow(d))
  shift <- u1[d$g1] + u2[d$g2]
  e2 <- exp(0.5 + 0.7 * d$x + shift)
  e3 <- exp(-0.3 + 1.1 * d$x + shift)
  den <- 1 + e2 + e3
  r <- runif(nrow(d))
  d$choice <- factor(ifelse(r < 1 / den, "a",
                            ifelse(r < (1 + e2) / den, "b", "c")))

  fit <- fit_multinomial(choice ~ x + (1 | g1) + (1 | g2), data = d)
  expect_true(fit$convergence$converged)
  expect_equal(unname(sqrt(fit$coefficients$random_var$g1[1, 1])), 0.9,
               tolerance = 0.2)
  expect_equal(unname(sqrt(fit$coefficients$random_var$g2[1, 1])), 0.5,
               tolerance = 0.25)
  expect_equal(unname(fit$coefficients$beta["b", "x"]), 0.7,
               tolerance = 0.15)

  # Nests the single-term model
  fit1 <- fit_multinomial(choice ~ x + (1 | g1), data = d)
  expect_gte(fit$logLik, fit1$logLik - 1e-6)
})

test_that("cdm() and lca(ordering) are reachable through gllamm()", {
  set.seed(3)
  Q <- rbind(diag(2), diag(2), c(1, 1))
  alpha <- matrix(rbinom(300 * 2, 1, 0.5), 300, 2)
  Y <- sapply(1:5, function(j) {
    m <- which(Q[j, ] == 1)
    rbinom(300, 1, 0.15 + 0.7 * (rowSums(alpha[, m, drop = FALSE]) ==
                                   length(m)))
  })
  f_cdm <- gllamm(Y, family = cdm(Q, model = "dina"))
  expect_s3_class(f_cdm, "gllamm_cdm")
  expect_equal(f_cdm$logLik, fit_cdm(Y, Q, model = "dina")$logLik,
               tolerance = 0.5)

  Yl <- sapply(1:4, function(j) rbinom(300, 1, c(0.2, 0.8)[1 + alpha[, 1]]))
  f_lca <- gllamm(Yl, family = lca(nclass = 2, ordering = "increasing"))
  expect_s3_class(f_lca, "gllamm_lca")
  expect_true(all(apply(f_lca$item_probs, 1, diff) >= -1e-8))
})
