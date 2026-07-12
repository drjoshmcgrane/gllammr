# Confirmatory model-based DIF (IRT-LR with latent-regression impact)

sim_dif_impact <- function(n = 500, ni = 8, seed = 55, dif_item = 4,
                           dif_size = 0.8, impact = -0.5,
                           two_pl = FALSE, kappa = 0) {
  set.seed(seed)
  g <- rep(c(0, 1), length.out = n)
  theta <- rnorm(n, mean = impact * g)
  b <- seq(-1.5, 1.5, length.out = ni)
  a <- if (two_pl) runif(ni, 0.7, 1.6) else rep(1, ni)
  resp <- sapply(seq_len(ni), function(j) {
    aj <- a[j] * exp(if (j == dif_item) kappa * g else 0)
    eta <- aj * (theta - b[j]) - (if (j == dif_item) dif_size * g else 0)
    rbinom(n, 1, plogis(eta))
  })
  list(resp = resp, grp = factor(ifelse(g == 1, "B", "A")), g = g)
}

test_that("uniform DIF and impact are separated; matches glmer exactly", {
  skip_if_not_installed("lme4")
  d <- sim_dif_impact(n = 800, ni = 8, dif_size = 0.8)
  res <- dif_irt(d$resp, dif = d$grp, model = "Rasch")

  expect_true(4 %in% res$flagged_items)
  # DIF effect on the logit metric near truth (-0.8)
  expect_equal(res$dif_results$delta_groupB[res$dif_results$item == 4],
               -0.8, tolerance = 0.3)
  # Impact recovered separately from DIF
  expect_equal(res$impact$gamma[1], -0.5, tolerance = 0.25)

  # Cross-walk: identical to the long-format GLMM under the same Laplace
  n <- nrow(d$resp); ni <- ncol(d$resp)
  long <- data.frame(y = as.vector(d$resp),
                     item = factor(rep(1:ni, each = n)),
                     id = factor(rep(1:n, times = ni)),
                     g = rep(d$g, ni))
  long$dif4 <- as.integer(long$item == 4) * long$g
  # bobyqa avoids the "Downdated VtV is not positive definite" breakdown the
  # default nloptwrap/Nelder-Mead PWRSS path can hit on some BLAS/Matrix
  # builds; if the reference fit still fails, skip (never on a gllammr fault).
  ctrl <- lme4::glmerControl(optimizer = "bobyqa")
  # Route through the shared ref_fit helper: it converts reference-fit
  # errors to skips AND gates lme4 off entirely on the Windows CI runner,
  # where lme4 2.0-1/Matrix 1.7-5 segfaults instead of erroring.
  ref_glmer <- function(form)
    ref_fit(lme4::glmer(form, data = long, family = binomial, nAGQ = 1,
                        control = ctrl))
  m1 <- ref_glmer(y ~ 0 + item + g + dif4 + (1 | id))
  m0 <- ref_glmer(y ~ 0 + item + g + (1 | id))
  lr_glmer <- 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0)))
  expect_equal(res$dif_results$chisq[res$dif_results$item == 4],
               lr_glmer, tolerance = 0.05)
})

test_that("pure impact yields no DIF flags", {
  # Group ability difference but NO item bias: the latent regression
  # must absorb it (the scenario where observed-score methods stumble)
  d <- sim_dif_impact(n = 800, ni = 8, seed = 77, dif_size = 0,
                      impact = -0.7, dif_item = 0)
  res <- dif_irt(d$resp, dif = d$grp, model = "Rasch")
  expect_lte(length(res$flagged_items), 1)
  expect_lt(res$impact$gamma[1], -0.4)
})

test_that("multiple DIF factors give multi-df tests", {
  set.seed(91)
  n <- 900; ni <- 8
  pd <- data.frame(gender = factor(sample(c("M", "F"), n, TRUE)),
                   lang = factor(sample(c("X", "Y"), n, TRUE)))
  theta <- rnorm(n) - 0.3 * (pd$gender == "M")
  b <- seq(-1.2, 1.2, length.out = ni)
  resp <- sapply(seq_len(ni), function(j) {
    eta <- theta - b[j]
    if (j == 3) eta <- eta - 0.9 * (pd$lang == "Y")
    rbinom(n, 1, plogis(eta))
  })
  res <- dif_irt(resp, dif = ~ gender + lang, person_data = pd,
                 model = "Rasch")
  expect_true(all(res$dif_results$df == 2))
  expect_true(3 %in% res$flagged_items)
  expect_equal(nrow(res$impact), 2)
})

test_that("wald with anchors agrees with LR on the same hypothesis", {
  d <- sim_dif_impact(n = 800, ni = 8, seed = 33, dif_size = 1.0)
  # Single studied item: the LR and Wald tests address the same model,
  # so the statistics agree closely
  res_lr <- dif_irt(d$resp, dif = d$grp, items = 4,
                    anchors = setdiff(1:8, 4), model = "Rasch")
  res_w <- dif_irt(d$resp, dif = d$grp, items = 4,
                   anchors = setdiff(1:8, 4), model = "Rasch",
                   method = "wald")
  expect_true(4 %in% res_lr$flagged_items)
  expect_true(4 %in% res_w$flagged_items)
  expect_equal(res_w$dif_results$chisq, res_lr$dif_results$chisq,
               tolerance = 0.15)

  # Joint Wald controls for the other studied items' DIF: with true DIF
  # only on item 4, item 5 is not flagged even though the one-at-a-time
  # LR test (item 4 constrained) leaks contamination into it
  res_wj <- dif_irt(d$resp, dif = d$grp, items = 3:5,
                    anchors = c(1, 2, 6, 7, 8), model = "Rasch",
                    method = "wald")
  expect_true(4 %in% res_wj$flagged_items)
  expect_false(5 %in% res_wj$flagged_items)

  expect_error(dif_irt(d$resp, dif = d$grp, method = "wald"),
               "anchors")
})

test_that("nonuniform DIF needs 2PL and is detected", {
  expect_error(
    dif_irt(matrix(rbinom(400, 1, .5), 50, 8),
            dif = factor(rep(0:1, 25)), type = "both"),
    "2PL")
  d <- sim_dif_impact(n = 1200, ni = 8, seed = 99, dif_size = 0,
                      two_pl = TRUE, kappa = 0.7)
  res <- dif_irt(d$resp, dif = d$grp, model = "2PL", type = "both")
  expect_true(4 %in% res$flagged_items)
  expect_true(all(res$dif_results$df == 2))
})

test_that("purified IRT-LR runs and stabilizes", {
  d <- sim_dif_impact(n = 700, ni = 8, seed = 13, dif_size = 1.0)
  res <- dif_irt(d$resp, dif = d$grp, model = "Rasch", purify = TRUE)
  expect_true(res$purification$converged)
  expect_true(4 %in% res$flagged_items)
  out <- capture.output(summary(res))
  expect_true(any(grepl("Impact", out)))
})
