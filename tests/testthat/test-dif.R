# Logistic-regression DIF: single/multiple factors, interactions,
# purification, polytomous items, difR agreement

sim_dif <- function(n = 1500, ni = 12, seed = 123,
                    dif_items = integer(0), dif_size = 0.8,
                    nonuniform_items = integer(0), nu_size = 0.8) {
  set.seed(seed)
  theta <- rnorm(n)
  g <- rep(c(0, 1), length.out = n)
  b <- seq(-1.5, 1.5, length.out = ni)
  resp <- matrix(0L, n, ni)
  for (j in seq_len(ni)) {
    eta <- theta - b[j]
    if (j %in% dif_items) eta <- eta - dif_size * g
    if (j %in% nonuniform_items) eta <- eta + nu_size * g * theta
    resp[, j] <- rbinom(n, 1, plogis(eta))
  }
  list(resp = resp, group = factor(ifelse(g == 1, "B", "A")),
       theta = theta)
}

test_that("clean data yields few flags; DIF items are detected", {
  d <- sim_dif(dif_items = c(3, 7))
  res <- dif_test(d$resp, dif = d$group, model = "Rasch")

  expect_s3_class(res, "dif_analysis")
  expect_equal(nrow(res$dif_results), 12)
  expect_true(all(c(3, 7) %in% res$flagged_items))
  # False positives controlled (allow at most 1 of the 10 clean items)
  expect_lte(length(setdiff(res$flagged_items, c(3, 7))), 1)
  # DIF items excluded from the final anchor
  expect_false(any(c(3, 7) %in% res$anchor_items))
})

test_that("uniform vs nonuniform types separate correctly", {
  d <- sim_dif(dif_items = 2, nonuniform_items = 9, seed = 7)
  res_u <- dif_test(d$resp, dif = d$group, type = "uniform",
                    model = "Rasch")
  res_n <- dif_test(d$resp, dif = d$group, type = "nonuniform",
                    model = "2PL")
  expect_true(2 %in% res_u$flagged_items)
  expect_true(9 %in% res_n$flagged_items)
  expect_false(2 %in% res_n$flagged_items)
})

test_that("multiple DIF factors are tested jointly", {
  skip_on_cran()  # large-n multi-factor dif_test fit; smoke kept above
  set.seed(31)
  n <- 2000; ni <- 12
  theta <- rnorm(n)
  pd <- data.frame(gender = factor(sample(c("M", "F"), n, TRUE)),
                   lang = factor(sample(c("X", "Y"), n, TRUE)))
  b <- seq(-1.5, 1.5, length.out = ni)
  resp <- sapply(seq_len(ni), function(j) {
    eta <- theta - b[j]
    if (j == 4) eta <- eta - 0.8 * (pd$gender == "M")   # gender DIF
    if (j == 8) eta <- eta - 0.8 * (pd$lang == "Y")     # language DIF
    rbinom(n, 1, plogis(eta))
  })

  res <- dif_test(resp, dif = ~ gender + lang, person_data = pd,
                  model = "Rasch")
  expect_true(all(c(4, 8) %in% res$flagged_items))
  expect_equal(length(res$dif_terms), 2)
  # df reflects the multi-factor test: both + 2 columns -> 4 df
  expect_true(all(res$dif_results$df == 4))
})

test_that("interaction DIF is detected only with the interaction term", {
  skip_on_cran()  # n=3000 interaction dif_test fit; smoke kept above
  set.seed(41)
  n <- 3000; ni <- 12
  theta <- rnorm(n)
  pd <- data.frame(gender = factor(sample(c("M", "F"), n, TRUE)),
                   lang = factor(sample(c("X", "Y"), n, TRUE)))
  b <- seq(-1.5, 1.5, length.out = ni)
  resp <- sapply(seq_len(ni), function(j) {
    eta <- theta - b[j]
    # Pure intersectional DIF on item 5: only male-Y speakers affected
    if (j == 5) eta <- eta - 1.0 * (pd$gender == "M") * (pd$lang == "Y")
    rbinom(n, 1, plogis(eta))
  })

  res_int <- dif_test(resp, dif = ~ gender * lang, person_data = pd,
                      model = "Rasch")
  expect_true(5 %in% res_int$flagged_items)
  expect_equal(length(res_int$dif_terms), 3)   # two mains + interaction
})

test_that("purification recovers a contaminated matching criterion", {
  # 3 of 14 items share same-direction DIF: matching is contaminated but
  # recoverable
  skip_on_cran()  # iterative purification refits; smoke kept above
  set.seed(11)
  n <- 2000; ni <- 14
  theta <- rnorm(n); g <- rep(c(0, 1), length.out = n)
  b <- seq(-1.5, 1.5, length.out = ni)
  resp <- sapply(seq_len(ni), function(j) {
    eta <- theta - b[j]
    if (j <= 3) eta <- eta - 0.8 * g
    rbinom(n, 1, plogis(eta))
  })
  grp <- factor(ifelse(g == 1, "B", "A"))

  res_p <- dif_test(resp, dif = grp, model = "Rasch", purify = TRUE)
  res_np <- dif_test(resp, dif = grp, model = "Rasch", purify = FALSE)

  expect_true(all(1:3 %in% res_p$flagged_items))
  expect_true(res_p$purification$converged)
  # Purified analysis should not flag more clean items than unpurified
  fp_p <- length(setdiff(res_p$flagged_items, 1:3))
  fp_np <- length(setdiff(res_np$flagged_items, 1:3))
  expect_lte(fp_p, max(fp_np, 1))
})

test_that("purification breakdown degrades gracefully", {
  # A third of the test with strong same-direction DIF: the DIF/impact
  # decomposition is unidentified and purification spirals; the analysis
  # must warn and return the last valid round, not crash
  skip_on_cran()  # n=2500 purification-breakdown refit loop; CI-only
  d <- sim_dif(n = 2500, dif_items = c(1, 2, 3, 4), dif_size = 1.0,
               seed = 11)
  expect_warning(
    res <- dif_test(d$resp, dif = d$group, model = "Rasch",
                    purify = TRUE),
    "anchor")
  expect_s3_class(res, "dif_analysis")
  expect_false(res$purification$converged)
  expect_true(all(is.finite(res$dif_results$p_value)))
})

test_that("score matching reproduces difR::difLogistic flags", {
  skip_if_not_installed("difR")
  skip_on_cran()  # cross-package agreement is CI-only
  d <- sim_dif(n = 1200, dif_items = c(3, 7), seed = 19)
  res <- dif_test(d$resp, dif = d$group, match = "score",
                  type = "both", purify = TRUE)
  ref <- difR::difLogistic(as.data.frame(d$resp), group = d$group,
                           focal.name = "B", type = "both",
                           purify = TRUE)
  ref_flagged <- ref$DIFitems
  if (identical(ref_flagged, "No DIF item detected")) {
    ref_flagged <- integer(0)
  }
  # Flag agreement and high rank correlation of the statistics
  expect_setequal(res$flagged_items, ref_flagged)
  expect_gt(cor(res$dif_results$chisq, ref$Logistik,
                method = "spearman"), 0.95)
})

test_that("anchors are respected and never tested", {
  d <- sim_dif(dif_items = 3, seed = 23)
  res <- dif_test(d$resp, dif = d$group, anchors = c(1, 2),
                  model = "Rasch")
  expect_false(any(c(1, 2) %in% res$dif_results$item))
  expect_true(all(c(1, 2) %in% res$anchor_items))
  expect_true(3 %in% res$flagged_items)
})

test_that("polytomous DIF via cumulative-logit regression", {
  skip_if_not_installed("MASS")
  set.seed(29)
  n <- 1500; ni <- 8
  theta <- rnorm(n)
  g <- rep(c(0, 1), length.out = n)
  resp <- sapply(seq_len(ni), function(j) {
    b0 <- (j - ni / 2) / 2
    eta <- theta - b0 - if (j == 4) 0.9 * g else 0
    cuts <- c(-0.8, 0.8)
    1L + rowSums(outer(eta + rlogis(n), cuts, ">"))
  })
  res <- dif_test(resp, dif = factor(g), model = "GRM")
  expect_true(4 %in% res$flagged_items)
  expect_lte(length(setdiff(res$flagged_items, 4)), 1)
})

test_that("p adjustment, print/summary, and plot work", {
  d <- sim_dif(dif_items = 3, seed = 37)
  res <- dif_test(d$resp, dif = d$group, model = "Rasch",
                  p_adjust = "BH")
  expect_true(all(res$dif_results$p_adj >= res$dif_results$p_value -
                    1e-12))
  out <- capture.output(summary(res))
  expect_true(any(grepl("delta-R2", out)))
  pdf(NULL)
  expect_silent(dif_plot(res, item = 3))
  dev.off()
})

test_that("deprecated wrapper still works", {
  d <- sim_dif(dif_items = 3, seed = 43, n = 800)
  expect_warning(
    res <- dif_test_with_data(d$resp, group = d$group, model = "Rasch"),
    "deprecated")
  expect_s3_class(res, "dif_analysis")
  expect_true(3 %in% res$flagged_items)
})

test_that("input validation", {
  d <- sim_dif(seed = 47, n = 200)
  expect_error(dif_test(d$resp, dif = d$group[1:10]), "must match")
  expect_error(dif_test(d$resp, dif = ~ g), "person_data is required")
  expect_error(dif_test(d$resp, dif = d$group, items = 99), "between")
  expect_error(dif_test(d$resp, dif = d$group, anchors = 99), "valid item")
})
