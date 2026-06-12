# Generic model comparison across gllammr model classes

test_that("compare_models builds a coherent table across model classes", {
  set.seed(5)
  np <- 400; ni <- 10
  theta <- rnorm(np)
  b <- seq(-1.5, 1.5, length.out = ni)
  a <- runif(ni, 0.6, 1.8)
  resp <- sapply(1:ni, function(j) {
    rbinom(np, 1, plogis(a[j] * (theta - b[j])))
  })

  rasch <- fit_irt(resp, model = "Rasch")
  twopl <- fit_irt(resp, model = "2PL")
  lca2 <- fit_lca(resp, nclass = 2, control = list(n_starts = 3))

  cmp <- compare_models(rasch = rasch, twopl = twopl, lca2 = lca2)
  expect_s3_class(cmp, "gllammr_model_comparison")
  expect_equal(cmp$model, c("rasch", "twopl", "lca2"))
  expect_equal(nrow(cmp), 3)
  # Akaike weights are a probability distribution
  expect_equal(sum(cmp$akaike_weight), 1, tolerance = 1e-10)
  # Parameter counts recovered from AIC identity
  expect_equal(cmp$n_params[cmp$model == "twopl"],
               cmp$n_params[cmp$model == "rasch"] + ni - 1)
  # 2PL generated the data: it should win AIC here
  expect_equal(cmp$model[which.min(cmp$AIC)], "twopl")
  # Deltas are zero at the minimum
  expect_equal(min(cmp$dAIC), 0)
  expect_equal(min(cmp$dBIC), 0)

  out <- capture.output(print(cmp))
  expect_true(any(grepl("Best by AIC", out)))
})

test_that("compare_models sorts and labels unnamed arguments", {
  set.seed(7)
  resp <- matrix(rbinom(2000, 1, 0.5), 200, 10)
  f1 <- fit_irt(resp, model = "Rasch")
  f2 <- fit_lca(resp, nclass = 2, control = list(n_starts = 2))
  cmp <- compare_models(f1, f2, sort_by = "BIC")
  expect_equal(cmp$BIC, sort(cmp$BIC))
  expect_true(all(nzchar(cmp$model)))
})

test_that("compare_models warns on differing observation counts", {
  set.seed(9)
  resp <- matrix(rbinom(3000, 1, 0.5), 300, 10)
  f_full <- fit_irt(resp, model = "Rasch")
  f_half <- fit_irt(resp[1:150, ], model = "Rasch")
  expect_warning(compare_models(full = f_full, half = f_half),
                 "different numbers of observations")
  expect_error(compare_models(f_full), "at least two")
})

test_that("latent_structure_comparison delegates to the generic table", {
  set.seed(41)
  theta <- c(-1, 0, 1)
  delta <- seq(1, -1, length.out = 6)
  cls <- sample(1:3, 600, TRUE)
  Y <- sapply(1:6, function(j) rbinom(600, 1, plogis(theta[cls] - delta[j])))
  cmp <- latent_structure_comparison(Y, nclass = 3, n_starts = 2)
  expect_s3_class(cmp, "gllammr_model_comparison")
  expect_true(all(c("dAIC", "dBIC") %in% names(cmp)))
})
