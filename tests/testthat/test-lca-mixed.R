# E10: LCA with polytomous and continuous indicators.

test_that("polytomous LCA matches poLCA", {
  skip_if_not_installed("poLCA")
  set.seed(101)
  n <- 800
  cls <- sample(1:2, n, TRUE, prob = c(0.6, 0.4))
  P1 <- list(c(.7, .2, .1), c(.1, .3, .6))
  Y <- sapply(1:4, function(j) {
    sapply(cls, function(k) sample(1:3, 1, prob = P1[[k]]))
  })
  colnames(Y) <- paste0("V", 1:4)

  fit <- fit_lca(Y, nclass = 2, control = list(n_starts = 5))

  dl <- as.data.frame(Y)
  f <- as.formula(paste0("cbind(", paste(names(dl), collapse = ","), ") ~ 1"))
  set.seed(1)
  ref <- poLCA::poLCA(f, data = dl, nclass = 2, nrep = 5, verbose = FALSE)

  expect_equal(fit$logLik, ref$llik, tolerance = 0.1)
  expect_equal(sort(unname(fit$class_probs)), sort(ref$P), tolerance = 0.01)
  expect_named(fit$cat_probs, paste0("V", 1:4))
  expect_equal(colSums(fit$cat_probs$V1), c(Class1 = 1, Class2 = 1),
               tolerance = 1e-8)
})

test_that("mixed binary + continuous LCA recovers parameters", {
  set.seed(102)
  n <- 800
  cls <- sample(1:2, n, TRUE, prob = c(0.6, 0.4))
  Y <- cbind(
    b1 = rbinom(n, 1, c(0.8, 0.2)[cls]),
    b2 = rbinom(n, 1, c(0.7, 0.25)[cls]),
    cont = rnorm(n, c(-1, 1.5)[cls], 0.8)
  )
  fit <- fit_lca(Y, nclass = 2, control = list(n_starts = 5))

  expect_true(fit$convergence$converged)
  expect_equal(sort(unname(fit$class_probs)), c(0.4, 0.6), tolerance = 0.06)
  expect_equal(sort(as.numeric(fit$gaussian_params$means)), c(-1, 1.5),
               tolerance = 0.15)
  expect_true(all(abs(fit$gaussian_params$sds - 0.8) < 0.15))
  # Binary item probabilities reported for binary rows only
  expect_true(all(is.na(fit$item_probs["cont", ])))
  expect_true(all(!is.na(fit$item_probs["b1", ])))
})

test_that("binary-only LCA is unchanged", {
  set.seed(103)
  n <- 400
  cls <- sample(1:2, n, TRUE)
  Y <- matrix(rbinom(n * 5, 1, c(0.8, 0.2)[cls]), n, 5)
  fit <- fit_lca(Y, nclass = 2, control = list(n_starts = 3))
  expect_true(fit$convergence$converged)
  expect_true(all(!is.na(fit$item_probs)))
  expect_null(fit$cat_probs)
  expect_null(fit$gaussian_params)
})


test_that("LCA EM (default) and TMB methods agree; EM matches poLCA", {
  skip_if_not_installed("poLCA")
  data("carcinoma", package = "poLCA", envir = environment())
  resp <- as.matrix(carcinoma) - 1L

  set.seed(1)
  fit_em <- fit_lca(resp, nclass = 2, control = list(n_starts = 5))
  set.seed(1)
  fit_tmb <- fit_lca(resp, nclass = 2, method = "tmb",
                     control = list(n_starts = 5))

  expect_equal(fit_em$method, "EM")
  expect_equal(fit_em$logLik, fit_tmb$logLik, tolerance = 1e-3)
  expect_equal(fit_em$logLik, -317.2568, tolerance = 1e-2)  # poLCA reference
  expect_equal(sort(unname(fit_em$class_probs)),
               sort(unname(fit_tmb$class_probs)), tolerance = 1e-2)
})
