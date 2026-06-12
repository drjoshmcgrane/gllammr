# Latent structure models: invariant item ordering, double monotonicity,
# located latent classes (Torres Irribarra & Diakow framework)

sim_located <- function(n = 1000, J = 8, K = 4, seed = 9) {
  set.seed(seed)
  theta <- seq(-1.5, 1.5, length.out = K)
  delta <- seq(1.4, -1.4, length.out = J)
  cls <- sample(1:K, n, TRUE)
  Y <- sapply(1:J, function(j) {
    rbinom(n, 1, plogis(theta[cls] - delta[j]))
  })
  list(Y = Y, theta = theta, delta = delta, cls = cls)
}

test_that("invariant item ordering constrains items within classes", {
  d <- sim_located()
  J <- ncol(d$Y)
  iio <- fit_lca(d$Y, nclass = 4,
                 item_ordering = lapply(1:(J - 1), function(k) c(k, k + 1)),
                 control = list(n_starts = 5))
  # Nondecreasing across items within every class
  expect_true(all(apply(iio$item_probs, 2,
                        function(col) all(diff(col) >= -1e-8))))
  # Classes themselves need not be ordered under IIO alone
  expect_s3_class(iio, "gllamm_lca")
})

test_that("'increasing' item_ordering shorthand uses column order", {
  d <- sim_located(seed = 21)
  iio <- fit_lca(d$Y, nclass = 3, item_ordering = "increasing",
                 control = list(n_starts = 3))
  expect_true(all(apply(iio$item_probs, 2,
                        function(col) all(diff(col) >= -1e-8))))
})

test_that("double monotonicity holds in both directions", {
  d <- sim_located(seed = 11)
  J <- ncol(d$Y)
  dm <- fit_lca(d$Y, nclass = 4, ordering = "increasing",
                item_ordering = lapply(1:(J - 1), function(k) c(k, k + 1)),
                control = list(n_starts = 5))
  expect_true(all(apply(dm$item_probs, 1,
                        function(row) all(diff(row) >= -1e-8))))
  expect_true(all(apply(dm$item_probs, 2,
                        function(col) all(diff(col) >= -1e-8))))
  # Nesting: DM can never beat the single-monotonicity models
  mon <- fit_lca(d$Y, nclass = 4, ordering = "increasing",
                 control = list(n_starts = 5))
  expect_lte(dm$logLik, mon$logLik + 1e-6)
})

test_that("located latent classes equal the NPML Rasch model", {
  d <- sim_located(n = 1000, seed = 9)
  J <- ncol(d$Y)
  lcr <- fit_lca(d$Y, nclass = 4, structure = "rasch",
                 control = list(n_starts = 5))

  # Lindsay-Clogg-Grego: the latent class Rasch model IS the Rasch model
  # with a k-point nonparametric ability distribution
  long <- data.frame(y = as.vector(d$Y),
                     item = factor(rep(1:J, each = nrow(d$Y))),
                     id = factor(rep(seq_len(nrow(d$Y)), J)))
  np <- fit_npml(y ~ 0 + item + (1 | id), data = long, k = 4,
                 family = binomial())
  expect_equal(lcr$logLik, np$logLik, tolerance = 0.1)

  # Locations sorted ascending, parameters recovered
  expect_true(all(diff(lcr$class_locations) > 0))
  expect_lt(max(abs(sort(lcr$item_difficulties) - sort(d$delta))), 0.3)
  expect_equal(unname(sum(lcr$item_difficulties)), 0, tolerance = 1e-6)
})

test_that("structure and ordering arguments are validated", {
  d <- sim_located(n = 300, seed = 31)
  expect_error(fit_lca(d$Y, nclass = 3, structure = "rasch",
                       ordering = "increasing"),
               "implies")
  Ymix <- cbind(d$Y[, 1:3], rnorm(300))
  expect_error(fit_lca(Ymix, nclass = 3, structure = "rasch"),
               "binary")
  expect_error(fit_lca(Ymix, nclass = 3, item_ordering = "increasing"),
               "binary")
  expect_error(fit_lca(d$Y, nclass = 3,
                       item_ordering = list(c(1, 2), c(2, 1))),
               "cycle")
})

test_that("latent_structure_comparison fits all six models coherently", {
  d <- sim_located(n = 800, J = 6, seed = 41)
  cmp <- latent_structure_comparison(d$Y, nclass = 3, n_starts = 3)

  expect_s3_class(cmp, "lca_structure_comparison")
  expect_equal(cmp$model, c("UN", "MON", "IIO", "DM", "LCR", "RM"))
  # Nesting in fit: UN >= MON >= DM and UN >= IIO >= DM >= LCR
  ll <- setNames(cmp$logLik, cmp$model)
  expect_gte(ll["UN"], ll["MON"] - 1e-6)
  expect_gte(ll["MON"], ll["DM"] - 1e-6)
  expect_gte(ll["IIO"], ll["DM"] - 1e-6)
  expect_gte(ll["DM"], ll["LCR"] - 1e-6)
  # UN..DM share the nominal parameter count; LCR and RM have fewer
  expect_equal(length(unique(cmp$n_params[1:4])), 1)
  expect_lt(cmp$n_params[5], cmp$n_params[1])
  expect_lt(cmp$n_params[6], cmp$n_params[5])
  # On located-class data a quantitative model wins
  expect_true(cmp$model[which.min(cmp$BIC)] %in% c("LCR", "RM"))
  out <- capture.output(print(cmp))
  expect_true(any(grepl("interval scale", out)))
})
