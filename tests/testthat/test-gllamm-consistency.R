# GLLAMM-framework internal consistency: the unification claims hold
# numerically within the package itself.

test_that("CDM with one attribute is exactly a 2-class LCA", {
  set.seed(11)
  n <- 600; J <- 6
  cls <- rbinom(n, 1, 0.45)
  pmat <- rbind(c(0.2, 0.8), c(0.3, 0.7), c(0.15, 0.85),
                c(0.25, 0.9), c(0.1, 0.6), c(0.35, 0.75))
  Y <- sapply(1:J, function(j) rbinom(n, 1, pmat[j, cls + 1]))
  f_cdm <- fit_cdm(Y, matrix(1, J, 1), model = "gdina", monotone = FALSE,
                   control = list(n_starts = 6, tol = 1e-9))
  f_lca <- fit_lca(Y, nclass = 2, control = list(n_starts = 6, tol = 1e-9))
  expect_equal(f_cdm$logLik, f_lca$logLik, tolerance = 1e-7)
})

test_that("EIRT with item residuals equals the crossed-RE GLLAMM", {
  skip_if_not_installed("lme4")
  data("VerbAgg", package = "lme4", envir = environment())
  VerbAgg$y <- as.integer(VerbAgg$r2 == "Y")
  resp <- with(VerbAgg, tapply(y, list(id, item), identity))
  resp <- matrix(as.integer(resp), nrow = nrow(resp))
  f_eirt <- fit_eirt(resp, data.frame(int = rep(1, ncol(resp))),
                     difficulty_formula = ~ 1, model = "Rasch",
                     item_residuals = TRUE)
  f_g <- gllamm(y ~ 1 + (1 | id) + (1 | item), data = VerbAgg,
                family = binomial())
  expect_equal(f_eirt$logLik, f_g$logLik, tolerance = 1e-3)
  sds <- sqrt(unlist(lapply(f_g$coefficients$random_var,
                            function(m) m[1, 1])))
  expect_equal(unname(f_eirt$ability_sd), unname(sds[1]), tolerance = 1e-3)
  expect_equal(unname(f_eirt$residual_sd$difficulty), unname(sds[2]),
               tolerance = 1e-3)
})
