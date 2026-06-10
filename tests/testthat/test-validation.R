# Cross-package validation: every case must pass its stated tolerance.
# Skipped on CRAN (reference packages are Suggests and the fits are slow).

test_that("GLLAMMR estimates agree with reference packages", {
  skip_on_cran()
  skip_if_not_installed("lme4")

  set.seed(2026)
  res <- gllammr_validate(verbose = FALSE)

  failures <- res[!res$pass, c("case", "statistic", "gllammr", "reference")]
  expect_true(all(res$pass),
              info = paste(capture.output(print(failures)), collapse = "\n"))
  expect_gte(nrow(res), 20)   # the suite actually ran, not all skipped
})
