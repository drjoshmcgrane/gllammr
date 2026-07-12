# Cross-package validation: every case must pass its stated tolerance.
# Skipped on CRAN (reference packages are Suggests and the fits are slow).

test_that("gllammr estimates agree with reference packages", {
  skip_on_cran()
  skip_if_not_installed("lme4")

  set.seed(2026)
  # verbose: per-case progress on unbuffered stderr, so a hard native crash
  # in one case (as on the Windows CI runner) names the case in the check log
  res <- gllammr_validate(verbose = TRUE)

  # Rows with pass == NA are reference-package skips (a numerical breakdown
  # inside lme4 etc. on some platforms) - not gllammr failures. Every row
  # that actually completed a comparison must pass its stated tolerance.
  completed <- !is.na(res$pass)
  failures <- res[completed & !res$pass, c("case", "statistic", "gllammr",
                                           "reference")]
  expect_true(all(res$pass[completed]),
              info = paste(capture.output(print(failures)), collapse = "\n"))

  # Errors must only ever originate in a reference package (surfaced as SKIP);
  # a genuine gllammr-side error is recorded as ERROR and must never occur.
  errors <- res[res$statistic == "ERROR", c("case", "note")]
  expect_true(nrow(errors) == 0,
              info = paste(capture.output(print(errors)), collapse = "\n"))

  expect_gte(sum(completed), 20)   # the suite actually ran, not all skipped
})
