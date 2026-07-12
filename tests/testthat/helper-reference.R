# Evaluate a reference-package model fit (lme4, ordinal, ...), converting any
# fitting *error* into a testthat skip. Some numerical breakdowns live entirely
# inside the reference package - e.g. lme4's "Downdated VtV is not positive
# definite" PWRSS failure, which surfaces on certain BLAS/Matrix builds - and
# are platform artefacts, not gllammr defects. This helper gates ONLY on the
# reference fit erroring; gllammr's own fit is never wrapped, so a genuine
# gllammr regression still fails the test loudly.
ref_fit <- function(expr) {
  tryCatch(
    force(expr),
    error = function(e)
      testthat::skip(paste("reference fit failed on this platform:",
                           conditionMessage(e))))
}
