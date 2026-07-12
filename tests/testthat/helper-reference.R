# Evaluate a reference-package model fit (lme4, ordinal, ...), converting any
# fitting *error* into a testthat skip. Some numerical breakdowns live entirely
# inside the reference package - e.g. lme4's "Downdated VtV is not positive
# definite" PWRSS failure, which surfaces on certain BLAS/Matrix builds - and
# are platform artefacts, not gllammr defects. This helper gates ONLY on the
# reference fit erroring; gllammr's own fit is never wrapped, so a genuine
# gllammr regression still fails the test loudly.
ref_fit <- function(expr) {
  # Cross-package numerical agreement is a continuous-integration concern, not a
  # CRAN one: an upstream release changing another package's estimates (or the
  # known lme4 2.0-1/Matrix 1.7-5 Windows segfault below) must never archive
  # gllammr. skip_on_cran() no-ops when NOT_CRAN=true, so CI still runs these.
  testthat::skip_on_cran()
  # lme4 2.0-1 with Matrix 1.7-5 segfaults (not errors) inside glmer/vcov on
  # the Windows GitHub runner, killing the R process before tryCatch can act
  # - confirmed with a from-source lme4 build, so it is an upstream bug, not
  # an ABI artefact. Skip reference fits there; they still run on the other
  # four CI platforms and everywhere locally.
  if (identical(Sys.getenv("GITHUB_ACTIONS"), "true") &&
      identical(.Platform$OS.type, "windows")) {
    testthat::skip("lme4 reference fits segfault on the Windows CI runner (upstream lme4/Matrix bug)")
  }
  tryCatch(
    force(expr),
    error = function(e)
      testthat::skip(paste("reference fit failed on this platform:",
                           conditionMessage(e))))
}
