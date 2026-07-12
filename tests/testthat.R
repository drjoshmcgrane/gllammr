library(testthat)
library(gllammr)

# Per-file progress markers on stderr (unbuffered), so that a hard native crash
# - e.g. a platform-specific segfault in compiled TMB/Matrix/Eigen code, which
# no R-level tryCatch can catch - leaves the name of the file that was running
# in the check log. A fully buffered stdout otherwise discards that context and
# hides the crash location (as on the Windows CI runner). Built defensively:
# any failure to construct the augmented reporter falls back to a plain run.
reporter <- tryCatch({
  marker <- R6::R6Class(
    "gllammrMarkerReporter", inherit = testthat::Reporter,
    public = list(
      start_file = function(filename) {
        cat(sprintf("[gllammr] running %s\n", filename), file = stderr())
        flush(stderr())
      }
    )
  )$new()
  testthat::MultiReporter$new(
    reporters = list(testthat::check_reporter(), marker))
}, error = function(e) NULL)

if (is.null(reporter)) {
  test_check("gllammr")
} else {
  test_check("gllammr", reporter = reporter)
}
