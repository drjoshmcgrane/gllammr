# CRAN comments for gllammr

## Submission type

New submission.

## Test environments

* Local: macOS 15 (Apple Silicon), R 4.5.1. The local machine lacks the
  gfortran toolchain `R CMD build`/`R CMD check` require on this R
  version, so local verification is tests-only: `testthat::test_local()`
  against the already-compiled package (all unit tests and the
  cross-package validation suite; see below).
* GitHub Actions (`R-CMD-check.yaml`), full `R CMD check --as-cran` on
  five platform/R-version combinations: macOS (release), Windows
  (release), Ubuntu (devel), Ubuntu (release), Ubuntu (oldrel-1).

## R CMD check results

Not run locally (see above); CI is the source of truth for `R CMD check`
on this submission.

Expected NOTEs on CRAN:

* "New submission" — first CRAN release.
* "installed size is ~14Mb (libs ~13Mb)" — the package compiles ~20 TMB
  model templates into a single shared object, which is typical for
  TMB-based packages (cf. glmmTMB).

## Windows compiler note

GCC >= 12 on Windows emits false-positive `-Warray-bounds` warnings from
Eigen's SSE intrinsics included via TMB
(https://gitlab.com/libeigen/eigen/-/issues/2506). The CI workflow
suppresses this on the Windows runner only, via `R_MAKEVARS_USER`
pointed at a generated Makevars that appends `-Wno-array-bounds`; the
flag is not shipped in the package's own `src/Makevars.win`, so it has
no effect on a CRAN build unless CRAN's own Windows toolchain hits the
same Eigen/GCC combination.

## Compilation

All C++ (TMB) templates compile once at install time into the package
shared object; there is no runtime compilation. A single translation
unit (src/gllammr.cpp, headers under src/include/) dispatches between
model templates, so compile time and memory are dominated by the one
TMB.hpp include (1-4 minutes, comparable to other TMB packages).

## Tests and validation

* Unit tests (`testthat::test_local()`): 3400+ tests, 0 failures, all
  skips justified (each remaining `skip()` guards a genuinely
  unavailable Suggests package or environment, e.g. `skip_if_not_installed()`).
* Cross-package validation suite (`validation/run_validation.R`,
  guarded by `skip_on_cran()` and Suggests-availability checks in the
  test suite): 93/93 checks pass against lme4, glmmTMB, ordinal, mirt,
  poLCA, npmlreg, and lavaan. Reference fits using the identical
  Laplace approximation (lme4 nAGQ = 1, ordinal::clmm) agree to between
  1e-5 and machine precision; packages using different integration
  (mirt, ltm) agree within stated, looser tolerances.

## Downstream dependencies

None (new package).
