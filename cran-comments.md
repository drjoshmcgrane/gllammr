# CRAN comments for gllammr

## Submission type

New submission. This is the first submission of gllammr to CRAN.

## Test environments

* Local: macOS 15 (Apple Silicon, arm64), R 4.5.1 — tests only. The
  local machine lacks the gfortran toolchain that `R CMD build`/
  `R CMD check` require on this R version, so local verification is
  `testthat::test_local()` against the already-compiled package (all
  unit tests plus the cross-package validation suite; see below).
* GitHub Actions (`R-CMD-check.yaml`), full `R CMD check --as-cran`
  (with `--run-donttest`) on five platform / R-version combinations:
  macOS (release), Windows (release), Ubuntu (devel), Ubuntu (release),
  and Ubuntu (oldrel-1).
* win-builder (R-devel): the tarball was uploaded to win-builder before
  submission; results were still queued at submission time. The GitHub
  Actions matrix above includes a full `--as-cran` check on Windows
  (release) which passes with 0 errors / 0 warnings.

## R CMD check results

0 errors | 0 warnings | 1 note (on CI, `--as-cran`).

Expected NOTEs on CRAN:

* "New submission" — this is the first CRAN release.
* "installed size is ~14Mb (libs ~13Mb)" — the package compiles ~20 TMB
  (Template Model Builder) model templates into a single shared object,
  which is typical for TMB-based packages (cf. glmmTMB).

## Comments

* **First submission.** No reverse dependencies.
* **Compiled code.** All C++ (TMB) templates compile once at install
  time into the package shared object; there is no run-time compilation
  and nothing is cached in the user's home directory. A single
  translation unit (`src/gllammr.cpp`, headers under `src/include/`)
  dispatches between model templates, so compile time and memory are
  dominated by the one `TMB.hpp` include (roughly 1-4 minutes, comparable
  to other TMB packages). OpenMP is used via the standard
  `$(SHLIB_OPENMP_CXXFLAGS)` convention and the package never forces a
  thread count.
* **Cross-package comparison tests are skipped on CRAN by design.** A
  subset of the test suite compares gllammr's numerical output against
  other packages (lme4, mirt, ltm, poLCA, ordinal, npmlreg, lavaan, CDM,
  difR, VGAM, nnet, survival). Their pass/fail depends on those packages'
  numerics, so an upstream release changing an estimate must not archive
  gllammr; and one of them (lme4 2.0-1 with Matrix 1.7-5) can segfault
  inside `glmer()` on Windows. These tests, and the cross-package
  validation suite under `validation/`, carry `skip_on_cran()` and run in
  full on GitHub Actions (where `NOT_CRAN=true`). Every model family
  still has fast unit tests and at least one end-to-end smoke fit that
  runs on CRAN, so the compiled code is exercised on every platform.
* **Test runtime.** With the cross-package and slow integration fits
  skipped, the CRAN-run subset completes well inside the check-time
  budget; the full suite runs on CI.

## Windows compiler note

GCC >= 12 on Windows emits false-positive `-Warray-bounds` warnings from
Eigen's SSE intrinsics included via TMB
(https://gitlab.com/libeigen/eigen/-/issues/2506). The CI workflow
suppresses this on the Windows runner only, via `R_MAKEVARS_USER`
pointed at a generated Makevars that appends `-Wno-array-bounds`; the
flag is not shipped in the package's own `src/Makevars.win`, so it has
no effect on a CRAN build unless CRAN's own Windows toolchain hits the
same Eigen/GCC combination.

## Downstream dependencies

None (new package).
</content>
