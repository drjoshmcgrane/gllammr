# CRAN comments for GLLAMMR 1.2.0

## Submission type

New submission.

## Test environments

* local: macOS (Apple Silicon), R 4.5.1
* R CMD check --as-cran run locally

## R CMD check results

(updated after each check run)

Expected NOTEs:

* "New submission" — first CRAN release.
* Possibly "installed size" — the package compiles ~20 TMB model templates
  into a single shared object; the resulting library is large, which is
  typical for TMB-based packages (cf. glmmTMB).

## Compilation

All C++ (TMB) templates are compiled once at install time into the package
shared object; there is no runtime compilation. A single translation unit
(src/GLLAMMR.cpp) dispatches between model templates, so compile time and
memory are dominated by the one TMB.hpp include (~1–4 minutes, comparable
to other TMB packages).

## Tests and validation

* 905 unit tests run unconditionally (no compilation-dependent skips);
  slow cross-package validation tests are guarded by skip_on_cran() and
  Suggests-availability checks.
* Estimates are cross-validated against lme4, glmmTMB, ordinal, mirt,
  poLCA, npmlreg, and lavaan (49 automated checks; see
  `gllammr_validate()`).

## Downstream dependencies

None (new package).
