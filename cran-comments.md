# CRAN comments for GLLAMMR 1.2.0

## Submission type

New submission.

## Test environments

* local: macOS 15 (Apple Silicon), R 4.5.1, clang 17
* R CMD check --as-cran (final local run: 2026-06-11)

## R CMD check results

0 ERRORs, 0 WARNINGs attributable to the package.

Local-toolchain artifacts on the development machine (absent on machines
with a TeX distribution / current HTML Tidy):

* "PDF version of manual" ERROR/WARNING — pdflatex is not installed
  locally; `checking Rd files ... OK` passes and all Rd content issues
  were resolved.
* "unable to verify current time", HTML Tidy version, leftover
  GLLAMMR-manual.tex — consequences of the same local toolchain.

Expected NOTEs on CRAN:

* "New submission" — first CRAN release.
* "installed size is 14.5Mb (libs 12.8Mb)" — the package compiles ~20 TMB
  model templates into a single shared object, which is typical for
  TMB-based packages (cf. glmmTMB).

## Compilation

All C++ (TMB) templates compile once at install time into the package
shared object; there is no runtime compilation. A single translation unit
(src/GLLAMMR.cpp, headers under src/include/) dispatches between model
templates, so compile time and memory are dominated by the one TMB.hpp
include (1–4 minutes, comparable to other TMB packages).

## Tests and validation

* 905 unit tests run unconditionally; the slow cross-package validation
  suite is guarded by skip_on_cran() and Suggests-availability checks.
* Estimates are cross-validated against lme4, glmmTMB, ordinal, mirt,
  poLCA, npmlreg, and lavaan — 49 automated agreement checks (all
  passing; see `gllammr_validate()`). Reference fits using the identical
  Laplace approximation (lme4 nAGQ = 1, ordinal::clmm) agree to between
  1e-5 and machine precision.

## Downstream dependencies

None (new package).
