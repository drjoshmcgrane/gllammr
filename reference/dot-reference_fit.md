# Evaluate a reference-package model fit, converting a fitting error into a platform skip. Numerical breakdowns inside a reference package (e.g. lme4's "Downdated VtV is not positive definite" on some BLAS/Matrix builds) then mark the case skipped, never failed. \`expr\` is evaluated lazily inside the handler.

Evaluate a reference-package model fit, converting a fitting error into
a platform skip. Numerical breakdowns inside a reference package (e.g.
lme4's "Downdated VtV is not positive definite" on some BLAS/Matrix
builds) then mark the case skipped, never failed. \`expr\` is evaluated
lazily inside the handler.

## Usage

``` r
.reference_fit(expr)
```
