# GLLAMMR 1.2.0

## Post-1.2.0 development

* `fit_eirt()` identification fixes: the explanatory GRM now expresses its
  ordered thresholds as sum-to-zero deviations around the item location
  (previously a free first threshold per item silently absorbed the
  difficulty regression, and collided with the item residuals under
  `item_residuals = TRUE`); `sigma_theta` is fixed at 1 whenever a
  discrimination level is estimated (2PL/GRM/GPCM - the `fit_irt`
  convention; the two traded off on a flat ridge); unused `step_param`
  cells are mapped off (previously left free, giving singular Hessians
  and NaN standard errors in polytomous EIRT models).
* New validation cases on the verbal aggression data: `eirt_verbagg`
  (De Boeck & Wilson LLTM+error vs the lme4 crossed-effects GLMM) and
  `eirt_verbagg_pcm` (location-explanatory PCM with random item effects
  vs Kim & Wilson 2019 published estimates).
* `fit_irt()` gains a Bock-Aitkin MML-EM estimation path
  (`method = "em"`), now the **default for single-level models**:
  20-50x faster than the Laplace path, matches mirt to correlation 1.0,
  and handles short tests where joint-Laplace 2PL diverges. The Laplace
  path remains the default whenever multi-level structure (`random`) or
  standard errors (`se = TRUE`) require it. EM abilities are EAP scores.
* `fit_irt(se = FALSE)` is the default (as in mirt); request SEs explicitly.
* Single-pass AD tapes in the gaussian/binomial/poisson templates.

## Build and infrastructure

* All TMB C++ templates are now compiled at install time into a single
  shared library (CRAN-compatible). No run-time compilation, no compiled
  artifacts cached in the user's home directory, and no C++ toolchain
  required after installation.

## Full GLLAMM parity features

* Random slopes `(x | g)` and uncorrelated slopes `(x || g)` with full
  variance-covariance estimation across Gaussian, binomial, Poisson, and
  gamma families.
* Crossed `(1 | g1) + (1 | g2)` and nested `(1 | g1/g2)` random effects in
  the general estimator.
* Multilevel explanatory IRT: person-level random effects
  (`random = ~ (1 | class)` etc.) and person-level random slopes in IRT and
  EIRT models, including nested and crossed structures.
* Gamma family (log link) for positive continuous responses.
* Parametric survival/frailty models (`fit_survival()`; exponential and
  Weibull with shared log-normal frailty).
* Joint mixed-response models sharing a random effect across Gaussian,
  binomial, and Poisson outcomes (`fit_mixed()`).
* Structural equation models with measurement and structural parts
  (`fit_sem()`).
* Rank-ordered (exploded) logit with optional random coefficients
  (`fit_rank()`).
* Level-specific survey weights via
  `weights = list(level1 = , level2 = )` (pseudo-likelihood, matching
  Stata GLLAMM `pweight()` behavior).
* Nonparametric maximum likelihood (NPML) mass-point random effects
  (`fit_npml()`), validated against npmlreg.
* Adaptive Gauss-Hermite quadrature via `integration = aghq(k)` as an
  alternative to the default Laplace approximation.
* Cluster-robust (sandwich) standard errors via
  `vcov(fit, type = "sandwich")`.
* Latent class analysis with mixed indicator types (binary, categorical,
  and continuous manifest variables).

## Unified interface

* All model classes are reachable through `gllamm()` family dispatch:
  `irt()`, `lca()`, `ordinal()`, `multinomial()`, and `binomial()` family
  objects route to the appropriate fitter, e.g.
  `gllamm(resp_matrix, family = irt("2PL"))` or
  `gllamm(y ~ x + (1 | g), data, family = ordinal("probit"))`.
* `predict(fit, type = "marginal")` population-averaged predictions and
  `simulate()` (returning a data frame) across model classes.

## Identification and correctness fixes

* Ordinal models: intercept absorbed into thresholds (standard
  identification), fixing previously shifted threshold estimates.
* IRT: latent variance fixed at 1 for 2PL/3PL/GRM/GPCM/NRM so
  discriminations are identified; Rasch frees the latent variance.
* LCA: corrected parameterization of class-membership and item-response
  logits; polytomous coding fixed.
* GRM: threshold orientation aligned with mirt conventions.

## Performance

* 40-70x speedups for standard GLMMs from dedicated TMB templates, sparse
  random-effects design matrices, and better starting values.

## Validation

* New cross-package validation suite, `gllammr_validate()`: 49 automated
  checks of estimates, variance components, and log-likelihoods against
  lme4, glmmTMB, ordinal, mirt, poLCA, npmlreg, and lavaan.

# GLLAMMR 1.1.0 and earlier

* Initial development releases: Gaussian/binomial/Poisson GLMMs, ordinal
  and multinomial models, dichotomous and polytomous IRT (Rasch, 2PL, 3PL,
  GRM, PCM, GPCM, NRM), DIF analysis, explanatory IRT, latent class
  analysis, marginal predictions, and frequency/probability weights.
