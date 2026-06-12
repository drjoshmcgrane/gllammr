# gllammr 1.2.0

## Post-1.2.0 development

* GLMM level-2 weight audit (follow-up to the EIRT weight finding).
  Level-2 (group) weights under Laplace previously scaled each group's
  likelihood-plus-prior contribution - approximate in the variance
  components and unbounded in principle (every weighted group
  contributes -(w-1)*log(sigma_u)). Integer level-2 frequency weights
  are now implemented by exact whole-group replication before fitting
  (weighted fits are identical to group-duplicated fits, verified for
  gaussian/poisson/binomial; replication happens before listwise
  deletion, so incomplete rows are handled exactly as duplicated data
  would be); non-integer level-2 weights under Laplace are rejected
  with guidance to `integration = aghq(k)`, which weights each group's
  log marginal likelihood outside the integral and is exact for
  arbitrary weights (verified identical to duplication). Binomial fits
  with list weights now route to the general engine instead of failing
  in `fit_binomial()`; ordinal models reject level-specific weights
  with a clear message; survival models likewise reject list weights
  clearly (previously an obscure coercion error).
* EIRT deep audit. **Saturated identities verified exactly:** an EIRT
  model with `difficulty_formula = ~ item` (and, where applicable,
  `discrimination_formula = ~ item`) reproduces the descriptive
  `fit_irt()` likelihood for Rasch, 2PL, PCM, GPCM, and GRM; multilevel
  Rasch EIRT matches `lme4::glmer` on the crossed person/school
  cross-walk (logLik, variance components, and standard errors).
  **Weight semantics fixed:** person-level weights on Laplace-based fits
  (`fit_eirt()`, and `fit_irt(method = "laplace")`) previously weighted
  only the response terms, leaving the ability prior unweighted - not a
  valid frequency- or sampling-weight likelihood. Weighting the prior
  instead is degenerate under the Laplace approximation (each weighted
  person contributes -(w-1)*log(sigma_theta), so the objective is
  unbounded and sigma_theta collapses to 0). Integer frequency weights
  are now implemented by exact replication of weighted persons, making
  weighted fits identical to duplicated-data fits (verified to 1e-6,
  single-level and multilevel, dichotomous and polytomous); non-integer
  person weights on Laplace paths are rejected with guidance to
  `method = "em"`, which weights each person's log marginal likelihood
  directly (already exact). **Multilevel method fixes:** `simulate()` on
  multilevel IRT/EIRT fits drew only the person deviation and ignored
  group random effects (school-level variance was ~4x too small);
  marginal predictions likewise integrated over `sigma_theta` only -
  both now use the full latent structure. **New guards:** rank-deficient
  `difficulty_formula`/`discrimination_formula`/`threshold_formula`
  designs error early naming the aliased columns (previously NaN
  standard errors); `threshold_formula` with GRM (or dichotomous
  models) errors instead of being silently ignored.
* CI: `nnet` added to Suggests (used by audit tests);
  `-Wno-array-bounds` moved out of the shipped `Makevars.win` (it
  triggered R CMD check's non-portable-flags warning) onto the Windows
  CI runner only.
* Package-wide audit (every S3 method probed on every fit variant;
  under-validated likelihoods checked externally). **Bug fixes:**
  on crossed/multi-term GLMMs, `simulate()` drew random effects for the
  first term only and marginal predictions integrated over only the
  first term's variance (both silently wrong; conditional newdata
  predictions also applied only the first term's BLUPs) - all three now
  loop over every term; `summary()` recycled mismatched covariance
  matrices into wrong standard errors for survival/rank/NPML fits (SEs
  are now shown only when a matching vcov exists); `vcov()` errored on
  SEM fits instead of returning the parameter covariance;
  `parse_formula()` mangled call-type responses like `Surv(t, d)`.
  **New methods:** `simulate()` for multinomial, survival (uncensored
  event times from the fitted hazard), IRT (parametric bootstrap,
  dichotomous + polytomous, EM and Laplace fits), EIRT, LCA (all
  indicator types), CDM, SEM (implied-distribution draws), NPML
  (mass-point draws), and mixed responses; `predict()` for polytomous
  IRT (category probabilities and expected scores); multinomial fits now
  store their random effects so `ranef()` works. `icc`/`VarCorr`/`ranef`
  refusals now point at the right accessor for each latent-variable
  class. **External validations added:** multinomial matches
  `nnet::multinom` in the fixed-effect limit (1e-4); Weibull frailty
  matches `survival::survreg` in the no-frailty limit (with the
  AFT-scale parameterization now documented: log hazard ratio =
  shape x beta); rank-ordered logit coefficients equal the exploded
  conditional logit (`survival::coxph`, Plackett-Luce equivalence,
  1e-4).
* Ordinal GLMM audit. **Bug fix:** the backward continuation-ratio link
  (`crl_backward`) was not a valid probability model - its category
  probabilities summed to 2*plogis(tau_max - eta) because the top
  category reused the last cumulative probability; it now uses the
  proper hazard-product form and, like `acl` and `crl_forward`, is
  validated against VGAM (acat/sratio; agreement ~1e-3). `predict()` and
  `simulate()` now support all six ordinal links (previously
  logit/probit only) through a shared category-probability helper
  mirroring the C++ likelihoods, and both - along with marginal
  predictions - handle crossed/multi-term random effects (marginal
  predictions previously integrated over only the first term's variance;
  simulate() errored). `icc()` gains a proper latent-response-scale
  branch for cumulative ordinal models (it previously fell through to
  the gaussian formula); threshold starting values are computed from raw
  spacings. Single-term ordinal fits run at ordinal::clmm speed parity
  with identical likelihoods. VGAM added to Suggests.
* New general `compare_models()`: a comparison table (logLik, parameter
  count, AIC/BIC with deltas, Akaike weights) for any set of fitted
  gllammr models, across model classes, with an n_obs comparability
  check; `latent_structure_comparison()` now delegates its table to it.
* Order-restricted LCA optimized: the isotonic poset projection runs in
  C (shared with the CDM engine), and the Ramsay acceleration trigger
  window widened to cover slowly converging EM (step ratio near 1; the
  logLik-revert safeguard makes this safe). MON/IIO/DM fits are 2.7-4.2x
  faster; the located class model now converges where it previously hit
  the iteration cap, and its EM is ~10x faster than fitting the same
  model through the equivalent NPML GLMM route (identical logLik).
* Latent structure analysis (Torres Irribarra & Diakow framework):
  `fit_lca()` gains `item_ordering` (invariant item ordering - item
  monotonicity within classes, via isotonic regression over the
  item-by-class grid; combined with `ordering` this is the double
  monotonicity model) and `structure = "rasch"` (located latent
  classes / latent class Rasch: logit pi_ic = theta_c - delta_i, classes
  on an interval scale, reported sorted by location with
  `class_locations` and `item_difficulties`). The located class model's
  likelihood is verified against `fit_npml()` on the long-format GLMM
  (the Lindsay-Clogg-Grego equivalence, to 0.02 logLik). New
  `latent_structure_comparison()` fits all six models of the framework
  (UN, MON, IIO, DM, LCR, RM) and returns the successive-comparison
  table for deciding between qualitative, ordinal, and quantitative
  latent structure.
* Vignette suite completed: four new evaluated vignettes -
  `explanatory-irt` (De Boeck & Wilson on the verbal aggression data,
  with the GLMM cross-walk and the Kim & Wilson polytomous extension),
  `cognitive-diagnosis` (Q-matrices, DINA/DINO/G-DINA, monotonicity,
  hierarchies), `sem-models` (CFA, fit indices, MIMIC, FIML,
  standardized solution), and `dif-analysis` (the two-stage
  screen-then-confirm DIF workflow). `latent-class` rewritten with
  evaluated code covering ordered and partially ordered classes;
  `advanced-features` rewritten with evaluated examples of mixed
  responses, frailty survival, rankings, NPML, AGHQ, and sandwich
  standard errors through the unified interface; `getting-started`
  gains the complete model-space table. All thirteen vignettes build
  in under two minutes. `fit_lca()` now rejects `nclass = 1` with an
  informative message.
* The unified `gllamm()` interface now reaches every model class: new
  family constructors `eirt()`, `sem()`, `mixed_response()`,
  `survival_family()`, and `ranking()`, plus `integration = npml(k)` for
  nonparametric (mass-point) latent distributions. Also fixed:
  `family = binomial()` with `integration = aghq(k)` or `npml(k)` was
  silently ignoring the integration request (the binomial fast path
  intercepted dispatch); it now routes through the requested engine.
* New `dif_irt()`: confirmatory model-based DIF (IRT-LR; Thissen et al.)
  as the companion to the `dif_test` screening tests. DIF parameters are
  item-by-covariate interactions inside the joint marginal-ML Rasch/2PL
  model, with a latent-regression impact term separating true ability
  differences from item bias; supports multiple DIF variables and
  interactions via the same formula interface, uniform and (2PL)
  nonuniform DIF, per-item LR tests or joint-model Wald tests with
  anchors, and purified IRT-LR. For the Rasch model with uniform DIF the
  tests are likelihood-identical to the De Boeck & Wilson long-format
  GLMM (`y ~ 0 + item + z + item_j:z + (1|person)`), verified against
  glmer (validation case `dif_irt_glmm`). Estimated DIF effects are
  reported on the logit metric with standard errors.
* DIF analysis rewritten (`dif_test`): logistic-regression DIF with a
  latent (EAP) or observed-score matching criterion, a formula interface
  for **multiple DIF variables and their interactions**
  (`dif = ~ gender * language`), **iterative purification** of the
  matching criterion (with graceful degradation and a warning when the
  DIF/impact decomposition is unidentified), uniform/nonuniform/joint
  LR tests with correct per-item degrees of freedom, Nagelkerke
  delta-R2 effect sizes with the Jodoin-Gierl A/B/C classification,
  optional multiplicity adjustment, anchor-item support, and
  cumulative-logit tests for polytomous items. With score matching the
  flags reproduce `difR::difLogistic` exactly (validation case
  `dif_logistic`). The previous implementation divided one global LR
  statistic equally across items and confounded DIF with impact via
  separate per-group calibrations; `dif_test_with_data` is deprecated
  and now wraps the new engine. `dif_plot` draws model-implied response
  curves by any DIF variable.
* SEM overhaul (`fit_sem`): exogenous latent variables now covary freely
  (previously silently orthogonal - a plain two-factor CFA was
  misspecified); MIMIC models (structural regressions on observed
  covariates, likelihood-equivalent to lavaan with `fixed.x = FALSE`);
  standard errors for every parameter (numerical observed information +
  delta method); the standard fit-index battery (chisq, CFI, TLI, RMSEA
  with 90% CI, SRMR - matching lavaan to 4 decimals); full-information
  maximum likelihood for missing data (`missing = "fiml"`, pattern-based,
  with an EM-estimated saturated model for fit indices); a standardized
  (std.all) solution; and a real `summary()` method. The legacy Laplace
  path warns that it treats exogenous factors as orthogonal.
* Package-wide missing-data policy, audited and tested: formula-based
  fitters (GLMM, ordinal, multinomial, survival, NPML, mixed responses)
  now listwise-delete rows with missing values in any model variable -
  with a warning and automatic weight alignment. Previously the response
  and fixed-effect design dropped NA rows while the random-effects
  design and grouping factor kept them, silently misaligning ordinal,
  survival, and NPML fits. Matrix-response latent variable models
  (IRT/EIRT/LCA/CDM) already used all observed responses (MAR) and are
  unchanged; `fit_rank` keeps its partial-ranking semantics for missing
  ranks.
* Multiple (crossed/nested) random-effects terms for ordinal and
  multinomial models: `fit_ordinal(y ~ x + (1 | rater) + (1 | item))`
  works for the logit, probit, adjacent-category, and continuation-ratio
  links (PPO remains single-term), matching `ordinal::clmm` with the
  same crossed structure to 1e-4; `fit_multinomial()` gains the same
  layout with the random effects acting as a common shifter across
  non-reference categories. The same `(1 | g1) + (1 | g2)` formulas now
  also work for `family = binomial()` through `gllamm()`.
* New `cdm()` family constructor and `lca(ordering = ...)` pass-through:
  cognitive diagnosis models and order-restricted latent class models
  are now reachable through the unified `gllamm()` interface.
* New `fit_cdm()`: cognitive diagnosis models for binary responses with
  a Q-matrix - saturated G-DINA (default), DINA, and DINO - with
  monotonicity in the attributes enforced by isotonic regression over
  the reduced-profile lattice (`monotone = TRUE`, the default) and
  optional attribute hierarchies (prerequisite relations prune the
  profile space). Closed-form accelerated EM; returns per-item kernel
  probabilities (guess/slip for DINA/DINO), profile prevalences, and
  per-person marginal attribute-mastery posteriors. Validated against
  CDM::din and CDM::gdina (logLik to ~1e-4 on simulated G-DINA; classic
  fraction-subtraction DINA matched on logLik and guess/slip). The EM
  loop runs in compiled C++ on BLAS level-3 kernels with SQUAREM
  acceleration (Varadhan & Roland 2008) - the 256-profile
  fraction-subtraction DINA fits in ~0.5s per start.
* `fit_lca()` gains order-restricted latent class models. Total orders
  (`ordering = "increasing"`; Croon's 1990 ordered LCM) constrain every
  binary item probability and gaussian indicator mean to be nondecreasing
  across classes. Partial orders (`ordering = list(c(1, 2), c(1, 3),
  c(2, 4), c(3, 4))` etc.) constrain only the specified class pairs,
  leaving unconnected classes incomparable - e.g. attribute-profile
  lattices with incomparable intermediate classes. The constrained M-step
  is a weighted isotonic regression over the class poset
  (pool-adjacent-violators on a chain, Dykstra's projection algorithm on
  a general DAG), so estimation remains closed-form EM with safeguarded
  Ramsay acceleration (extrapolations are projected back into the
  constraint set). A total order removes label switching by construction.
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

# gllammr 1.1.0 and earlier

* Initial development releases: Gaussian/binomial/Poisson GLMMs, ordinal
  and multinomial models, dichotomous and polytomous IRT (Rasch, 2PL, 3PL,
  GRM, PCM, GPCM, NRM), DIF analysis, explanatory IRT, latent class
  analysis, marginal predictions, and frequency/probability weights.
