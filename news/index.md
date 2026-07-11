# Changelog

## gllammr 1.2.0

### Post-1.2.0 development

- **Test/validation robustness against reference-package numerical
  failures.** The cross-walk reference fits that back several IRT/DIF
  tests and validation cases (`rasch_lsat`, `eirt_verbagg`,
  `dif_irt_glmm`, and the `dif_irt` glmer cross-walk) now fit `lme4`
  with the `bobyqa` optimizer, which avoids the “Downdated VtV is not
  positive definite” breakdown the default `nloptwrap`/Nelder-Mead PWRSS
  path can hit on some BLAS/Matrix builds (seen on Linux CI runners with
  Matrix 1.7-5). Estimates are unchanged within every test’s existing
  tolerance. As a belt-and-braces fallback, when a *reference-package*
  fit still errors on a platform the test now `skip()`s with the reason
  and the validation harness records the case as skipped (`pass = NA`)
  rather than failed; a genuine gllammr-side error is still reported and
  fails the suite. No tolerance was weakened.

- **Category-probability predictions (`predict(fit, type = "probs")`)
  now work for Laplace-fitted polytomous IRT models, not just EM.**
  Previously the Laplace path (`fit_irt(..., method = "laplace")`, and
  the `se = TRUE` default via `gllamm(family = irt())`) never assigned
  the `gllamm_irt_poly` S3 class, so the `predict.gllamm_irt_poly`
  method was unreachable and `type = "probs"`/`"expected"` errored.
  Laplace polytomous fits (PCM, GPCM, GRM, NRM), including multi-level
  ones, now carry that class and return the same category-probability
  structure as EM fits.

- **Behavior change:
  [`fit_irt()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_irt.md)
  (and `gllamm(family = irt())`) now compute standard errors by default
  (`se = TRUE`),** aligning the IRT fitter with
  [`fit_sem()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_sem.md)/[`sem()`](https://drjoshmcgrane.github.io/gllammr/reference/sem.md)
  and standard R practice. Because SEs require the Laplace path,
  `method = "auto"` now resolves to `"laplace"` under the default; pass
  `se = FALSE` to skip SE computation and recover the previous fast
  Bock-Aitkin EM default for single-level models. `method = "auto"`
  still selects `"em"` whenever non-integer person weights (an EM-only
  feature) are supplied, in which case the default `se = TRUE` is
  skipped silently (an explicit `se = TRUE` warns, as before). The
  [`irt()`](https://drjoshmcgrane.github.io/gllammr/reference/irt.md)
  family constructor gains a matching `se` argument that
  [`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
  passes through.

- **CI: coverage, lint, pkgdown, and scheduled-validation workflows.**
  Added `.github/workflows/test-coverage.yaml` (covr + Codecov, with a
  build-artifact fallback), `.github/workflows/lint.yaml` (non-blocking
  `lintr::lint_package()` against a permissive `.lintr` baseline),
  `.github/workflows/pkgdown.yaml` (builds and deploys the pkgdown site,
  with a `_pkgdown.yml` reference index covering every exported and
  internal documented topic), and `.github/workflows/validation.yaml`
  (weekly cron plus manual dispatch running
  `validation/run_validation.R` against the full Suggests reference
  packages, failing the job if any cross-package check regresses).

- **Bug fix:
  [`fit()`](https://drjoshmcgrane.github.io/gllammr/reference/fit.md) on
  LCA models.**
  [`fit.gllamm_lca()`](https://drjoshmcgrane.github.io/gllammr/reference/fit.gllamm_lca.md)
  always errored (`non-numeric argument to mathematical function`)
  because it read `object$n_classes`/`object$posterior_probs`, fields
  that don’t exist on a `gllamm_lca` fit (the real fields are
  `nclass`/`posterior`). Entropy, class proportions, and average
  posterior probabilities (APPA) now compute correctly; caught by a new
  test exercising
  [`fit()`](https://drjoshmcgrane.github.io/gllammr/reference/fit.md) on
  an actual fitted LCA model.

- **Test-suite hardening.** Seven
  `tryCatch(gllamm(...), error = ... skip(...))` wrappers in the
  marginal-prediction tests were dead scaffolding from a since-fixed bug
  (verified never to fire) and now call
  [`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
  directly, so a future regression fails the suite instead of silently
  skipping. The two `step_formula`/`step_data` validation guards in
  [`fit_eirt()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_eirt.md)
  (unsupported model; missing `step_data`) and the shared step/threshold
  “at least 3 categories” guard are now covered by `expect_error()`
  tests, and `step_formula` reachability through the
  `gllamm(family = eirt(...))` unified interface is locked down by an
  equivalence test against
  [`fit_eirt()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_eirt.md)
  directly. Former empty placeholder tests for
  [`fit()`](https://drjoshmcgrane.github.io/gllammr/reference/fit.md) on
  LCA/IRT models and IRT/LCA S3
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods are
  now real fits exercising the documented structure. Package version
  bumped to 1.2.0.9000 (development version) ahead of the next release.

- **Packaging hygiene.** `.Rbuildignore` now excludes the built tarball
  and `R CMD check` directory under either package-name casing
  (`GLLAMMR_*.tar.gz`/`gllammr_*.tar.gz`, `GLLAMMR.Rcheck`/
  `gllammr.Rcheck`); no compiled objects, tarballs, or check directories
  were tracked in git.

- **Vectorized marginal-prediction integrator.** The Monte Carlo
  integrator behind `predict(type = "marginal")` no longer loops in R
  over the `n_sim` draws. All population-level random-effects draws are
  generated up front and the inverse-link probabilities reduced
  column-wise in a couple of matrix operations, with a memory guard that
  processes the draws in column blocks once `n_obs * n_sim` exceeds 5e7.
  Draws are produced in the exact random-number order the former
  per-replicate loop consumed, so a fixed seed reproduces prior results
  bit-for-bit (a new equivalence test locks this to 1e-10 end-to-end and
  1e-12 against the old Welford integrator). Measured speedups range
  from ~1.2x (large `n_obs x n_sim`, bounded by the inverse-link
  evaluations) to ~3x when the draw count dominates the observation
  count; random-slope fits gain most because the covariance is now
  Cholesky-factored once rather than on every replicate. Cook’s-distance
  and DIF purification refit loops were audited for cacheable model
  matrices: their invariant design matrices are already hoisted out of
  the loops and the loops are refit-bound (`glm.fit`/TMB), so no further
  change was warranted.

- **Internal TMB engine consolidation.** The legacy v1 TMB interface
  (`fit_tmb_gllamm()`, `R/tmb_interface.R`) has been removed; it was
  unreachable dead code
  ([`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
  always dispatched to the v2 engine when available, which is always).
  [`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
  now calls
  [`fit_tmb_gllamm_v2()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_tmb_gllamm_v2.md)
  unconditionally for the non-AGHQ path. No user-facing change.

- **Fit-time robustness diagnostics.** Shared internal helpers now give
  consistent, informative diagnostics across every fitter:

  - Standard errors are validated after
    [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html): a
    warning is emitted when the Hessian is not positive definite
    (standard errors unreliable, e.g. an over-parameterized model or a
    variance component near zero), and each fitted object carries an
    `se_ok` flag.
  - [`coef()`](https://rdrr.io/r/stats/coef.html),
    [`vcov()`](https://rdrr.io/r/stats/vcov.html), and
    [`predict()`](https://rdrr.io/r/stats/predict.html) warn when called
    on a model that did not converge (“estimates may be unreliable”).
  - Optimization failures (an nlminb error or a non-finite objective)
    now raise a single informative error suggesting different starting
    values or a simpler model.
  - The marginal-prediction sampler falls back to the nearest
    positive-definite matrix
    ([`Matrix::nearPD`](https://rdrr.io/pkg/Matrix/man/nearPD.html))
    when the random-effects covariance is not positive definite, instead
    of erroring.
  - [`solve()`](https://rdrr.io/r/base/solve.html) is guarded in the
    sandwich (robust) variance estimator and in Cook’s-distance
    diagnostics: a singular matrix now warns and returns an `NA`-filled
    result of the expected shape rather than aborting.

- **Step-level predictors in explanatory IRT.**
  [`fit_eirt()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_eirt.md)
  (and the
  [`eirt()`](https://drjoshmcgrane.github.io/gllammr/reference/eirt.md)
  family) gain `step_formula`/`step_data` for PCM/GPCM: covariates that
  vary WITHIN an item across its steps, each with a single common
  coefficient (delta_im = b_i + xi_0m + sum_k eta_k x_imk + e_im) - the
  step-property models of Kim & Wilson (2019). Combines freely with
  difficulty_formula (item level) and threshold_formula (item properties
  with step-specific effects); a step covariate constant within items
  reproduces the equivalent item-covariate fit exactly. A combined
  item/threshold/step rank check rejects cross-level collinearity by
  name.

- **LPCM identification fix.** The threshold regression previously
  shared a flat ridge with the difficulty regression whenever a column
  (including the intercept) appeared in both designs - all standard
  errors NaN. Threshold coefficients are now sum-to-zero deviations
  across an item’s thresholds (matching the package’s PCM/GRM
  conventions), with item-level main effects carried by
  difficulty_formula: both location effects and step deviations are
  recovered with finite SEs.

- Newdata-prediction audit across all formula-based classes (predictions
  on newdata equal to in-sample fits for seen groups; unseen groups get
  population-level predictions). **Bug fixes:** binomial fits from the
  dedicated single-term path ignored the estimated group effects on
  newdata (silently fixed-effects-only) and refused
  [`ranef()`](https://drjoshmcgrane.github.io/gllammr/reference/ranef.md) -
  BLUPs are now stored on the fit; NPML fits had no working
  [`predict()`](https://rdrr.io/r/stats/predict.html) at all
  (non-conformable error on newdata, NULL in-sample) -
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) now use
  posterior-mean mass-point intercepts per group, with the prior mean
  for unseen groups.

- Sandwich-SE audit. Cluster-robust fixed-effect covariances match
  `clubSandwich::vcovCR(type = "CR0")` on the same ML linear mixed model
  (within 1%), and a 200-replication Monte Carlo under heteroskedastic
  misspecification confirms the sandwich SE tracks the true sampling SD
  (0.084 vs 0.090) where the model-based SE is badly anticonservative
  (0.057). clubSandwich added to Suggests.

- Mixed-response audit. Verified: single-outcome
  [`fit_mixed()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_mixed.md)
  equals the corresponding GLMM exactly; with genuinely shared random
  effects the joint fit dominates the sum of separate fits; with
  independent outcomes the forced loading-1 coupling degrades fit as
  expected (outcome-specific loadings remain documented as unsupported).
  **New methods:** [`predict()`](https://rdrr.io/r/stats/predict.html)
  and [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) for
  `gllamm_mixed` (per-outcome conditional predictions, link or response
  scale, newdata supported) -
  [`predict()`](https://rdrr.io/r/stats/predict.html) previously
  returned NULL silently.

- GLMM level-2 weight audit (follow-up to the EIRT weight finding).
  Level-2 (group) weights under Laplace previously scaled each group’s
  likelihood-plus-prior contribution - approximate in the variance
  components and unbounded in principle (every weighted group
  contributes -(w-1)\*log(sigma_u)). Integer level-2 frequency weights
  are now implemented by exact whole-group replication before fitting
  (weighted fits are identical to group-duplicated fits, verified for
  gaussian/poisson/binomial; replication happens before listwise
  deletion, so incomplete rows are handled exactly as duplicated data
  would be); non-integer level-2 weights under Laplace are rejected with
  guidance to `integration = aghq(k)`, which weights each group’s log
  marginal likelihood outside the integral and is exact for arbitrary
  weights (verified identical to duplication). Binomial fits with list
  weights now route to the general engine instead of failing in
  [`fit_binomial()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_binomial.md);
  ordinal models reject level-specific weights with a clear message;
  survival models likewise reject list weights clearly (previously an
  obscure coercion error).

- EIRT deep audit. **Saturated identities verified exactly:** an EIRT
  model with `difficulty_formula = ~ item` (and, where applicable,
  `discrimination_formula = ~ item`) reproduces the descriptive
  [`fit_irt()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_irt.md)
  likelihood for Rasch, 2PL, PCM, GPCM, and GRM; multilevel Rasch EIRT
  matches [`lme4::glmer`](https://rdrr.io/pkg/lme4/man/glmer.html) on
  the crossed person/school cross-walk (logLik, variance components, and
  standard errors). **Weight semantics fixed:** person-level weights on
  Laplace-based fits
  ([`fit_eirt()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_eirt.md),
  and `fit_irt(method = "laplace")`) previously weighted only the
  response terms, leaving the ability prior unweighted - not a valid
  frequency- or sampling-weight likelihood. Weighting the prior instead
  is degenerate under the Laplace approximation (each weighted person
  contributes -(w-1)\*log(sigma_theta), so the objective is unbounded
  and sigma_theta collapses to 0). Integer frequency weights are now
  implemented by exact replication of weighted persons, making weighted
  fits identical to duplicated-data fits (verified to 1e-6, single-level
  and multilevel, dichotomous and polytomous); non-integer person
  weights on Laplace paths are rejected with guidance to
  `method = "em"`, which weights each person’s log marginal likelihood
  directly (already exact). **Multilevel method fixes:**
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) on multilevel
  IRT/EIRT fits drew only the person deviation and ignored group random
  effects (school-level variance was ~4x too small); marginal
  predictions likewise integrated over `sigma_theta` only - both now use
  the full latent structure. **New guards:** rank-deficient
  `difficulty_formula`/`discrimination_formula`/`threshold_formula`
  designs error early naming the aliased columns (previously NaN
  standard errors); `threshold_formula` with GRM (or dichotomous models)
  errors instead of being silently ignored.

- CI: `nnet` added to Suggests (used by audit tests);
  `-Wno-array-bounds` moved out of the shipped `Makevars.win` (it
  triggered R CMD check’s non-portable-flags warning) onto the Windows
  CI runner only.

- Package-wide audit (every S3 method probed on every fit variant;
  under-validated likelihoods checked externally). **Bug fixes:** on
  crossed/multi-term GLMMs,
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) drew random
  effects for the first term only and marginal predictions integrated
  over only the first term’s variance (both silently wrong; conditional
  newdata predictions also applied only the first term’s BLUPs) - all
  three now loop over every term;
  [`summary()`](https://rdrr.io/r/base/summary.html) recycled mismatched
  covariance matrices into wrong standard errors for survival/rank/NPML
  fits (SEs are now shown only when a matching vcov exists);
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) errored on SEM fits
  instead of returning the parameter covariance;
  [`parse_formula()`](https://drjoshmcgrane.github.io/gllammr/reference/parse_formula.md)
  mangled call-type responses like `Surv(t, d)`. **New methods:**
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) for multinomial,
  survival (uncensored event times from the fitted hazard), IRT
  (parametric bootstrap, dichotomous + polytomous, EM and Laplace fits),
  EIRT, LCA (all indicator types), CDM, SEM (implied-distribution
  draws), NPML (mass-point draws), and mixed responses;
  [`predict()`](https://rdrr.io/r/stats/predict.html) for polytomous IRT
  (category probabilities and expected scores); multinomial fits now
  store their random effects so
  [`ranef()`](https://drjoshmcgrane.github.io/gllammr/reference/ranef.md)
  works. `icc`/`VarCorr`/`ranef` refusals now point at the right
  accessor for each latent-variable class. **External validations
  added:** multinomial matches
  [`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html) in the
  fixed-effect limit (1e-4); Weibull frailty matches
  [`survival::survreg`](https://rdrr.io/pkg/survival/man/survreg.html)
  in the no-frailty limit (with the AFT-scale parameterization now
  documented: log hazard ratio = shape x beta); rank-ordered logit
  coefficients equal the exploded conditional logit
  ([`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html),
  Plackett-Luce equivalence, 1e-4).

- Ordinal GLMM audit. **Bug fix:** the backward continuation-ratio link
  (`crl_backward`) was not a valid probability model - its category
  probabilities summed to 2\*plogis(tau_max - eta) because the top
  category reused the last cumulative probability; it now uses the
  proper hazard-product form and, like `acl` and `crl_forward`, is
  validated against VGAM (acat/sratio; agreement ~1e-3).
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) now support all
  six ordinal links (previously logit/probit only) through a shared
  category-probability helper mirroring the C++ likelihoods, and both -
  along with marginal predictions - handle crossed/multi-term random
  effects (marginal predictions previously integrated over only the
  first term’s variance; simulate() errored).
  [`icc()`](https://drjoshmcgrane.github.io/gllammr/reference/icc.md)
  gains a proper latent-response-scale branch for cumulative ordinal
  models (it previously fell through to the gaussian formula); threshold
  starting values are computed from raw spacings. Single-term ordinal
  fits run at ordinal::clmm speed parity with identical likelihoods.
  VGAM added to Suggests.

- New general
  [`compare_models()`](https://drjoshmcgrane.github.io/gllammr/reference/compare_models.md):
  a comparison table (logLik, parameter count, AIC/BIC with deltas,
  Akaike weights) for any set of fitted gllammr models, across model
  classes, with an n_obs comparability check;
  [`latent_structure_comparison()`](https://drjoshmcgrane.github.io/gllammr/reference/latent_structure_comparison.md)
  now delegates its table to it.

- Order-restricted LCA optimized: the isotonic poset projection runs in
  C (shared with the CDM engine), and the Ramsay acceleration trigger
  window widened to cover slowly converging EM (step ratio near 1; the
  logLik-revert safeguard makes this safe). MON/IIO/DM fits are 2.7-4.2x
  faster; the located class model now converges where it previously hit
  the iteration cap, and its EM is ~10x faster than fitting the same
  model through the equivalent NPML GLMM route (identical logLik).

- Latent structure analysis (Torres Irribarra & Diakow framework):
  [`fit_lca()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_lca.md)
  gains `item_ordering` (invariant item ordering - item monotonicity
  within classes, via isotonic regression over the item-by-class grid;
  combined with `ordering` this is the double monotonicity model) and
  `structure = "rasch"` (located latent classes / latent class Rasch:
  logit pi_ic = theta_c - delta_i, classes on an interval scale,
  reported sorted by location with `class_locations` and
  `item_difficulties`). The located class model’s likelihood is verified
  against
  [`fit_npml()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_npml.md)
  on the long-format GLMM (the Lindsay-Clogg-Grego equivalence, to 0.02
  logLik). New
  [`latent_structure_comparison()`](https://drjoshmcgrane.github.io/gllammr/reference/latent_structure_comparison.md)
  fits all six models of the framework (UN, MON, IIO, DM, LCR, RM) and
  returns the successive-comparison table for deciding between
  qualitative, ordinal, and quantitative latent structure.

- Vignette suite completed: four new evaluated vignettes -
  `explanatory-irt` (De Boeck & Wilson on the verbal aggression data,
  with the GLMM cross-walk and the Kim & Wilson polytomous extension),
  `cognitive-diagnosis` (Q-matrices, DINA/DINO/G-DINA, monotonicity,
  hierarchies), `sem-models` (CFA, fit indices, MIMIC, FIML,
  standardized solution), and `dif-analysis` (the two-stage
  screen-then-confirm DIF workflow). `latent-class` rewritten with
  evaluated code covering ordered and partially ordered classes;
  `advanced-features` rewritten with evaluated examples of mixed
  responses, frailty survival, rankings, NPML, AGHQ, and sandwich
  standard errors through the unified interface; `getting-started` gains
  the complete model-space table. All thirteen vignettes build in under
  two minutes.
  [`fit_lca()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_lca.md)
  now rejects `nclass = 1` with an informative message.

- The unified
  [`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
  interface now reaches every model class: new family constructors
  [`eirt()`](https://drjoshmcgrane.github.io/gllammr/reference/eirt.md),
  [`sem()`](https://drjoshmcgrane.github.io/gllammr/reference/sem.md),
  [`mixed_response()`](https://drjoshmcgrane.github.io/gllammr/reference/mixed_response.md),
  [`survival_family()`](https://drjoshmcgrane.github.io/gllammr/reference/survival_family.md),
  and
  [`ranking()`](https://drjoshmcgrane.github.io/gllammr/reference/ranking.md),
  plus `integration = npml(k)` for nonparametric (mass-point) latent
  distributions. Also fixed: `family = binomial()` with
  `integration = aghq(k)` or `npml(k)` was silently ignoring the
  integration request (the binomial fast path intercepted dispatch); it
  now routes through the requested engine.

- New
  [`dif_irt()`](https://drjoshmcgrane.github.io/gllammr/reference/dif_irt.md):
  confirmatory model-based DIF (IRT-LR; Thissen et al.) as the companion
  to the `dif_test` screening tests. DIF parameters are
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

- DIF analysis rewritten (`dif_test`): logistic-regression DIF with a
  latent (EAP) or observed-score matching criterion, a formula interface
  for **multiple DIF variables and their interactions**
  (`dif = ~ gender * language`), **iterative purification** of the
  matching criterion (with graceful degradation and a warning when the
  DIF/impact decomposition is unidentified), uniform/nonuniform/joint LR
  tests with correct per-item degrees of freedom, Nagelkerke delta-R2
  effect sizes with the Jodoin-Gierl A/B/C classification, optional
  multiplicity adjustment, anchor-item support, and cumulative-logit
  tests for polytomous items. With score matching the flags reproduce
  [`difR::difLogistic`](https://rdrr.io/pkg/difR/man/difLogistic.html)
  exactly (validation case `dif_logistic`). The previous implementation
  divided one global LR statistic equally across items and confounded
  DIF with impact via separate per-group calibrations;
  `dif_test_with_data` is deprecated and now wraps the new engine.
  `dif_plot` draws model-implied response curves by any DIF variable.

- SEM overhaul (`fit_sem`): exogenous latent variables now covary freely
  (previously silently orthogonal - a plain two-factor CFA was
  misspecified); MIMIC models (structural regressions on observed
  covariates, likelihood-equivalent to lavaan with `fixed.x = FALSE`);
  standard errors for every parameter (numerical observed information +
  delta method); the standard fit-index battery (chisq, CFI, TLI, RMSEA
  with 90% CI, SRMR - matching lavaan to 4 decimals); full-information
  maximum likelihood for missing data (`missing = "fiml"`,
  pattern-based, with an EM-estimated saturated model for fit indices);
  a standardized (std.all) solution; and a real
  [`summary()`](https://rdrr.io/r/base/summary.html) method. The legacy
  Laplace path warns that it treats exogenous factors as orthogonal.

- Package-wide missing-data policy, audited and tested: formula-based
  fitters (GLMM, ordinal, multinomial, survival, NPML, mixed responses)
  now listwise-delete rows with missing values in any model variable -
  with a warning and automatic weight alignment. Previously the response
  and fixed-effect design dropped NA rows while the random-effects
  design and grouping factor kept them, silently misaligning ordinal,
  survival, and NPML fits. Matrix-response latent variable models
  (IRT/EIRT/LCA/CDM) already used all observed responses (MAR) and are
  unchanged; `fit_rank` keeps its partial-ranking semantics for missing
  ranks.

- Multiple (crossed/nested) random-effects terms for ordinal and
  multinomial models: `fit_ordinal(y ~ x + (1 | rater) + (1 | item))`
  works for the logit, probit, adjacent-category, and continuation-ratio
  links (PPO remains single-term), matching
  [`ordinal::clmm`](https://rdrr.io/pkg/ordinal/man/clmm.html) with the
  same crossed structure to 1e-4;
  [`fit_multinomial()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_multinomial.md)
  gains the same layout with the random effects acting as a common
  shifter across non-reference categories. The same
  `(1 | g1) + (1 | g2)` formulas now also work for `family = binomial()`
  through
  [`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md).

- New
  [`cdm()`](https://drjoshmcgrane.github.io/gllammr/reference/cdm.md)
  family constructor and `lca(ordering = ...)` pass-through: cognitive
  diagnosis models and order-restricted latent class models are now
  reachable through the unified
  [`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
  interface.

- New
  [`fit_cdm()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_cdm.md):
  cognitive diagnosis models for binary responses with a Q-matrix -
  saturated G-DINA (default), DINA, and DINO - with monotonicity in the
  attributes enforced by isotonic regression over the reduced-profile
  lattice (`monotone = TRUE`, the default) and optional attribute
  hierarchies (prerequisite relations prune the profile space).
  Closed-form accelerated EM; returns per-item kernel probabilities
  (guess/slip for DINA/DINO), profile prevalences, and per-person
  marginal attribute-mastery posteriors. Validated against CDM::din and
  CDM::gdina (logLik to ~1e-4 on simulated G-DINA; classic
  fraction-subtraction DINA matched on logLik and guess/slip). The EM
  loop runs in compiled C++ on BLAS level-3 kernels with SQUAREM
  acceleration (Varadhan & Roland 2008) - the 256-profile
  fraction-subtraction DINA fits in ~0.5s per start.

- [`fit_lca()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_lca.md)
  gains order-restricted latent class models. Total orders
  (`ordering = "increasing"`; Croon’s 1990 ordered LCM) constrain every
  binary item probability and gaussian indicator mean to be
  nondecreasing across classes. Partial orders
  (`ordering = list(c(1, 2), c(1, 3), c(2, 4), c(3, 4))` etc.) constrain
  only the specified class pairs, leaving unconnected classes
  incomparable - e.g. attribute-profile lattices with incomparable
  intermediate classes. The constrained M-step is a weighted isotonic
  regression over the class poset (pool-adjacent-violators on a chain,
  Dykstra’s projection algorithm on a general DAG), so estimation
  remains closed-form EM with safeguarded Ramsay acceleration
  (extrapolations are projected back into the constraint set). A total
  order removes label switching by construction.

- [`fit_eirt()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_eirt.md)
  identification fixes: the explanatory GRM now expresses its ordered
  thresholds as sum-to-zero deviations around the item location
  (previously a free first threshold per item silently absorbed the
  difficulty regression, and collided with the item residuals under
  `item_residuals = TRUE`); `sigma_theta` is fixed at 1 whenever a
  discrimination level is estimated (2PL/GRM/GPCM - the `fit_irt`
  convention; the two traded off on a flat ridge); unused `step_param`
  cells are mapped off (previously left free, giving singular Hessians
  and NaN standard errors in polytomous EIRT models).

- New validation cases on the verbal aggression data: `eirt_verbagg` (De
  Boeck & Wilson LLTM+error vs the lme4 crossed-effects GLMM) and
  `eirt_verbagg_pcm` (location-explanatory PCM with random item effects
  vs Kim & Wilson 2019 published estimates).

- [`fit_irt()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_irt.md)
  gains a Bock-Aitkin MML-EM estimation path (`method = "em"`), now the
  **default for single-level models**: 20-50x faster than the Laplace
  path, matches mirt to correlation 1.0, and handles short tests where
  joint-Laplace 2PL diverges. The Laplace path remains the default
  whenever multi-level structure (`random`) or standard errors
  (`se = TRUE`) require it. EM abilities are EAP scores.

- `fit_irt(se = FALSE)` is the default (as in mirt); request SEs
  explicitly.

- Single-pass AD tapes in the gaussian/binomial/poisson templates.

### Build and infrastructure

- All TMB C++ templates are now compiled at install time into a single
  shared library (CRAN-compatible). No run-time compilation, no compiled
  artifacts cached in the user’s home directory, and no C++ toolchain
  required after installation.

### Full GLLAMM parity features

- Random slopes `(x | g)` and uncorrelated slopes `(x || g)` with full
  variance-covariance estimation across Gaussian, binomial, Poisson, and
  gamma families.
- Crossed `(1 | g1) + (1 | g2)` and nested `(1 | g1/g2)` random effects
  in the general estimator.
- Multilevel explanatory IRT: person-level random effects
  (`random = ~ (1 | class)` etc.) and person-level random slopes in IRT
  and EIRT models, including nested and crossed structures.
- Gamma family (log link) for positive continuous responses.
- Parametric survival/frailty models
  ([`fit_survival()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_survival.md);
  exponential and Weibull with shared log-normal frailty).
- Joint mixed-response models sharing a random effect across Gaussian,
  binomial, and Poisson outcomes
  ([`fit_mixed()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_mixed.md)).
- Structural equation models with measurement and structural parts
  ([`fit_sem()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_sem.md)).
- Rank-ordered (exploded) logit with optional random coefficients
  ([`fit_rank()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_rank.md)).
- Level-specific survey weights via
  `weights = list(level1 = , level2 = )` (pseudo-likelihood, matching
  Stata GLLAMM `pweight()` behavior).
- Nonparametric maximum likelihood (NPML) mass-point random effects
  ([`fit_npml()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_npml.md)),
  validated against npmlreg.
- Adaptive Gauss-Hermite quadrature via `integration = aghq(k)` as an
  alternative to the default Laplace approximation.
- Cluster-robust (sandwich) standard errors via
  `vcov(fit, type = "sandwich")`.
- Latent class analysis with mixed indicator types (binary, categorical,
  and continuous manifest variables).

### Unified interface

- All model classes are reachable through
  [`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
  family dispatch:
  [`irt()`](https://drjoshmcgrane.github.io/gllammr/reference/irt.md),
  [`lca()`](https://drjoshmcgrane.github.io/gllammr/reference/lca.md),
  [`ordinal()`](https://drjoshmcgrane.github.io/gllammr/reference/ordinal.md),
  [`multinomial()`](https://drjoshmcgrane.github.io/gllammr/reference/multinomial.md),
  and
  [`binomial()`](https://drjoshmcgrane.github.io/gllammr/reference/binomial.md)
  family objects route to the appropriate fitter, e.g.
  `gllamm(resp_matrix, family = irt("2PL"))` or
  `gllamm(y ~ x + (1 | g), data, family = ordinal("probit"))`.
- `predict(fit, type = "marginal")` population-averaged predictions and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) (returning a
  data frame) across model classes.

### Identification and correctness fixes

- Ordinal models: intercept absorbed into thresholds (standard
  identification), fixing previously shifted threshold estimates.
- IRT: latent variance fixed at 1 for 2PL/3PL/GRM/GPCM/NRM so
  discriminations are identified; Rasch frees the latent variance.
- LCA: corrected parameterization of class-membership and item-response
  logits; polytomous coding fixed.
- GRM: threshold orientation aligned with mirt conventions.

### Performance

- 40-70x speedups for standard GLMMs from dedicated TMB templates,
  sparse random-effects design matrices, and better starting values.

### Validation

- New cross-package validation suite,
  [`gllammr_validate()`](https://drjoshmcgrane.github.io/gllammr/reference/gllammr_validate.md):
  49 automated checks of estimates, variance components, and
  log-likelihoods against lme4, glmmTMB, ordinal, mirt, poLCA, npmlreg,
  and lavaan.

## gllammr 1.1.0 and earlier

- Initial development releases: Gaussian/binomial/Poisson GLMMs, ordinal
  and multinomial models, dichotomous and polytomous IRT (Rasch, 2PL,
  3PL, GRM, PCM, GPCM, NRM), DIF analysis, explanatory IRT, latent class
  analysis, marginal predictions, and frequency/probability weights.
