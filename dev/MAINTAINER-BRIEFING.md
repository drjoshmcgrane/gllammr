# gllammr 1.3.0 — Maintainer's Briefing

Prepared to let the maintainer answer CRAN-reviewer and psychometrics-expert
questions about design decisions from the code, not from memory. Every claim
below is derived from reading the source in this tree; file:line references are
given where a reviewer is likely to probe. Where a decision is genuinely
debatable, it is flagged as such rather than defended.

This file is developer-only. It is excluded from the build via `^dev$` in
`.Rbuildignore` and never ships in the tarball.

---

## 1. Architecture in one page

The package is a thin, opinionated R layer over a single TMB shared object plus
three hand-written C/C++ numerical cores. The data flow is:

**(a) Family / spec constructors** (`R/families.R`). Every model class is
reachable through a family object: `ordinal()`, `binomial()`, `irt()`,
`eirt()`, `lca()`, `cdm()`, `sem()`, `mixed_response()`, `ranking()`,
`survival_family()`, `multinomial()`. Each returns a small S3 list carrying its
options and a class tag (`irt_family`, `sem_family`, …). `binomial()` is built
on top of `stats::binomial()` so the object still works when `gllammr::binomial`
masks `stats::binomial` (`families.R:188-201`); the S4 old-class registrations
at `families.R:10-14` keep these usable by S4 consumers such as lme4 slots.
Integration is chosen orthogonally with `aghq(k)` / `npml(k)` (`R/aghq.R`).

**(b) `gllamm()` dispatch** (`R/gllamm.R`). The single user-facing entry point
inspects `inherits(family, "<type>_family")` and routes to a specialised
`fit_*` engine (`gllamm.R:178-307`). Crucially, the matrix-response families
(`irt`, `eirt`, `lca`, `cdm`) and the data-frame family (`sem`) dispatch
*before* formula validation, because for those the first positional argument is
not a formula — it is the response matrix or data frame. All of the `fit_*`
engines are also exported and callable directly, so `gllamm(Y, family = irt("2PL"))`
and `fit_irt(Y, model = "2PL")` are two doors to the same code.

**(c) Computational back ends.** Three exist, and which one serves a family is
the single most important thing to have straight:

| Family / model | Engine | Back end |
|---|---|---|
| GLMM (gaussian/binomial/poisson/Gamma), ordinal, multinomial, IRT (default), EIRT, SEM, survival frailty, rank, mixed-response, LCA (TMB path), NPML | R driver → `TMB::MakeADFun(DLL="gllammr")` | Laplace via TMB, one C++ translation unit |
| Two-level single-random-intercept GLMM with `integration = aghq(k)` | `fit_tmb_gllamm_aghq` (`R/aghq.R`) | Adaptive Gauss–Hermite, TMB `glmm_aghq` template |
| Polytomous IRT under `method = "em"` (GRM/PCM/GPCM/NRM) | `fit_irt_em` → `.Call(C_em_poly)` | Bock–Aitkin MML-EM, `src/em_poly.cpp` |
| CDM (G-DINA/DINA/DINO) | `fit_cdm` → `.Call(C_em_cdm)` | MML-EM on BLAS kernels, `src/em_cdm.cpp` |
| LCA (default `method = "em"`) | `fit_lca_em` (`R/lca_em.R`) | Closed-form R EM (poLCA algorithm), order restrictions via `src/isotonic.cpp` / PAVA |

**(d) TMB dispatch layer.** `src/gllammr.cpp` is a *single* translation unit.
It `#include`s the 22 model headers under `src/include/` and defines one
`objective_function::operator()`, which switches on `DATA_STRING(model_name)` to
the right `gllamm_<name>()` template (`gllammr.cpp:30-60`). Exactly one
`operator()` is permitted per DLL, hence the string-dispatch design. The package
provides its own `R_init_gllammr` (`gllammr.cpp:79-85`) rather than
`-DTMB_LIB_INIT` so that the plain-C EM entry points register alongside TMB's
routines in the same `.so`: `C_em_poly` (9 args), `C_em_cdm` (10),
`C_isotonic_poset` (3) at `gllammr.cpp:71-76`. This is why the ~14 MB `libs`
size is unavoidable — every template shares one `TMB.hpp` include.

**(e) Post-estimation S3 layer.** Each engine returns a classed list
(`gllamm`, `gllamm_irt`, `gllamm_irt_poly`, `gllamm_lca`, `gllamm_cdm`,
`gllamm_sem`, …) with `coef`/`summary`/`predict`/`ranef`/`VarCorr`/`simulate`
methods spread across `R/predict_*.R`, `R/*_methods.R`, `R/plot_*.R`. 141
`@export`s in total.

**File map a maintainer needs:** dispatch `R/gllamm.R`; specs `R/families.R`,
`R/aghq.R`; robustness helpers `R/fit_checks.R`; dichotomous/polytomous IRT
`R/irt.R` + `R/irt_em.R` + `src/em_poly.cpp`; EIRT `R/eirt.R` +
`src/include/gllamm_eirt.hpp`; DIF `R/dif.R` + `R/dif_irt.R` +
`src/include/gllamm_irt_dif.hpp`; LCA `R/latent_class.R` + `R/lca_em.R`; CDM
`R/cdm.R` + `src/em_cdm.cpp`; SEM `R/sem.R` + `R/sem_ml.R`; validation
`R/validation.R` + `validation/RESULTS.md`; tests `tests/testthat/` (57 files).

---

## 2. Estimation design decisions

**Laplace via TMB is the default and the spine of the package.** The random
effects are integrated out by TMB's automatic Laplace approximation; parameter
optimisation of the resulting marginal likelihood is done by `stats::nlminb`.
The choice is deliberate: it is exact for Gaussian responses, fast, and — this
is the point that makes the validation suite tight — *bit-for-bit the same
approximation lme4 uses at `nAGQ = 1` and `ordinal::clmm` uses*. That is why
those references are expected to agree to numerical precision rather than merely
within tolerance (see §4).

**`nlminb` and the guard rails.** Every optimiser call goes through
`safe_nlminb()` (`R/fit_checks.R:83-101`), which turns an nlminb error or a
non-finite objective into one informative error instead of an opaque failure.
Standard errors come from `TMB::sdreport`, and every engine funnels the result
through `check_sdreport()` (`R/fit_checks.R:27-46`). `se_ok` is `TRUE` **only
when** `sdreport` succeeded *and* `sdr$pdHess` is `TRUE`; a non-positive-definite
Hessian yields `se_ok = FALSE` and a warning that the model may be
over-parameterised or a variance component near zero. These two helpers are
called from every fitter (grep `check_sdreport` across `R/`: aghq, binomial,
dif_irt, eirt, irt, latent_class, npml, ordinal all use it), so the SE-trust
signal is uniform across families.

**Adaptive Gauss–Hermite quadrature (`R/aghq.R`) is opt-in.** `aghq(k)` selects
the `glmm_aghq` template through `fit_tmb_gllamm_aghq`. It is deliberately
narrow: two-level, a *single* random intercept, family in {gaussian, binomial,
poisson} (`aghq.R:71-79`). The driver implements classic adaptive quadrature —
it alternates parameter optimisation with per-group re-centring on the
posterior mode and re-scaling by the posterior curvature (`aghq.R:151-187`),
iterating until the adaptation points move less than `adapt_tol`. It matters
exactly where Laplace is weakest: small clusters with a large random-effect
variance. The validation case `aghq_binomial` is constructed for precisely that
regime (6 observations/group, `sigma_u = 2`) and checks agreement with
`glmer(nAGQ = 15)` to `logLik` tolerance `0.05` (`validation.R:611-642`).
Laplace is `aghq(1)`; the documentation says so (`aghq.R:7`).

**Bock–Aitkin MML-EM cores exist *alongside* Laplace, not instead of it.** Two
reasons.

1. *Numerical robustness for item models on short tests.* The C++ polytomous
   core (`src/em_poly.cpp`) runs the whole EM loop in C++ with fixed
   Gauss–Hermite quadrature: the E-step accumulates person-by-node posteriors
   and expected category counts; M-steps take damped Newton steps per item with
   finite-difference derivatives, monitoring the true marginal log-likelihood
   (`em_poly.cpp:1-15`). Parameterisations mirror `R/irt_em.R` exactly —
   GRM `[tau_1, log-spacings (K-2), log a]`, PCM `[delta_1..delta_{K-1}]`,
   GPCM adds `log a`, NRM `[c_1..c_{K-1}, log a]`. This path handles the 5-item
   LSAT 2PL that diverges under joint Laplace (see §4).

2. *Non-integer weights are EM-only, and the reason is a real identification
   pathology, not laziness.* Read the comments at `R/irt.R:173-181` and
   `R/gllamm.R:370-378`. Under Laplace, integer frequency weights are
   implemented by **exact replication** of the weighted person/group — each copy
   gets its own random effect — which reproduces a duplicated-data fit exactly.
   The tempting shortcut of instead *scaling* each unit's
   (likelihood + log-prior) contribution by its weight is only approximate under
   the Laplace normalisation, and worse, its objective is **unbounded**: a unit
   with weight `w` contributes `-(w - 1)·log(sigma_theta)` (equivalently
   `log(sigma_u)` for GLMMs), so as the variance component is driven to zero the
   penalised objective runs to `+∞` in the direction of collapse. Scaling would
   therefore drive `sigma_theta → 0`. The EM path sidesteps this because it
   weights each person's *log marginal* likelihood (with the random effect
   already integrated out), which is exact for arbitrary non-negative weights.
   The same argument is why `aghq` can pass level-2 weights straight through
   (`gllamm.R:379-381`): it too weights the log marginal outside the integral.

**NPML mass points** (`npml(k)`, `R/npml.R`) replace the normal latent
distribution with `k` estimated locations and masses (Aitkin 1999), validated
against `npmlreg::allvc` (`validation.R:573-608`). Level-specific weights are
rejected under NPML (`gllamm.R:436-439`).

**`method = "auto"` resolution in `fit_irt`** (`R/irt.R:133-156`). The rule,
verbatim in logic: if a `random` (multilevel) formula is present → `laplace`;
else if non-integer weights are present → `em`; else if `se = TRUE` (the
default) → `laplace`; else → `em`. The `se = TRUE` **default** is the tail that
wags this dog, and it is intentional: standard errors require the Laplace path
(EM does not yet produce them), and the package-wide convention is that fits
report SEs unless asked not to. So the default user experience — `fit_irt(Y)`
or `gllamm(Y, family = irt())` — is Laplace with SEs. `se = FALSE` buys the
faster EM path for single-level models (roughly halving fit time on large
person samples). If a user explicitly asks for `se = TRUE` on the EM path it is
declined with a warning (`irt.R:148-155`); the default `se = TRUE` is dropped
silently when `auto` has already routed to EM for another reason.

---

## 3. Identification constraints, per family

This is the section a psychometrics reviewer will press on. Each constraint
below is read out of the actual template or engine.

**Rasch / 2PL / 3PL scaling** (`src/include/gllamm_irt.hpp`). Abilities carry a
`theta_p ~ N(0, sigma_theta^2)` prior (`gllamm_irt.hpp:40-42`) and
`sigma_theta = exp(log_sigma_theta)` is a **free** parameter (`:32`, `:90`). All
item difficulties are free; there is no anchored item and no sum-to-zero on the
difficulties. The metric is therefore fixed the GLLAMM/lme4 way: **ability mean
at 0 (via the prior), ability SD estimated**, item locations free. This is the
key contrast with mirt/ltm, which fix `theta ~ N(0,1)` and free the metric
through the discriminations; the GRM validation case explicitly converts between
the two conventions before comparing (`b* = b/sigma_theta`, `a* = a·sigma_theta`,
`validation.R:385-401`). For 2PL the discriminations are likewise free with no
prior. **Debatable point to own:** in the *DIF* template only, `log_sigma_theta`
is mapped to 0 (i.e. `sigma_theta = 1`) for the 2PL variant while it stays free
for Rasch (`gllamm_irt_dif.hpp:38`) — a deliberate choice to make the 2PL DIF
model identified against impact, but an asymmetry worth being able to explain.
**3PL is the least-constrained model in the package:** the guessing parameter
enters as `prob = c + (1 - c)·invlogit(eta)` with `c = guessing(item)` used
directly on the probability scale, with no logit/bound and no prior
(`gllamm_irt.hpp:62-70`). Combined with the Laplace-on-short-tests issue below,
3PL is the model to be most cautious about; it is offered but not
heavily validated.

**GRM ordered thresholds** (`src/include/gllamm_eirt.hpp:269-303`, and the
descriptive `em_poly.cpp` parameterisation). Ordered thresholds are built as
sum-to-zero deviations around the item location `b_i`:
`u_1 = 0`, `u_m = u_{m-1} + exp(step_param(item, m-2))` for `m ≥ 2`, then
`tau_m = b_i + (u_m − mean(u))`. The `exp()` guarantees strictly increasing
thresholds (order constraint); centring on `mean(u)` makes the spacings
orthogonal to the location so the difficulty regression is identified — a free
first threshold per item would absorb any change in `gamma`
(`gllamm_eirt.hpp:270-275`). This uses `K−2` free spacing parameters per item,
the same budget as PCM/GPCM.

**PCM / GPCM step sum-to-zero** (`gllamm_eirt.hpp:71-72, 305-311`). Step
deviations `s_{im}` satisfy a sum-to-zero constraint with `K−2` free deviations
per item and the last set to `s_{i,K-2} = −(s_{i,0} + … + s_{i,K-3})`. The item
location `b_i` carries the overall difficulty (`delta_{im} = b_i + s_{im}`, the
MFRM two-fold decomposition), and GPCM adds an item-specific discrimination. In
the standalone EM core the PCM runs internally as a GPCM with a shared slope on
the `N(0,1)` scale and is back-transformed in R (the free-ability-SD
equivalence, `em_poly.cpp:13-15`).

**LPCM threshold regression and its ridge** (`gllamm_eirt.hpp:75-81`). This is
the subtle one. `xi_dev` holds free step-deviation weights of shape
`[p_thresh × (K−2)]`, and the full `xi` row is constrained to **sum to zero
across thresholds** (last column = −rowsum). The stated reason is exactly right
and worth quoting to a reviewer: threshold covariates supply step-specific
*deviations*, while item-level main effects belong in `difficulty_formula`;
without the sum-to-zero constraint "the difficulty and threshold regressions
share a flat ridge whenever a column (including the intercept) appears in both
designs." In other words the constraint removes an aliasing between the location
regression and the threshold regression, not a substantive restriction on
the model.

**Step-level predictors** (`gllamm_eirt.hpp:32-33, 83-84`). `eta_step` are
*common* coefficients (shared across items) on step-level covariates `W_step`,
whose layout is one row per item-step cell, `row = item·(max_categories−1) +
(m−1)`. Step residuals `e_step` around the step regression carry a
`N(0, sigma_e_step^2)` prior (`:132-140`).

**EIRT item residuals — LLTM vs LLTM+error** (`gllamm_eirt.hpp:41, 120-130`).
With `item_residuals = 0` the model is a pure LLTM: item parameters are exact
functions of the covariates. With `item_residuals = 1` (the default) it is
LLTM-plus-error: `epsilon_b ~ N(0, sigma_epsilon_b^2)` residuals sit around the
difficulty regression, and (for models that read discrimination)
`epsilon_a ~ N(0, sigma_epsilon_a^2)` around the log-discrimination regression.
A neat correctness anchor: the epsilon_a prior is *skipped* when the model never
reads discrimination, because evaluating a prior on a fixed value only shifts
the log-likelihood by a constant (`:120-129`). The **saturated** EIRT — a
difficulty formula with one coefficient per item — is by construction identical
to descriptive IRT, and that identity is a regression test (§4).

**DIF model** (`src/include/gllamm_irt_dif.hpp`, header comment `:1-12`). This
is the IRT-LR DIF formulation (Thissen–Steinberg–Wainer / De Boeck–Wilson).
Impact is a **latent regression with no intercept**: `theta_p ~ N(z_p'gamma,
sigma^2)` so the reference profile has latent mean 0 (`:5-7, 49-53`). Uniform
DIF adds `z_p'delta_i` to the logit; nonuniform (2PL) DIF scales the
discrimination by `exp(z_p'kappa_i)`. **Anchors identify DIF against impact**:
items flagged `dif_item = -1` carry no DIF parameters (`:30-31`) and pin the
metric so that DIF and group-mean-ability differences are separable. Without
anchors these are confounded; the design forces the user to declare them.

**Multilevel padding-cell trick** (`gllamm_eirt.hpp:149-154`). The group-level
random-effect matrix `u_random` is stored as `[max_n_groups × n_random_effects]`,
so levels with fewer groups than the maximum have unused trailing cells. Those
padding cells are given a `N(0,1)` prior. This is a deliberate device: it keeps
the Laplace Hessian positive definite (an unpenalised free cell would give a
zero-curvature direction) **without changing the marginal likelihood**, since
the cells enter no observation's linear predictor. It is the kind of thing a
sharp reviewer might query as "phantom parameters with priors"; the honest
answer is that they are inert placeholders and the prior is a numerical
regulariser, not a modelling assumption.

**LCA / CDM identification** (`R/latent_class.R`, `R/lca_em.R`,
`src/em_cdm.cpp`, `src/isotonic.cpp`). Plain LCA is a finite mixture and suffers
label switching; the package resolves it structurally rather than by
post-hoc relabelling. `ordering = "increasing"` (Croon 1990) imposes a **total
order** on the classes and "resolves label switching by construction"
(`latent_class.R:31`); a partial order (a list/matrix of class-index pairs)
resolves it up to the poset's automorphisms. The order-restricted M-step is a
weighted isotonic regression — PAVA for a chain (`lca_em.R` `.pava_weighted`),
Dykstra's cyclic projection over a general DAG for a poset
(`src/isotonic.cpp`, `C_isotonic_poset`). An `item_ordering` argument (Croon
1991 invariant item ordering) can additionally constrain item difficulties, and
when both are set the classes are reported *sorted by location* on a shared
interval scale (`latent_class.R:44-50`). CDM uses `src/em_cdm.cpp` for
G-DINA/DINA/DINO: a BLAS-level-3 EM (the E-step is a `dgemm`, the M-step
closed-form weighted proportions per reduced-attribute-profile group),
isotonically projected over the profile lattice by Dykstra when `monotone`
(default), and SQUAREM-accelerated (`em_cdm.cpp:1-21`). Monotonicity is what
identifies which mastery pattern is "higher". The cross-walk that one-attribute
G-DINA equals a 2-class LCA to `1e-7` (§4) is the sharpest single check that the
two independent EM cores agree.

**SEM scaling** (`R/sem.R`, `R/sem_ml.R`). Marker-variable identification: the
first indicator's loading on each latent is fixed at 1 (`sem.R:4-5`), encoded in
`lambda_pattern` with the code `2 = fixed-at-1`, `1 = free`, `0 = zero`
(`sem.R:120-127`, applied at `sem.R:292-293`). Covariate pseudo-indicators are
fixed likewise. Estimation is covariance-based ML (Wishart/FIML likelihood,
`sem_ml.R:1-4`) on the sample covariance and mean vector; a `standardized`
(std.all) solution is derived afterward (`sem.R:176-198`). FIML for missing data
is `missing = "fiml"`, validated against lavaan (§4).

---

## 4. Correctness evidence

**Cross-package validation harness** (`R/validation.R`, exported as
`gllammr_validate()`). Twenty-one canonical-dataset cases plus a four-case
large-scale tier (n in the tens of thousands, long item batteries — the sizes
where quadrature grids and tolerances fail silently), listed at
`validation.R:41-56`. The design principle is explicit in the roxygen
(`validation.R:1-7`) and in each case's tolerance: references that use the
**same Laplace approximation** are expected to agree to numerical precision —
lme4 at `nAGQ = 1` (gaussian sleepstudy `1e-4` on betas, binomial toenail,
poisson grouseticks, survival-as-Poisson, EIRT-as-crossed-glmer, IRT-LR DIF),
`ordinal::clmm` (`1e-2`) — while references that use **different integration**
agree only within tolerance — mirt/ltm EM-quadrature for 2PL/GRM (correlations
`> 0.995`, mean-abs-difference bands), poLCA, CDM, npmlreg, lavaan, difR, VGAM,
nnet. Reference fits are wrapped in `.reference_fit()` (`validation.R:112-127`)
so that a numerical breakdown *inside a reference package* (e.g. lme4's
"Downdated VtV is not positive definite" on some BLAS/Matrix builds) is recorded
as a skip (`pass = NA`), never as a gllammr failure; gllammr's own fit is never
wrapped, so a genuine regression still fails.

**Internal cross-walk identities** (in `tests/testthat/`, gllammr-vs-gllammr, no
external package). These are the strongest evidence because they are exact by
construction:

- **One-attribute CDM == 2-class LCA** to `tolerance = 1e-7`
  (`test-gllamm-consistency.R:14`) — ties the two independent EM cores together.
- **EIRT with item residuals == crossed random-effects GLLAMM** on VerbAgg to
  `1e-3` on logLik and both random-effect SDs
  (`test-gllamm-consistency.R:28-33`); a separate multilevel-EIRT-vs-`glmer`
  cross-walk sits at `0.02` (`test-eirt-audit.R:86`).
- **Saturated EIRT == descriptive IRT** (`difficulty_formula = ~ item`) to
  `1e-5` on logLik (`test-eirt-audit.R:27-28`).
- **Frequency weights == duplicated data** to `1e-6` on the Laplace/EIRT route
  and `1e-4` on the EM route (`test-eirt-audit.R:44-54`), with further
  level-1/level-2/aghq/LCA/CDM duplication identities in `test-level-weights.R`
  and `test-irt-em.R:87`. This is the empirical backstop for the weight-handling
  design in §2.

**Current status.** The cross-package validation table reports **93/93 checks
pass**, regenerated 2026-07-20 against gllammr **1.3.0** / R 4.5.1 / TMB 1.9.17
(`validation/RESULTS.md:3-5,11`). The testthat suite is 57 files; it is what
runs under `R CMD check` (the cross-package subset gated off CRAN — §5).

**Known limitations (own them before the reviewer finds them).** 2PL/3PL
discrimination on very short tests (≈5 items) can diverge under the Laplace
approximation — a limitation *shared with joint ML*, documented at
`validation/RESULTS.md:109-113` and in `R/irt.R`'s roxygen. The mitigations are
in the box already: the EM path (`method = "em"`) handles the 5-item LSAT 2PL
the way ltm's EM does (validated, `twopl_lsat_em`, `validation.R:645-667`), and
the 2PL Laplace path is validated instead on a 20-item simulated test where each
person carries enough information (`twopl_simulated`, correlations `> 0.995`,
`validation.R:310-345`). 3PL, as noted in §3, is the least-constrained and
least-validated model.

---

## 5. Anticipated reviewer questions

**Q: Why does this package exist next to lme4 / glmmTMB / mirt / TAM / lavaan /
poLCA?**
A: It is one unified GLLAMM (Rabe-Hesketh–Skrondal–Pickles) framework with one
call — `gllamm(formula, data, family = ...)` — and one weight/prediction/
diagnostic vocabulary spanning GLMMs, IRT, EIRT, DIF, LCA, CDM, SEM,
mixed-response and frailty-survival models. The value is (i) a migration path
for Stata GLLAMM users, and (ii) the models that live *at the intersections*,
which no single reference package covers: explanatory IRT with item covariates
and random item residuals (LLTM+error), order-restricted / cognitive-diagnostic
LCA, IRT-LR DIF with a latent-regression impact term, joint mixed-type
outcomes. The validation suite deliberately shows that on the *overlap* with
each specialist package the numbers agree, so the framework is a superset, not
a competitor with different answers.

**Q: Why is the installed size ~14 MB?**
A: TMB template compilation. Twenty-two model templates compile into one shared
object that shares a single `TMB.hpp`/Eigen/CppAD include; `libs` is ~13 MB.
This is typical of TMB packages (cf. glmmTMB) and is noted as an expected NOTE
in `cran-comments.md:30-32`. There is no run-time compilation and nothing is
cached in the user's home directory (`cran-comments.md:37-46`).

**Q: Why are cross-package tests skipped on CRAN?**
A: Two reasons, both in `cran-comments.md:46-59`. (1) *Upstream-drift / archival
risk*: those tests assert agreement with other packages' numerics, so an
upstream release that changes an estimate could fail the suite and threaten
archival for something outside our control. (2) *A real segfault*: lme4 2.0-1
with Matrix 1.7-5 can segfault (not error — unrecoverable) inside `glmer()` on
Windows. So every such test carries `skip_on_cran()`, centralised in
`tests/testthat/helper-reference.R`'s `ref_fit()` wrapper, which additionally
skips reference fits outright when `GITHUB_ACTIONS == "true"` and the platform
is Windows, and converts any other reference-fit error into a `skip()`. These
run in full on GitHub Actions where `NOT_CRAN = true`. Every family still has
fast unit tests and at least one end-to-end smoke fit that runs on CRAN, so the
compiled code is exercised on every platform.

**Q: Why 22 Suggests?**
A: They are the reference packages for the validation cross-walks and the
vignette examples — lme4, mirt, TAM, MASS, lavaan, npmlreg, glmmTMB, ordinal,
poLCA, survival, HSAUR3, ltm, CDM, difR, VGAM, nnet, clubSandwich, plus the
tooling (testthat, R6, knitr, rmarkdown, ggplot2). Every use is guarded by
`requireNamespace(..., quietly = TRUE)` (see the top of each `.validate_*`
function) or `skip_if_not_installed`, so none is required to install, check-lite,
or use the package.

**Q: Licence / copyright — is any third-party code bundled?**
A: GPL-3 (`DESCRIPTION`, `LICENSE`), sole author/maintainer Josh McGrane. The
`src/` tree contains only the package's own code: `gllammr.cpp`, the 22
`gllamm_*.hpp` templates, and the three hand-written cores `em_poly.cpp`,
`em_cdm.cpp`, `isotonic.cpp`. No third-party source is vendored — Eigen, TMB and
CppAD enter *only* through `LinkingTo: TMB, RcppEigen` headers (the standard TMB
mechanism), and the only non-package includes are `<TMB.hpp>`, `<R.h>`,
`<Rinternals.h>`, `<R_ext/BLAS.h>`, and standard-library headers. The
algorithms are original implementations that *mirror* published methods (poLCA
EM, Dykstra isotonic projection, SQUAREM S3), not copied code.

**Q: Test runtime on CRAN?**
A: With the cross-package and slow integration fits skipped, the CRAN subset is
designed to complete well inside the check-time budget; the full suite runs on
CI (`cran-comments.md:57-59`). Measured at release prep (2026-07-12, Apple
Silicon, under concurrent CPU load): the CRAN-simulated subset
(`NOT_CRAN=false`) ran 3,223 assertions in ~101 s versus ~4 min for the full
suite — roughly 2-3.5 min projected on CRAN hardware. If a reviewer wants an
authoritative number, measure fresh with `R CMD check` on the tarball.

**Q (procedural): what if the reviewer asks for changes?**
A: Apply the fix, bump the version to 1.3.1 in `DESCRIPTION`, add a NEWS entry,
regenerate any affected `man/*.Rd` with roxygen, update `cran-comments.md`
("Resubmission" + a short list of what changed and why), rebuild the tarball via
the `build-tarball` workflow (§6 — the local machine cannot build it), and
resubmit through the same flow.

---

## 6. Operational appendix

**Local dev quirks.** The maintainer's macOS machine lacks the gfortran
toolchain that `R CMD build` / `R CMD check` need on this R version, so local
verification is `testthat::test_local()` against the already-compiled package,
*not* a local `R CMD build`/`check` (`cran-comments.md:9-13`). Build the CRAN
tarball via the **`build-tarball`** GitHub Actions workflow
(`workflow_dispatch`-only; it runs `R CMD build .` with all deps and uploads the
source tarball as an artifact), never locally.

**CI workflow map** (`.github/workflows/`):
- `R-CMD-check.yaml` — `R CMD check --as-cran` on push/PR across macOS-release,
  windows-release, ubuntu-devel/release/oldrel-1 (`fail-fast: false`). Contains
  the Windows-only Eigen `-Warray-bounds` suppression, injected at runtime via
  `R_MAKEVARS_USER` — deliberately *not* shipped in `src/Makevars.win`
  (`cran-comments.md:62-70`) so it can't affect a CRAN build. `NOT_CRAN=true`
  here is what makes the `skip_on_cran()` cross-package tests run.
- `build-tarball.yaml` — manual; produces the submission tarball.
- `validation.yaml` — weekly cron + manual; runs `validation/run_validation.R`
  against the full reference-package set and fails if fewer than all checks pass
  or any row errored.
- `test-coverage.yaml` — covr → Codecov on push/PR.
- `lint.yaml` — `lintr::lint_package()`, report-only (never blocks).
- `pkgdown.yaml` — builds/deploys the docs site.

**`src/Makevars` and `src/Makevars.win`** are identical: OpenMP via the standard
`$(SHLIB_OPENMP_CXXFLAGS)` in both `PKG_CXXFLAGS` and `PKG_LIBS`,
`-DTMB_EIGEN_DISABLE_WARNINGS` in `PKG_CPPFLAGS`, and the standard
`LAPACK_LIBS`/`BLAS_LIBS`/`FLIBS`. No `-Wno-array-bounds` in either — that lives
only on the CI Windows runner. The package never forces a thread count.

**Where the lme4-Windows guard lives.** Two mirror copies of the same policy:
`tests/testthat/helper-reference.R` `ref_fit()` for the test suite, and
`.reference_fit()` / `.reference_skip()` in `R/validation.R:100-127` for the
validation harness. Both skip lme4 reference fits when
`GITHUB_ACTIONS == "true"` on Windows, and convert any other reference-fit error
into a skip. If lme4/Matrix upstream is fixed, this guard can be relaxed.

**Test / validation gate commands.** `testthat::test_local()` (or
`devtools::test()`) for the unit suite; `gllammr_validate()` for the
cross-package numerical cross-walks (writes the table behind
`validation/RESULTS.md`); `gllammr_validate(scale = "large")` for the
large-n tier.
