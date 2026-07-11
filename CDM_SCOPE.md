# Scope: Cognitive Diagnosis Models via item-specific partial orders

Status: phases 1-2 IMPLEMENTED 2026-06-11 (R/cdm.R, tests/testthat/
test-cdm.R, validation case cdm_fraction_dina). Phase 3 items remain
deferred. Implementation notes vs scope: dedicated
[`fit_cdm_em()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_cdm_em.md)
(reusing
[`.isotonic_poset()`](https://drjoshmcgrane.github.io/gllammr/reference/dot-isotonic_poset.md)/[`.topological_order()`](https://drjoshmcgrane.github.io/gllammr/reference/dot-topological_order.md))
instead of generalizing `fit_lca_em` - simpler and zero risk to LCA
paths; complete-data fast path (one BLAS product per E-step, denominator
= class counts); pi floored at 1e-12 (2^A profiles legitimately hit zero
prevalence); defaults n_starts = 3, tol = 1e-7. PHASE-3 SPEED ITEM DONE:
full EM loop in C++ (src/em_cdm.cpp) on dgemm kernels with SQUAREM
acceleration (the iteration count, not the arithmetic, was the wall:
~490 EM steps -\> ~115 E-steps). Fraction-subtraction: 29.7s (initial R)
-\> 1.4s for 3 starts (~0.45s per start; CDM::din 0.14s at much looser
convergence and 1 start, and our logLik is 0.07 better). Large G-DINA
(20000 x 30, A = 5): 2.7s for 3 starts vs CDM::gdina 1.9s for 1.
Prerequisite work: ordered LCA (commit 32716f3), partially ordered LCA
(commit 149aa17) — the isotonic-poset M-step machinery this builds on.

## The bridge from where we are

The partially ordered LCA constrains ALL items by ONE class poset. A CDM
is the natural endpoint of relaxing that in two steps:

1.  **Item-specific posets**: each item gets its own set of class-pair
    constraints. Mechanically trivial —
    [`.isotonic_poset()`](https://drjoshmcgrane.github.io/gllammr/reference/dot-isotonic_poset.md)
    already runs per item row; only the API and bookkeeping change.
2.  **Structure from a Q-matrix**: classes become attribute profiles (A
    binary attributes -\> K = 2^A classes), and each item’s constraints
    are *derived* rather than user-listed: item j depends only on the
    attributes it measures (equality across classes sharing the item’s
    reduced profile) and is monotone in those attributes (isotonic over
    the reduced sub-lattice). That is exactly the monotone saturated
    G-DINA / LCDM, with DINA and DINO as collapsed special cases.

Step 2 subsumes step 1 for every use case we care about, so the
deliverable is a CDM interface, not a raw item-specific-poset API.

## Statistical core (all closed-form EM, no new estimation theory)

- E-step: identical to the existing LCA E-step (N x K posterior weights,
  BLAS matrix products), K = 2^A (A \<= ~8 practical).
- M-step per item j with k_j measured attributes:
  1.  collapse the K posterior columns into 2^{k_j} groups by the item’s
      reduced profile (a weighted column-sum — closed form);
  2.  weighted success proportions per group;
  3.  isotonic projection over the reduced attribute lattice using the
      existing
      [`.isotonic_poset()`](https://drjoshmcgrane.github.io/gllammr/reference/dot-isotonic_poset.md)
      (subset order on {0,1}^{k_j});
  4.  expand back to K columns. DINA/DINO: same with exactly 2 groups
      (all-required-mastered vs not / any-measured-mastered vs none),
      giving guess g_j and slip s_j with monotonicity g_j \<= 1 - s_j
      enforced by the same projection.
- Profile prevalences pi: saturated multinomial M-step (as now).
- Ramsay acceleration + isotonic projection of extrapolations carries
  over unchanged.
- Attribute hierarchies (linear, convergent, divergent) = pruning
  disallowed profiles from the class space before fitting — cheap and
  high-value; include in phase 1.

## API sketch

``` R
fit_cdm(Y, Q,
        model = c("gdina", "dina", "dino"),   # gdina = monotone saturated
        hierarchy = NULL,    # optional list of attribute pairs a -> b
                             # (mastery of b requires a); prunes profiles
        monotone = TRUE,     # FALSE = unconstrained saturated G-DINA
        weights = NULL, control = list())
```

Returns `gllamm_cdm` (inheriting `gllamm`): - `item_params`: per item,
P(Y=1 \| reduced profile) table (for DINA/DINO also g/s), with the
reduced-profile labels - `profile_probs` (named by profile bit-pattern),
`posterior` (N x K), - `attribute_posteriors` (N x A marginal mastery
probabilities), `modal_profile`, `modal_attributes` - logLik / AIC / BIC
(nominal df; chi-bar-square caveat documented), EM convergence info -
print/summary: item tables, attribute prevalence, classification
certainty (mean max posterior)

Internals: generalize
[`fit_lca_em()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_lca_em.md)
with an optional per-item grouping matrix + per-item edge list rather
than forking a new EM.

## Validation plan (comparator: CDM package — already installed)

| Check | Data | Reference | Tolerance |
|----|----|----|----|
| DINA g/s + logLik | fraction.subtraction (de la Torre 2009 classic) | CDM::din() | logLik \< 0.5 abs; g/s \< 1e-2 |
| DINO | simulated | CDM::din(rule=“DINO”) | same |
| Saturated G-DINA | simulated (A=3, J=15) | CDM::gdina() | item params \< 1e-2 |
| Monotone G-DINA | simulated with non-monotone noise | internal: logLik \<= unconstrained; recovery | recovery |
| Attribute classification | fraction.subtraction | CDM::din() MAP profiles | agreement \> 95% |
| Hierarchy pruning | simulated linear hierarchy | profile space + recovery | recovery |

Add 1-2 cases to
[`gllammr_validate()`](https://drjoshmcgrane.github.io/gllammr/reference/gllammr_validate.md)
(CDM goes in Suggests) + a testthat file mirroring the ordered-LCA one
(M-step equivalences: DINA == 2-group collapse; G-DINA with q_j =
all-ones == saturated LCA column; monotone == unconstrained when
constraints don’t bind).

## Phasing and effort

- **Phase 1 — core** (one working session): fit_cdm with
  dina/dino/gdina + monotone + hierarchy; print/summary; testthat suite.
  Risk: low — every piece reuses proven machinery.
- **Phase 2 — validation** (half session): CDM-package cross-checks,
  harness cases, NEWS/docs/vignette section.
- **Phase 3 — deferred** (explicitly out of scope now):
  - polytomous responses (sequential G-DINA)
  - person covariates on attribute membership (structural loglinear /
    logistic model instead of saturated pi)
  - M2 / RMSEA2 limited-information fit statistics
  - Q-matrix validation/discovery (GDINA::Qval style)
  - C++ E-step for the large tier (only if A \>= 8 with N \>= 50k shows
    up in benchmarks; the LCA EM experience says R/BLAS will hold to
    there)

## Decisions taken in this scope (flag if you disagree)

1.  Separate
    [`fit_cdm()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_cdm.md)
    rather than overloading
    [`fit_lca()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_lca.md)
    — profile semantics, Q-matrix, and outputs differ too much; shared
    internals.
2.  Validate against CDM (installed) not GDINA (not installed); GDINA
    can be added later for a second opinion.
3.  Monotonicity ON by default (`monotone = TRUE`) — interpretability
    first, matching the package’s ordering-restriction direction.
