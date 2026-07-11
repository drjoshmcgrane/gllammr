# Fit Cognitive Diagnosis Models

Fits cognitive diagnosis models (CDMs) for binary responses: latent
classes are attribute profiles (combinations of A binary skills), and a
Q-matrix declares which attributes each item measures. Item response
probabilities depend only on the item's measured attributes (its
"reduced profile"), giving the saturated G-DINA / LCDM family, with DINA
and DINO as two-group special cases.

## Usage

``` r
fit_cdm(
  Y,
  Q,
  model = c("gdina", "dina", "dino"),
  hierarchy = NULL,
  monotone = TRUE,
  weights = NULL,
  control = list()
)
```

## Arguments

- Y:

  Binary response matrix (persons x items; NA allowed)

- Q:

  Binary Q-matrix (items x attributes): `Q[j, a] = 1` if item j measures
  attribute a. Every row must have at least one 1.

- model:

  Item model:

  - `"gdina"` (default): saturated G-DINA - one free response
    probability per distinct reduced profile of the item's measured
    attributes.

  - `"dina"`: conjunctive - two probabilities per item (guessing for
    profiles missing any required attribute, 1 - slip for profiles
    mastering all of them).

  - `"dino"`: disjunctive - mastery of any measured attribute suffices.

- hierarchy:

  Optional attribute hierarchy: a list of attribute index pairs (or
  two-column matrix), each pair `c(a, b)` meaning attribute `a` is a
  prerequisite of attribute `b`; profiles with `b` mastered but not `a`
  are removed from the latent space. Attribute names (colnames of Q) may
  be used instead of indices.

- monotone:

  Constrain response probabilities to be nondecreasing in the attributes
  (default TRUE): mastering more of the measured attributes can never
  lower the success probability. Implemented as weighted isotonic
  regression over the reduced-profile lattice in the M-step, so
  estimation remains closed-form EM.

- weights:

  Optional vector of case weights (one per person)

- control:

  List: `n_starts` (default 3), `max_iter` (default 2000), `tol`
  (default 1e-7, absolute logLik change)

## Value

An object of class `gllamm_cdm` with components including `item_params`
(per item, P(Y = 1) by reduced profile; for DINA/DINO also
`guess`/`slip`), `profile_probs`, `posterior` (persons x profiles),
`attribute_posteriors` (persons x attributes marginal mastery
probabilities), `modal_profile`, and `logLik`/`AIC`/`BIC`.

## Details

Estimation is marginal maximum likelihood via EM with closed-form
M-steps (weighted proportions pooled over the item's reduced-profile
groups, isotonically projected when `monotone = TRUE`) and safeguarded
Ramsay acceleration. AIC/BIC use the nominal parameter count;
likelihood-ratio tests against less constrained models have non-standard
null distributions when monotonicity binds.

## Examples

``` r
if (FALSE) { # \dontrun{
Q <- rbind(c(1, 0), c(0, 1), c(1, 1), c(1, 0), c(0, 1))
fit <- fit_cdm(Y, Q, model = "dina")
summary(fit)
fit$attribute_posteriors  # P(mastery) per person and attribute
} # }
```
