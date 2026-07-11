# Fit Explanatory Item Response Theory Models

Fit IRT models where item parameters are modeled as functions of item
covariates. For polytomous models, supports GRM (cumulative logit), PCM
(adjacent-categories), and GPCM (adjacent-categories with
discrimination). PCM can include optional threshold-level predictors via
threshold_formula.

## Usage

``` r
fit_eirt(
  response_matrix,
  item_data,
  difficulty_formula = ~1,
  discrimination_formula = ~1,
  threshold_formula = NULL,
  step_formula = NULL,
  step_data = NULL,
  person_data = NULL,
  random = NULL,
  weights = NULL,
  model = c("Rasch", "2PL", "GRM", "PCM", "GPCM"),
  item_residuals = TRUE,
  start = NULL,
  control = list()
)
```

## Arguments

- response_matrix:

  Matrix of item responses (persons x items)

- item_data:

  Data frame of item-level covariates (must have n_items rows)

- difficulty_formula:

  Formula for item location/difficulty regression (e.g., ~ word_freq).
  Model: b_i = W_diff %\*% gamma + epsilon_b (if item_residuals = TRUE)

- discrimination_formula:

  Formula for discrimination regression (e.g., ~ item_type). Applies to
  2PL, GRM, and GPCM models. Model: log(a_i) = W_disc %\*% delta +
  epsilon_a (if item_residuals = TRUE)

- threshold_formula:

  Formula for threshold-specific regression (e.g., ~ abstractness).
  Coefficients are step-specific DEVIATIONS that sum to zero across an
  item's thresholds; the item-level main effect of the same covariate
  belongs in difficulty_formula (the two are jointly identified this
  way). Only used for PCM/GPCM models. When specified, enables
  threshold-difficulty regression (Kim & Wilson 2019 LPCM framework):
  delta_im = b_i + sum_k xi_km \* x_ik + e_im

- step_formula:

  Formula for step-level covariates (PCM/GPCM only): predictors that
  vary WITHIN an item across its steps, each with a single common
  coefficient (\\\delta\_{im} = b_i + \xi\_{0m} + \sum_k \eta_k
  x\_{imk} + e\_{im}\\). The intercept is dropped; per-step baselines
  come from the threshold intercepts. Combine freely with
  difficulty_formula (item level) and threshold_formula (item properties
  with step-specific effects).

- step_data:

  Data frame for step_formula with one row per item-step cell in
  item-major order (item 1 steps 1..K-1, then item 2, ...; n_items x
  (max_categories - 1) rows; pad items with fewer categories with NA
  rows - those cells are never read).

- person_data:

  Optional data frame with person-level variables for multi-level models

- random:

  Optional random effects formula (e.g., ~ (1 \| class))

- weights:

  Optional vector of integer person-level frequency weights (length =
  number of persons). Implemented by exact replication of weighted
  persons, so results are identical to fitting the duplicated data.
  Non-integer weights are not supported under the Laplace approximation.
  Length must equal number of persons. Default: all observations
  weighted equally. Weights are expanded to observation-level
  (person-item) internally.

- model:

  IRT model type: "Rasch" or "2PL" (dichotomous), or "GRM", "PCM",
  "GPCM" (polytomous)

- item_residuals:

  Logical. If TRUE (default), includes item-specific residuals (LLTM +
  error). If FALSE, uses pure LLTM where item parameters are exactly
  predicted by covariates with no residuals.

- start:

  Optional starting values list

- control:

  Control parameters for nlminb optimization

## Value

An object of class `gllamm_eirt`

## Details

\*\*Dichotomous models:\*\* - \*\*Rasch\*\*: P(Y=1) = logit^(-1)(theta -
b_i), where b_i = W_diff %\*% gamma \[+ epsilon_b\] - \*\*2PL\*\*:
P(Y=1) = logit^(-1)(a_i \* (theta - b_i)), where log(a_i) = W_disc %\*%
delta \[+ epsilon_a\]

\*\*Polytomous models (Kim & Wilson 2019; De Boeck & Wilson 2004
framework):\*\* - \*\*GRM\*\*: Cumulative logit with ordered thresholds
expressed as sum-to-zero deviations around the item location b_i (so the
difficulty regression is identified). Supports discrimination_formula. -
\*\*PCM\*\*: Adjacent-categories logit, two-fold parameterization (MFRM
approach): delta_im = b_i + s_im, where sum(s_im) = 0 across steps for
each item. When threshold_formula is specified, uses threshold
regression: delta_im = b_i + sum_k xi_km \* x_ik + e_im (Kim & Wilson
LPCM framework) - \*\*GPCM\*\*: As PCM with item-specific discrimination
a_i. Supports both discrimination_formula and threshold_formula.

\*\*Relation to the GLLAMM framework:\*\* difficulties enter the linear
predictor linearly, exactly as in the canonical GLLAMM formulation.
Discriminations are GLLAMM factor loadings; the discrimination
regression models them on the log scale (a_i = exp(W %\*% delta \[+
epsilon_a\])), a reparameterization that keeps loadings positive and
makes covariate effects multiplicative. A model with
discrimination_formula = ~ 1 and item_residuals = FALSE has constant
loadings and is GLLAMM-canonical; covariate-structured or random
loadings are a (standard, benign) nonlinear extension.

\*\*Item residuals:\*\* When item_residuals = TRUE (default), item
parameters include residual terms epsilon_b and epsilon_a (LLTM +
error). When FALSE, uses pure LLTM where parameters are exactly
determined by covariates.

## Examples

``` r
if (FALSE) { # \dontrun{
# Dichotomous EIRT with item and discrimination predictors
item_chars <- data.frame(
  word_frequency = rnorm(20),
  item_length = rpois(20, 5),
  item_type = factor(sample(c("concrete", "abstract"), 20, replace = TRUE))
)
fit_2pl <- fit_eirt(responses, item_data = item_chars,
                    difficulty_formula = ~ word_frequency + item_length,
                    discrimination_formula = ~ item_type,
                    model = "2PL")

# Pure LLTM (no residuals)
fit_lltm <- fit_eirt(responses, item_data = item_chars,
                     difficulty_formula = ~ word_frequency,
                     model = "Rasch",
                     item_residuals = FALSE)

# Polytomous PCM (adjacent-categories logit, Rasch family)
fit_pcm <- fit_eirt(poly_responses, item_data = item_chars,
                    difficulty_formula = ~ abstractness,
                    model = "PCM")

# PCM with threshold predictors (LPCM framework)
fit_pcm_thresh <- fit_eirt(poly_responses, item_data = item_chars,
                           difficulty_formula = ~ abstractness,
                           threshold_formula = ~ cognitive_level,
                           model = "PCM")

# GPCM with all predictors
fit_gpcm <- fit_eirt(poly_responses, item_data = item_chars,
                     difficulty_formula = ~ word_frequency,
                     discrimination_formula = ~ item_type,
                     threshold_formula = ~ cognitive_level,
                     model = "GPCM")
} # }
```
