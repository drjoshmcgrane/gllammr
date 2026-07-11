# Confirmatory model-based DIF (IRT likelihood-ratio tests)

Tests differential item functioning inside the joint marginal-ML IRT
model (IRT-LR DIF; Thissen, Steinberg & Wainer 1993), the De Boeck &
Wilson item-by-covariate-interaction formulation. Group differences in
ability (impact) are modeled by a latent regression \\\theta_p \sim
N(z_p'\gamma, \sigma^2)\\, so DIF parameters measure item bias beyond
true ability differences - the model-based companion to the
observed-criterion screening tests in
[`dif_test`](https://drjoshmcgrane.github.io/gllammr/reference/dif_test.md).

## Usage

``` r
dif_irt(
  response_matrix,
  dif,
  person_data = NULL,
  model = c("Rasch", "2PL"),
  items = NULL,
  anchors = NULL,
  type = c("uniform", "both"),
  method = c("lr", "wald"),
  purify = FALSE,
  alpha = 0.05,
  p_adjust = "none",
  max_iter = 10,
  control = list()
)
```

## Arguments

- response_matrix:

  Binary item response matrix (persons x items; NA allowed)

- dif:

  DIF specification: a grouping vector or a one-sided formula over
  `person_data` (multiple variables and interactions supported, as in
  [`dif_test`](https://drjoshmcgrane.github.io/gllammr/reference/dif_test.md))

- person_data:

  Data frame with the DIF variables (required when `dif` is a formula)

- model:

  "Rasch" (default) or "2PL"

- items:

  Item indices to test (default: all non-anchor items)

- anchors:

  Optional indices of items constrained DIF-free

- type:

  "uniform" (default; covariate shifts of the item logit,
  \\z'\delta_i\\) or "both" (additionally nonuniform DIF: covariate
  scaling of the discrimination, \\a_i e^{z'\kappa_i}\\; 2PL only)

- method:

  "lr" (default; per-item likelihood-ratio tests, refitting the joint
  model) or "wald" (one joint fit with all studied items free - requires
  explicit `anchors` for identification - and per-item block Wald tests)

- purify:

  Purified IRT-LR (default FALSE): after each round, items flagged so
  far keep free DIF parameters in both compared models, so their misfit
  cannot contaminate the impact estimate; iterate until the flag set
  stabilizes

- alpha:

  Significance level

- p_adjust:

  Multiple-testing correction (`p.adjust` method)

- max_iter:

  Maximum purification rounds

- control:

  Optimization control list

## Value

An object of class `dif_irt`: `dif_results` (per item: LR or Wald
chi-square, df, p, adjusted p, flag, and the estimated uniform DIF
effects \\\delta\\ on the logit metric with standard errors), `impact`
(latent regression coefficients \\\gamma\\ with SEs - the ability
difference attributable to the covariates), `flagged_items`, and
purification details.

## Details

For a studied item the compared models differ only in that item's DIF
parameters; all item difficulties (and discriminations) are estimated
jointly with the latent regression, so the test is a genuine
marginal-likelihood ratio with q (uniform) or 2q ("both") degrees of
freedom, q the number of DIF design columns. For the Rasch model with
uniform DIF this is likelihood-equivalent to the long-format GLMM
`y ~ 0 + item + z + item_j:z + (1 | person)` under the same Laplace
approximation (verified in the test suite).

## Examples

``` r
if (FALSE) { # \dontrun{
# Screen, then confirm
screen <- dif_test(resp, dif = ~ gender * language, person_data = pd)
confirm <- dif_irt(resp, dif = ~ gender * language, person_data = pd,
                   items = screen$flagged_items)
summary(confirm)
} # }
```
