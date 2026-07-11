# Test for Differential Item Functioning (DIF)

Logistic-regression DIF (Swaminathan & Rogers 1990; Zumbo 1999) with
latent-trait or observed-score matching, generalized to multiple DIF
variables and their interactions through a formula interface, with
iterative purification of the matching criterion. Polytomous items are
tested with cumulative-logit (proportional odds) regression.

## Usage

``` r
dif_test(
  response_matrix,
  dif,
  person_data = NULL,
  model = c("auto", "Rasch", "2PL", "GRM", "PCM", "GPCM"),
  match = c("theta", "score"),
  items = NULL,
  type = c("both", "uniform", "nonuniform"),
  purify = TRUE,
  anchors = NULL,
  alpha = 0.05,
  p_adjust = "none",
  max_iter = 10
)
```

## Arguments

- response_matrix:

  Item response matrix (persons x items): binary 0/1 or polytomous 1..K
  (NA allowed)

- dif:

  The DIF specification: either a single grouping vector, or a one-sided
  formula over columns of `person_data` - e.g. `~ gender`,
  `~ gender + language`, or `~ gender * language` (the interaction tests
  whether DIF for one factor differs by the level of the other).

- person_data:

  Data frame with the DIF variables (required when `dif` is a formula)

- model:

  Matching (measurement) model for the latent criterion: "auto"
  (default: 2PL for dichotomous, GRM for polytomous), or any of "Rasch",
  "2PL", "GRM", "PCM", "GPCM"

- match:

  "theta" (default; EAP score from the anchor items under `model`) or
  "score" (observed anchor-set total score, the classical
  Swaminathan-Rogers criterion, comparable to
  [`difR::difLogistic`](https://rdrr.io/pkg/difR/man/difLogistic.html))

- items:

  Item indices to test (default: all)

- type:

  "both" (default; joint test of uniform + nonuniform DIF), "uniform"
  (group effects given the criterion), or "nonuniform" (criterion x
  group interactions given the uniform terms)

- purify:

  Iteratively purify the matching criterion (default TRUE): re-derive it
  from the currently unflagged items and re-test, until the flagged set
  stabilizes

- anchors:

  Optional item indices guaranteed DIF-free; they are always part of the
  matching criterion and never tested

- alpha:

  Significance level for flagging (default 0.05)

- p_adjust:

  Multiple-testing correction passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) ("none" default;
  e.g. "BH", "holm")

- max_iter:

  Maximum purification iterations (default 10)

## Value

An object of class `dif_analysis`: `dif_results` (per item: LR
chi-square, df, p, adjusted p, Nagelkerke \\\Delta R^2\\ effect size
with the Jodoin-Gierl A/B/C classification, flag), `flagged_items`,
`anchor_items`, `purification` (iterations, history, converged), the
matching scores, and per-item full-model coefficients for plotting.

## Details

For each studied item the nested models \$\$M_0: y ~ m,\quad M_1: y ~
m + Z,\quad M_2: y ~ m + Z + m:Z\$\$ are fitted, where m is the matching
criterion and Z the design matrix of the DIF formula. Uniform DIF is
\\M_1\\ vs \\M_0\\, nonuniform \\M_2\\ vs \\M_1\\, and "both" \\M_2\\ vs
\\M_0\\. With multiple DIF variables each test has as many degrees of
freedom as Z has columns, and an interaction in the formula (e.g.
`~ g1 * g2`) tests intersectional DIF beyond the additive effects.
Effect sizes are Nagelkerke \\\Delta R^2\\ between the compared models
(A \< 0.035, B \< 0.07, C otherwise; Jodoin & Gierl 2001).

Purification (Lord 1980; Candell & Drasgow 1988): items flagged in one
round are removed from the matching criterion for the next, so DIF items
do not contaminate the score against which DIF is judged.

## Examples

``` r
if (FALSE) { # \dontrun{
# Single factor, purified
res <- dif_test(resp, dif = gender)

# Two factors plus their interaction, latent matching
res <- dif_test(resp, dif = ~ gender * language, person_data = persons)
summary(res)
dif_plot(res, item = 3, by = "gender")
} # }
```
