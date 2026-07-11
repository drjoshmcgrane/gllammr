# Fit Rank-Ordered Logit Models with Taste Heterogeneity

Fits an exploded (rank-ordered) logit model for ranking data, with an
optional group-level random coefficient on one alternative attribute (a
"taste shifter"). A constant-within-case random intercept cancels in
ranking likelihoods, so heterogeneity must attach to an attribute that
varies across alternatives.

## Usage

``` r
fit_rank(
  formula,
  case,
  data,
  random = NULL,
  weights = NULL,
  start = NULL,
  control = list()
)
```

## Arguments

- formula:

  Model formula `rank ~ x1 + x2` where the response is the rank of each
  alternative within its case (1 = most preferred). Unranked
  alternatives may be coded `NA`; they remain in every choice set
  without contributing a stage.

- case:

  One-sided formula naming the case identifier, e.g. `~ id`

- data:

  Data frame in long format (one row per alternative per case)

- random:

  Optional one-sided formula `~ (0 + attribute | group)` giving the
  attribute carrying the random coefficient and the grouping variable.
  NULL fits a fixed-effects rank-ordered logit (the random SD is fixed
  near zero).

- weights:

  Optional one weight per case (matched via the case id)

- start:

  Optional starting values

- control:

  Optimization control list

## Value

An object of class `gllamm_rank`

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_rank(rank ~ price + quality, case = ~ subject,
                random = ~ (0 + price | region), data = d)
} # }
```
