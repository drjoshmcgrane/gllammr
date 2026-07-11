# Test for DIF with explicit response data (deprecated)

Deprecated single-factor wrapper kept for backward compatibility; use
[`dif_test`](https://drjoshmcgrane.github.io/gllammr/reference/dif_test.md),
which supports multiple DIF variables, interactions, and iterative
purification.

## Usage

``` r
dif_test_with_data(
  response_matrix,
  group,
  model = c("Rasch", "2PL", "3PL", "GRM", "PCM", "GPCM"),
  items = NULL,
  type = c("both", "uniform", "nonuniform"),
  method = NULL,
  alpha = 0.05
)
```

## Arguments

- response_matrix:

  Matrix of item responses (persons x items)

- group:

  Grouping vector

- model:

  Matching model

- items:

  Items to test (default: all)

- type:

  Type of DIF

- method:

  Ignored (kept for compatibility)

- alpha:

  Significance level

## Value

Object of class `dif_analysis`
