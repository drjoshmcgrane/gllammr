# Plot Individual Classification Uncertainty

Visualize posterior probabilities for individual cases

## Usage

``` r
plot_classification_uncertainty(
  x,
  cases = 1:min(20, nrow(x$posterior)),
  sort_by = c("entropy", "modal", "index"),
  ...
)
```

## Arguments

- x:

  A gllamm_lca object

- cases:

  Which cases to plot (default: first 20)

- sort_by:

  Sort cases by: "entropy" (default), "modal", or "index"

- ...:

  Additional arguments

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_lca(data, nclass = 3)
plot_classification_uncertainty(fit, cases = 1:30)
} # }
```
