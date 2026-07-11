# Plot Ordinal Model Effects for Multiple Covariates

Compare effects across multiple covariates in an ordinal model

## Usage

``` r
plot_ordinal_effects(
  object,
  covariates = NULL,
  sort_by = c("magnitude", "name", "none"),
  ...
)
```

## Arguments

- object:

  A gllamm_ordinal object

- covariates:

  Vector of covariate names (default: all non-intercept)

- sort_by:

  Sort covariates by: "magnitude", "name", or "none"

- ...:

  Additional arguments

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_ordinal(rating ~ temp + contact + (1 | judge), data = wine)
plot_ordinal_effects(fit)
} # }
```
