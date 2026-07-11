# Plot Item Covariate Effects

Visualize the relationship between item covariates and item parameters

## Usage

``` r
plot_item_covariates(
  object,
  covariate,
  parameter = c("difficulty", "discrimination"),
  ...
)
```

## Arguments

- object:

  A gllamm_eirt object

- covariate:

  Name of covariate to plot

- parameter:

  Which parameter to plot: "difficulty" or "discrimination"

- ...:

  Additional arguments passed to plot

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_eirt(responses, item_data,
                difficulty_formula = ~ word_freq)

plot_item_covariates(fit, covariate = "word_freq")
} # }
```
