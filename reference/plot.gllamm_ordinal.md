# Plot Ordinal Regression Model Diagnostics

Create diagnostic plots for ordinal models including cumulative
probabilities, category probabilities, threshold parameters, and
covariate effects

## Usage

``` r
# S3 method for class 'gllamm_ordinal'
plot(x, which = 1:3, covariate = NULL, covariate_values = NULL, ...)
```

## Arguments

- x:

  A gllamm_ordinal object

- which:

  Which plots to produce (1=cumulative, 2=category, 3=thresholds,
  4=effects)

- covariate:

  Name of covariate to plot (default: first non-intercept covariate)

- covariate_values:

  Optional vector of covariate values to plot (default: -2 to 2)

- ...:

  Additional arguments passed to plotting functions

## Details

Plot types:

- `which = 1`: Cumulative Probabilities P(Y \<= k) vs covariate

- `which = 2`: Category Probabilities P(Y = k) vs covariate

- `which = 3`: Threshold Parameters on latent scale

- `which = 4`: Covariate Effects (shows non-proportional effects for
  PPO)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit ordinal model
fit <- fit_ordinal(rating ~ temp + contact + (1 | judge),
                   data = wine, link = "logit")

# Plot all diagnostics for 'temp' covariate
plot(fit, which = 1:4, covariate = "temp")

# Plot only cumulative probabilities
plot(fit, which = 1, covariate = "contact")
} # }
```
