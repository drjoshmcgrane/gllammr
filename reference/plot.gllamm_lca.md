# Plot Latent Class Analysis Results

Create diagnostic plots for latent class models including class
profiles, item probability heatmaps, and classification summaries

## Usage

``` r
# S3 method for class 'gllamm_lca'
plot(x, which = 1:3, ...)
```

## Arguments

- x:

  A gllamm_lca object

- which:

  Which plots to produce (1=profiles, 2=heatmap, 3=classification)

- ...:

  Additional arguments passed to plotting functions

## Details

Plot types:

- `which = 1`: Class Profiles (line plot of item probabilities by class)

- `which = 2`: Item Probability Heatmap

- `which = 3`: Classification Summary (barplot of class assignments)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit LCA model
fit <- fit_lca(indicators, nclass = 3)

# Plot all diagnostics
plot(fit, which = 1:3)

# Plot only class profiles
plot(fit, which = 1)
} # }
```
