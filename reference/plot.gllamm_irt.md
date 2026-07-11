# Plot IRT Model Diagnostics

Create diagnostic plots for IRT models including item characteristic
curves, item information functions, test information, and ability
distributions

## Usage

``` r
# S3 method for class 'gllamm_irt'
plot(x, which = 1:4, items = NULL, ...)
```

## Arguments

- x:

  A gllamm_irt object

- which:

  Which plots to produce (1=ICC, 2=IIF, 3=TIF, 4=Ability)

- items:

  Which items to plot (default: first 6 items)

- ...:

  Additional arguments passed to plotting functions

## Details

Plot types:

- `which = 1`: Item Characteristic Curves (ICC)

- `which = 2`: Item Information Functions (IIF)

- `which = 3`: Test Information Function (TIF)

- `which = 4`: Person Ability Distribution

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit 2PL model
fit <- fit_irt(responses, model = "2PL")

# Plot all diagnostics for items 1-3
plot(fit, which = 1:4, items = 1:3)

# Plot only ICCs
plot(fit, which = 1, items = 1:5)
} # }
```
