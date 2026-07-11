# Fit Multinomial Regression Models with Random Effects

Fit baseline category logit models for nominal (unordered) responses

## Usage

``` r
fit_multinomial(
  formula,
  data,
  reference = NULL,
  start = NULL,
  control = list()
)
```

## Arguments

- formula:

  Formula with syntax: y ~ x + (terms \| group)

- data:

  Data frame

- reference:

  Reference category (default: first level)

- start:

  Optional starting values

- control:

  Control parameters

## Value

An object of class `gllamm_multinomial`

## Details

For nominal response Y with K categories, using baseline category logit:

\$\$P(Y = k \| x) = \frac{\exp(x'\beta_k)}{1 + \sum\_{j=1}^{K-1}
\exp(x'\beta_j)}\$\$

where category 0 is the reference with \\\beta_0 = 0\\.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate multinomial data
data$choice <- factor(sample(c("A", "B", "C"), 100, replace = TRUE))

# Fit multinomial model
fit <- fit_multinomial(choice ~ price + quality + (1 | person),
                       data = data)
summary(fit)
} # }
```
