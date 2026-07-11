# Multinomial Family for Unordered Categorical Outcomes

Create a family object for baseline-category multinomial logit models
through the unified
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
interface.

## Usage

``` r
multinomial(reference = NULL)
```

## Arguments

- reference:

  Reference category (default: first level)

## Value

A family object of class `multinomial_family`

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- gllamm(choice ~ x + (1 | region), data = d, family = multinomial())
} # }
```
