# Mixed-Response Family for Joint Outcome Models

Create a family object for joint models of mixed-type outcomes sharing a
random effect, through the unified
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
interface. The first argument of
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
is the shared random-effects formula (e.g. `~ 1 | group`). See
[`fit_mixed`](https://drjoshmcgrane.github.io/gllammr/reference/fit_mixed.md).

## Usage

``` r
mixed_response(...)
```

## Arguments

- ...:

  Named outcome formulas: `gaussian = y1 ~ x`, `binomial = y2 ~ x`,
  `poisson = y3 ~ x` (any subset)

## Value

An object of class `mixed_family`

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- gllamm(~ 1 | clinic, data = d,
              family = mixed_response(gaussian = severity ~ age,
                                      binomial = dropout ~ age))
} # }
```
