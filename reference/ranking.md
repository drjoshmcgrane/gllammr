# Rank-Ordered Logit Family

Create a family object for rank-ordered (exploded) logit models through
the unified
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
interface. See
[`fit_rank`](https://drjoshmcgrane.github.io/gllammr/reference/fit_rank.md).

## Usage

``` r
ranking(case)
```

## Arguments

- case:

  Case (chooser) identifier: a one-sided formula (`~ chooser`) or a
  variable name

## Value

An object of class `rank_family`

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- gllamm(rank ~ price + quality, data = d,
              family = ranking(case = ~ chooser),
              random = ~ (1 | chooser))
} # }
```
