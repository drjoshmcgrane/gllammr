# Extract Person Abilities from IRT Models

Extract estimated person abilities (theta) from IRT models. For
multi-level models, this returns the person-level deviations (theta_0).
Use `composite_theta` to get total abilities including random effects.

## Usage

``` r
abilities(object, ...)
```

## Arguments

- object:

  A fitted IRT model (theta_0 + random effects) instead of just theta_0.
  Default FALSE.

- ...:

  Additional arguments (not used)

## Value

A named vector of person abilities

## Examples

``` r
if (FALSE) { # \dontrun{
# Person-level deviations
abilities(fit)

# Total abilities (including class effects)
abilities(fit, composite = TRUE)
} # }
```
