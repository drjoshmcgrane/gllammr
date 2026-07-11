# Parse GLLAMM formula

Parses a formula into fixed effects and random effects components.
Supports lme4-style syntax without requiring lme4.

## Usage

``` r
parse_formula(formula, data)
```

## Arguments

- formula:

  A formula object with syntax: y ~ x + (terms \| group)

- data:

  A data frame containing the variables

## Value

A list with components:

- fixed_formula:

  Formula for fixed effects

- random_terms:

  List of random effects specifications

- response_name:

  Name of response variable
