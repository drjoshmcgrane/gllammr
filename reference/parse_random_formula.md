# Parse Random Effects Formula

Internal function to parse lme4-style random effects formulas

## Usage

``` r
parse_random_formula(formula, data)
```

## Arguments

- formula:

  Random effects formula (e.g., ~ (1 \| class) + (1 \| school))

- data:

  Data frame containing grouping variables

## Value

List of parsed random effects terms
