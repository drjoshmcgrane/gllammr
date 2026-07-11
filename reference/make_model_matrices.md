# Extract model matrices

Create design matrices for fixed and random effects

## Usage

``` r
make_model_matrices(parsed_formula, data)
```

## Arguments

- parsed_formula:

  Parsed formula object from parse_formula()

- data:

  Data frame

## Value

List with X (fixed effects), Z (random effects), and grouping info
