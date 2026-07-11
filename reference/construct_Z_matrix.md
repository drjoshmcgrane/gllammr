# Construct random effects design matrix for newdata

Construct random effects design matrix for newdata

## Usage

``` r
construct_Z_matrix(newdata, random_terms, group_var = NULL)
```

## Arguments

- newdata:

  New data frame

- random_terms:

  Parsed random effects terms from original model

- group_var:

  Name of grouping variable

## Value

Random effects design matrix Z
