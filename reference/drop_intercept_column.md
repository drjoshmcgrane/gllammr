# Drop the intercept column from a fixed-effects design matrix

Cumulative-link (ordinal) models absorb the location into the
thresholds; keeping a free intercept alongside free thresholds leaves
the model unidentified (only their difference enters the likelihood).

## Usage

``` r
drop_intercept_column(X)
```

## Arguments

- X:

  Design matrix from model.matrix()

## Value

X without its "(Intercept)" column
