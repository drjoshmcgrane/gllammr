# Compute Intraclass Correlation Coefficients

Calculate ICCs for multi-level IRT models

## Usage

``` r
icc(x, ...)
```

## Arguments

- x:

  A fitted multi-level IRT model

- ...:

  Additional arguments (not used; methods may add `level`)

## Value

A named vector of ICCs, or a single ICC if level specified

## Examples

``` r
if (FALSE) { # \dontrun{
# All ICCs
icc(fit)

# Specific level
icc(fit, level = "class")
} # }
```
