# Plot item response curves by DIF group

Plots the model-implied probability of a correct/positive response (or
the expected score, for polytomous items) against the matching
criterion, for each level of one DIF variable, from the per-item full
DIF model. Other DIF variables are held at their reference level.

## Usage

``` r
dif_plot(dif_result, item, by = NULL, ...)
```

## Arguments

- dif_result:

  Object from
  [`dif_test`](https://drjoshmcgrane.github.io/gllammr/reference/dif_test.md)

- item:

  Item index (position among the tested items)

- by:

  Name of the DIF variable to display (default: the first one)

- ...:

  Additional graphical parameters
