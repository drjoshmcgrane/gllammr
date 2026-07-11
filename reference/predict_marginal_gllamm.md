# Internal function for marginal predictions

Computes population-averaged predictions by integrating over random
effects

## Usage

``` r
predict_marginal_gllamm(object, newdata = NULL, n_sim = 1000, se.fit = FALSE)
```

## Arguments

- object:

  Fitted gllamm object

- newdata:

  New data for predictions (NULL = use original data)

- n_sim:

  Number of Monte Carlo samples

- se.fit:

  Return standard errors?

## Value

Vector of marginal predictions, or list with fit and se.fit
