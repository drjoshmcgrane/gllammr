# Internal function for marginal IRT predictions

Computes E\[P(Y=1\|θ)\] where θ ~ N(0, σ²_θ)

## Usage

``` r
predict_marginal_irt(object, items, n_sim = 1000)
```

## Arguments

- object:

  Fitted IRT model

- items:

  Item indices to predict

- n_sim:

  Number of MC samples

## Value

Vector of marginal probabilities (one per item)
