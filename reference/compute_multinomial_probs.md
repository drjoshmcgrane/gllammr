# Compute multinomial probabilities

Compute multinomial probabilities

## Usage

``` r
compute_multinomial_probs(X, beta, eta_random, n_categories)
```

## Arguments

- X:

  Fixed effects design matrix

- beta:

  Beta matrix (K-1) × p

- eta_random:

  Random effects contribution (vector of length n_obs)

- n_categories:

  Number of categories

## Value

Matrix of probabilities (n_obs × n_categories)
