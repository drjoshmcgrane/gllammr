# Monte Carlo integration for marginal predictions

Computes E\[g^(-1)(X'beta + Z'u)\] by averaging over draws from u ~ N(0,
Sigma_u)

## Usage

``` r
mc_integrate_fixed_samples(X, Z, beta, u_samples, inv_link_fn)
```

## Arguments

- X:

  Fixed effects design matrix (n x p)

- Z:

  Random effects design matrix (n x q) for a single observation

- beta:

  Fixed effects coefficients (p x 1)

- u_samples:

  Matrix of random effects samples (n_sim x q)

- inv_link_fn:

  Inverse link function

## Value

Vector of marginal predictions (length n)
