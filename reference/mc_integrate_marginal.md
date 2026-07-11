# Monte Carlo integration for marginal predictions (one sample at a time)

More memory-efficient version that processes samples sequentially

## Usage

``` r
mc_integrate_marginal(X, Z, beta, Sigma_u, inv_link_fn, n_sim = 1000)
```

## Arguments

- X:

  Fixed effects design matrix

- Z:

  Random effects design matrix

- beta:

  Fixed effects coefficients

- Sigma_u:

  Random effects variance-covariance matrix

- inv_link_fn:

  Inverse link function

- n_sim:

  Number of Monte Carlo samples

## Value

List with fit and se
