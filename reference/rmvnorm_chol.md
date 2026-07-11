# Draw samples from a zero-mean multivariate normal

Cholesky-based sampler so MASS can remain in Suggests.

## Usage

``` r
rmvnorm_chol(n, Sigma)
```

## Arguments

- n:

  Number of samples

- Sigma:

  Variance-covariance matrix

## Value

n x ncol(Sigma) matrix of draws
