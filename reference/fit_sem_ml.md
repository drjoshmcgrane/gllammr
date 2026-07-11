# Covariance-based ML estimation for SEM (Wishart / FIML likelihood)

Fits the measurement + recursive structural model. With complete data
the data reduce to the sample covariance matrix and mean vector (Wishart
ML, the lavaan/LISREL approach - fitting cost independent of N); with
`missing = "fiml"` the casewise (missing-pattern) normal likelihood is
maximized directly, with means as free parameters.

## Usage

``` r
fit_sem_ml(
  Y,
  lambda_pattern,
  beta_pattern,
  theta_zero = NULL,
  missing = c("listwise", "fiml"),
  se = TRUE,
  control = list()
)
```

## Arguments

- Y:

  Observed-variable matrix (indicators first, then any observed
  structural covariates appended by the caller)

- lambda_pattern:

  p x q loading pattern (0 zero / 1 free / 2 fixed 1)

- beta_pattern:

  q x q structural pattern (1 = free path row ~ col)

- theta_zero:

  Logical p-vector: TRUE for rows whose residual variance is fixed at 0
  (covariate pseudo-indicators)

- missing:

  "listwise" or "fiml"

- se:

  Compute standard errors (numerical Hessian; default TRUE)

- control:

  Optimization control list

## Details

Model structure: y = nu + Lambda eta + epsilon; eta = B eta + zeta.
Exogenous latent variables (no incoming B paths - including the
pseudo-latents that carry observed structural covariates) have a free
covariance matrix parameterized by its Cholesky factor; endogenous
disturbances are uncorrelated.
