# EM for the saturated multivariate-normal model under missingness

Standard EM on the expected sufficient statistics: the E-step fills in
conditional means (and adds the conditional covariance) per missing
pattern; the M-step is the sample mean / covariance of the completed
statistics. Used for the FIML saturated log-likelihood (fit indices).

## Usage

``` r
.mvn_em_saturated(Y, max_iter = 500, tol = 1e-08)
```
