# MML-EM estimation for IRT models (Bock-Aitkin)

Marginal maximum likelihood via the EM algorithm with fixed
Gauss-Hermite quadrature over the ability distribution - the algorithm
class used by mirt and TAM. The E-step computes person-by-node posterior
weights in one matrix product; M-steps are independent small
optimizations per item.

## Usage

``` r
fit_irt_em(
  response_matrix,
  model,
  weights = NULL,
  mc_items = NULL,
  quad_points = 61,
  max_iter = 500,
  tol = 1e-04,
  control = list()
)
```

## Details

Identification matches the Laplace path: 2PL/3PL/GRM/GPCM/NRM fix the
ability SD at 1; Rasch and PCM keep a free ability SD via the
common-slope equivalence (theta ~ N(0, sigma) with unit slope is theta ~
N(0,1) with shared slope sigma; difficulties back-transform as b = d \*
sigma).
