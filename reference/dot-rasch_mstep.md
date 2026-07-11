# Weighted logistic M-step for the located latent class (Rasch) model

Given expected success counts and totals per item x class cell, fits
logit pi_jc = theta_c - delta_j by weighted binomial GLM (the
complete-data M-step of the latent class Rasch model of Lindsay, Clogg &
Grego 1991) and returns the fitted probability matrix.

## Usage

``` r
.rasch_mstep(num, den)
```
