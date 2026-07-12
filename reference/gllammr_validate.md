# Cross-package validation of gllammr estimates

Fits canonical benchmark datasets with gllammr and with established
reference packages, and reports the agreement. Reference packages that
use the same Laplace approximation (lme4 with nAGQ = 1, ordinal::clmm)
should agree to numerical precision; packages using different
integration schemes (mirt, ltm EM quadrature) agree within small
tolerances.

## Usage

``` r
gllammr_validate(
  cases = "all",
  scale = c("standard", "large", "all"),
  verbose = TRUE
)
```

## Arguments

- cases:

  Character vector of case names to run, or "all" (default). Available:
  "gaussian_sleepstudy", "binomial_toenail", "poisson_grouseticks",
  "ordinal_wine", "rasch_lsat", "twopl_simulated", "lca_carcinoma",
  "grm_science", "gamma_simulated", "survival_exponential",
  "sem_lavaan", "lca_polytomous", "npml_binomial", "aghq_binomial",
  "twopl_lsat_em", "eirt_verbagg", "eirt_verbagg_pcm",
  "cdm_fraction_dina", "ordinal_crossed", "dif_logistic",
  "dif_irt_glmm".

- scale:

  "standard" (default) runs the canonical-dataset cases; "large" runs
  the large-scale tier (n in the tens of thousands, long item
  batteries - sizes where quadrature grids and tolerances can fail
  silently); "all" runs both.

- verbose:

  Print progress messages (default TRUE)

## Value

Data frame with one row per compared statistic: case, statistic, gllammr
value, reference value, absolute and relative difference, tolerance, and
pass/fail.

## Details

All reference packages are Suggests; cases whose reference package is
not installed are skipped.

## Examples

``` r
# \donttest{
if (requireNamespace("lme4", quietly = TRUE)) {
  gllammr_validate(cases = "gaussian_sleepstudy")
}
#> Validating: gaussian_sleepstudy
#>                  case      statistic    gllammr  reference     abs_diff
#> 1 gaussian_sleepstudy beta_intercept  251.40510  251.40510 9.094947e-13
#> 2 gaussian_sleepstudy      beta_Days   10.46729   10.46729 2.433609e-13
#> 3 gaussian_sleepstudy         logLik -875.96967 -875.96967 5.456968e-12
#> 4 gaussian_sleepstudy  var_intercept  565.51532  565.51521 1.046647e-04
#> 5 gaussian_sleepstudy      var_slope   32.68217   32.68219 1.690561e-05
#> 6 gaussian_sleepstudy  cov_int_slope   11.05537   11.05537 2.590172e-06
#>       rel_diff tolerance pass note
#> 1 3.617646e-15     1e-04 TRUE     
#> 2 2.324966e-14     1e-04 TRUE     
#> 3 6.229631e-15     1e-03 TRUE     
#> 4 1.850786e-07     1e-02 TRUE     
#> 5 5.172728e-07     1e-02 TRUE     
#> 6 2.342908e-07     5e-02 TRUE     
# }
```
