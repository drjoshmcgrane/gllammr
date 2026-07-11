# Structural Equation Models

## Overview

[`fit_sem()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_sem.md)
(equivalently `gllamm(d, family = sem(...))`) fits structural equation
models with continuous indicators: confirmatory factor analysis,
recursive path models among latent variables, MIMIC models (latent
variables regressed on observed covariates), and full-information
maximum likelihood for missing data — with the standard fit-index
battery and lavaan-matching results (the validation suite pins
log-likelihoods, estimates, standard errors, and CFI/TLI/RMSEA/SRMR
against lavaan).

Estimation is Wishart maximum likelihood on the sample covariance matrix
for complete data — cost independent of sample size — and a casewise
(missing-pattern) likelihood under FIML.

## Confirmatory factor analysis

Two correlated factors, three indicators each. The first indicator of
each factor is the marker (loading fixed at 1); exogenous factors covary
freely, as in lavaan.

``` r

library(gllammr)
#> 
#> Attaching package: 'gllammr'
#> The following object is masked from 'package:stats':
#> 
#>     binomial
set.seed(2026)
n <- 500
f1 <- rnorm(n)
f2 <- 0.6 * f1 + sqrt(1 - 0.36) * rnorm(n)   # factor correlation 0.6
d <- data.frame(
  x1 = f1 + rnorm(n, 0, .6),  x2 = 0.8 * f1 + rnorm(n, 0, .6),
  x3 = 1.2 * f1 + rnorm(n, 0, .6),
  y1 = f2 + rnorm(n, 0, .5),  y2 = 0.9 * f2 + rnorm(n, 0, .5),
  y3 = 1.1 * f2 + rnorm(n, 0, .5))

cfa <- fit_sem(
  measurement = list(verbal = ~ x1 + x2 + x3,
                     quant  = ~ y1 + y2 + y3),
  data = d)
print(cfa)
#> Structural Equation Model (continuous indicators)
#> 
#> Observations: 500
#> 
#> Loadings (marker indicators fixed at 1):
#>    verbal  quant
#> x1 1.0000 0.0000
#> x2 0.8152 0.0000
#> x3 1.1982 0.0000
#> y1 0.0000 1.0000
#> y2 0.0000 0.8834
#> y3 0.0000 1.0878
#> 
#> Latent (co)variances:
#>        verbal  quant
#> verbal 1.0572 0.5716
#> quant  0.5716 0.9402
#> 
#> chisq(8) = 3.692, p = 0.884 | CFI = 1.000 | TLI = 1.004
#> RMSEA = 0.000 [0.000, 0.025] | SRMR = 0.009
#> 
#> Log-likelihood: -3503.99
```

The header reports the global fit: a non-significant chi-square, CFI/TLI
near 1, RMSEA near 0 (with its 90% interval), and a small SRMR all
indicate the two-factor structure reproduces the observed covariances.
[`summary()`](https://rdrr.io/r/base/summary.html) adds the full
parameter table — estimates, standard errors, z statistics, and the
standardized (std.all) solution, where loadings become indicator-factor
correlations and `verbal~~quant` the factor correlation:

``` r

summary(cfa)
#> Structural Equation Model (continuous indicators)
#> 
#> Observations: 500
#> 
#> Loadings (marker indicators fixed at 1):
#>    verbal  quant
#> x1 1.0000 0.0000
#> x2 0.8152 0.0000
#> x3 1.1982 0.0000
#> y1 0.0000 1.0000
#> y2 0.0000 0.8834
#> y3 0.0000 1.0878
#> 
#> Latent (co)variances:
#>        verbal  quant
#> verbal 1.0572 0.5716
#> quant  0.5716 0.9402
#> 
#> chisq(8) = 3.692, p = 0.884 | CFI = 1.000 | TLI = 1.004
#> RMSEA = 0.000 [0.000, 0.025] | SRMR = 0.009
#> 
#> Log-likelihood: -3503.99 
#> 
#> Parameter estimates:
#>           label    est     se       z pvalue est_std
#>      verbal=~x2 0.8152 0.0351 23.1923      0  0.8377
#>      verbal=~x3 1.1982 0.0475 25.2432      0  0.9055
#>       quant=~y2 0.8834 0.0327 27.0284      0  0.8720
#>       quant=~y3 1.0878 0.0379 28.7172      0  0.9054
#>  verbal~~verbal 1.0572 0.0909 11.6275      0  1.0000
#>   verbal~~quant 0.5716 0.0582  9.8154      0  0.5733
#>    quant~~quant 0.9402 0.0752 12.5017      0  1.0000
#>          x1~~x1 0.3724 0.0347 10.7460      0  0.2605
#>          x2~~x2 0.2985 0.0253 11.7970      0  0.2982
#>          x3~~x3 0.3332 0.0414  8.0395      0  0.1800
#>          y1~~y1 0.2377 0.0236 10.0632      0  0.2018
#>          y2~~y2 0.2313 0.0204 11.3380      0  0.2396
#>          y3~~y3 0.2446 0.0263  9.2959      0  0.1802
```

## Structural and MIMIC models

Regressions among latent variables, and of latent variables on observed
covariates, go in `structural`. Mixing the two is the MIMIC model:

``` r

d$ses <- rnorm(n)
f2b <- 0.5 * f1 + 0.3 * d$ses + rnorm(n, 0, .75)   # quant depends on ses
d$y1 <- f2b + rnorm(n, 0, .5)
d$y2 <- 0.9 * f2b + rnorm(n, 0, .5)
d$y3 <- 1.1 * f2b + rnorm(n, 0, .5)

mimic <- fit_sem(
  measurement = list(verbal = ~ x1 + x2 + x3,
                     quant  = ~ y1 + y2 + y3),
  structural = list(quant ~ verbal + ses),
  data = d)
mimic$param_table[grepl("~", mimic$param_table$label) &
                    !grepl("~~|=~", mimic$param_table$label), ]
#>          label       est         se        z       pvalue
#> 5 quant~verbal 0.3856450 0.04070434 9.474298 2.685616e-21
#> 6    quant~ses 0.2782425 0.03821528 7.280922 3.315461e-13
```

Observed covariates are handled in the joint-normal formulation, so the
results are likelihood-equivalent to lavaan with `fixed.x = FALSE`.

## Missing data: FIML

The default is listwise deletion (with a warning). With
`missing = "fiml"` every observed value contributes through the
missing-pattern likelihood — the assumption is missing at random:

``` r

d_na <- d
set.seed(1)
d_na$x1[sample(n, 60)] <- NA
d_na$y2[sample(n, 60)] <- NA

fiml <- fit_sem(
  measurement = list(verbal = ~ x1 + x2 + x3,
                     quant  = ~ y1 + y2 + y3),
  structural = list(quant ~ verbal),
  data = d_na, missing = "fiml")
c(n_used = fiml$n_obs, converged = fiml$convergence$converged)
#>    n_used converged 
#>       500         1
fiml$fit_measures[c("chisq", "df", "cfi", "rmsea")]
#>       chisq          df         cfi       rmsea 
#> 15.59460725  8.00000000  0.99591441  0.04357352
```

All 500 rows are retained, factor scores are computed casewise from each
person’s observed variables, and the fit indices use an EM-estimated
saturated model. (One convention note: under FIML our SRMR uses the
EM-saturated covariance, which can differ from lavaan’s in the third
decimal; everything else matches exactly.)

## Reading the fit measures

``` r

round(cfa$fit_measures, 4)
#>          chisq             df         pvalue            cfi            tli 
#>         3.6925         8.0000         0.8837         1.0000         1.0037 
#>          rmsea rmsea_ci_lower rmsea_ci_upper           srmr 
#>         0.0000         0.0000         0.0251         0.0091
```

- `chisq`/`df`/`pvalue`: exact-fit test against the saturated model.
- `cfi`, `tli`: incremental fit vs the independence baseline (rules of
  thumb: \> 0.95 good).
- `rmsea` with 90% CI: misfit per degree of freedom (\< 0.06 good).
- `srmr`: average standardized covariance residual (\< 0.08 good).

Inequality-free SEMs here are standard likelihood territory, so AIC/BIC
comparisons across non-nested structures are also reported.

## Notes on scope

Indicators are continuous; binary/ordinal measurement models belong to
the IRT side of the package
([`vignette("irt-models")`](https://drjoshmcgrane.github.io/gllammr/articles/irt-models.md)),
and joint models mixing indicator types share latent variables through
[`mixed_response()`](https://drjoshmcgrane.github.io/gllammr/reference/mixed_response.md).
The structural model must be recursive (no feedback loops). The legacy
full-data Laplace path (`method = "laplace"`) remains for latent-only
models but treats exogenous factors as uncorrelated and is no longer the
default.

## References

- Bollen, K. A. (1989). *Structural Equations with Latent Variables*.
  Wiley.
- Skrondal, A., & Rabe-Hesketh, S. (2004). *Generalized Latent Variable
  Modeling*. Chapman & Hall/CRC.
- Rosseel, Y. (2012). lavaan: an R package for structural equation
  modeling. *Journal of Statistical Software*, 48(2).
