# gllammr: Generalized Linear Latent and Mixed Models

Comprehensive implementation of Generalized Linear Latent and Mixed
Models following the 'Stata' 'GLLAMM' framework by Rabe-Hesketh,
Skrondal, and Pickles (2004)
[doi:10.1007/BF02295939](https://doi.org/10.1007/BF02295939) . Supports
multilevel generalized linear models with Gaussian, binomial, Poisson,
ordinal, and multinomial responses; item response theory models
including dichotomous (Rasch, 2PL, 3PL) and polytomous (GRM, PCM, GPCM,
NRM) models; differential item functioning analysis; explanatory item
response models with item covariates; latent class analysis; structural
equation models; mixed response models; and survival models. All models
support frequency and probability weights for survey data and aggregated
observations. Provides marginal (population-averaged) predictions via
Monte Carlo integration. Uses 'TMB' (Template Model Builder) for
efficient computation via automatic differentiation and Laplace
approximation. Provides comprehensive diagnostics, visualization, and
model comparison tools.

## See also

Useful links:

- <https://github.com/drjoshmcgrane/gllammr>

- <https://drjoshmcgrane.github.io/gllammr/>

- Report bugs at <https://github.com/drjoshmcgrane/gllammr/issues>

## Author

**Maintainer**: Josh McGrane <drjoshmcgrane@gmail.com>
