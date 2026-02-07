# GLLAMMR 0.1.0 (Development)

## Phase 1: Foundation & Basic GLMM (Week 1 - Current)

### New Features

* Initial package setup with TMB backend
* Basic 2-level Gaussian GLMM with random intercepts
* lme4-style formula parsing (without lme4 dependency)
* S3 class system with `print()`, `summary()`, `coef()`, `vcov()`, `logLik()` methods
* Extract functions: `fixef()`, `ranef()`, `VarCorr()`
* Basic `predict()` and `simulate()` methods
* Comprehensive test suite with formula parsing, basic fitting, and simulation-recovery tests

### Infrastructure

* TMB C++ templates for Gaussian likelihood
* Formula parser supporting `(term | group)` syntax
* Model matrix construction
* TMB interface with optimization via `nlminb()`
* Input validation

### Documentation

* README with quick start guide
* Getting Started vignette stub
* Function documentation via roxygen2
* Implementation plan document

### Known Limitations (Current Phase)

* Only Gaussian family with identity link supported
* Only single random effects term (random intercept)
* No random slopes yet
* No other GLM families (binomial, Poisson, etc.)
* No multilevel (3+) models yet
* No factor models, IRT, or latent class models yet

## Roadmap

### Phase 2 (Weeks 5-8) - Enhanced GLMM
* Random slopes with `(x | group)` syntax
* Binomial and Poisson families
* Multiple levels (3+)
* Crossed random effects
* Adaptive quadrature
* Weights and offsets

### Phase 3 (Weeks 9-10) - Categorical Responses
* Ordinal responses (proportional odds, probit)
* Multinomial responses
* Rankings

### Phase 4 (Weeks 11-14) - Factor Models & IRT
* Factor analysis (CFA)
* IRT models: Rasch, 2PL, 3PL, GRM, PCM
* Multidimensional IRT
* Differential item functioning (DIF)

### Phase 5 (Weeks 15-17) - Latent Class Models
* Latent class analysis
* Growth mixture models
* Latent class IRT
* Nonparametric random effects

### Phase 6 (Weeks 18-20) - Mixed Response & SEM
* Mixed response models
* Structural equation models
* Joint models (longitudinal + survival)
* Higher-order factors

### Phase 7 (Weeks 21-22) - Advanced Features
* Censored/survival data
* Parameter constraints
* Missing data handling
* Advanced inference methods

### Phase 8 (Week 23) - Prediction & Post-Estimation
* Enhanced prediction methods
* Simulation utilities
* Diagnostic tools
* Visualization functions

### Phase 9 (Week 24) - Documentation & Release
* Comprehensive vignettes
* Complete function documentation
* CRAN submission
* Official release

---

# Version History Format

Future releases will follow semantic versioning (MAJOR.MINOR.PATCH):

* MAJOR version for incompatible API changes
* MINOR version for added functionality (backward compatible)
* PATCH version for bug fixes (backward compatible)
