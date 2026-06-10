# GLLAMMR 0.2.0 (Development - February 2026)

## Major New Features

### Multi-Level IRT Models ✨✨✨

* **Full support for hierarchical/clustered IRT data**
  - Nested structures: students > classes > schools > districts
  - Crossed effects: students × time, students × raters
  - Partially nested: some students in classes, others not (via NA)

* **lme4-style formula syntax**
  ```r
  # Single level
  fit <- fit_irt(responses, model = "2PL",
                 person_data = data, random = ~ (1 | class))

  # Nested
  fit <- fit_irt(responses, model = "Rasch",
                 person_data = data, random = ~ (1 | school/class))

  # Multiple explicit levels
  fit <- fit_irt(responses, model = "2PL",
                 person_data = data,
                 random = ~ (1 | state) + (1 | district) + (1 | school))

  # Crossed effects (longitudinal)
  fit <- fit_irt(responses, model = "2PL",
                 person_data = long_data,
                 random = ~ (1 | student) + (1 | time))
  ```

* **Comprehensive S3 methods**
  - `VarCorr()` - Extract variance components
  - `icc()` - Compute intraclass correlations
  - `ranef()` - Extract random effects by level
  - `abilities(composite = TRUE)` - Total abilities (person + group effects)
  - `coef(type = "random")` - Access random effects

* **Model support**
  - ✅ Dichotomous: Rasch, 2PL, 3PL with random effects
  - ✅ Polytomous: GRM, PCM, GPCM, NRM with random effects
  - 🔄 EIRT: Template ready, R integration in progress

* **Enhanced print output**
  - Displays variance components table
  - Shows intraclass correlations (ICCs)
  - Lists grouping variables and group counts

* **Example output:**
  ```
  Multi-Level IRT Model ( Rasch )

  Random Effects:
    Grouping variables: class_id
    Number of groups: 20

  Variance Components:
     Groups Variance Std.Dev
   class_id   0.3558  0.5965
     Person   0.7093  0.8422
   Residual   3.2899  1.8138

  Intraclass Correlations:
       Level       ICC
    class_id 0.0817
      Person 0.1629
  ```

* **Comprehensive test suite**
  - 50+ unit tests for parsing, fitting, and methods
  - Parameter recovery validation
  - Edge case handling (partial nesting, crossed effects)

## API Changes

### EIRT Model Specification (Breaking Change)

* **Removed `model = "LPCM"`** - Merged into PCM for consistency
  - Old: `fit_eirt(..., model = "LPCM", threshold_formula = ~ x)`
  - New: `fit_eirt(..., model = "PCM", threshold_formula = ~ x)`
  - PCM with `threshold_formula` automatically uses threshold regression (LPCM framework)
  - GPCM now also supports `threshold_formula`

* **Added `item_residuals` parameter** - Control item-level residuals
  - `item_residuals = TRUE` (default): LLTM + error (b_i = W×γ + ε_b)
  - `item_residuals = FALSE`: Pure LLTM (b_i = W×γ, no residuals)
  - Allows testing whether residuals are needed via LRT

* **Better documentation** for `discrimination_formula`
  - Already worked for 2PL, GRM, GPCM but was undocumented
  - Now fully documented with examples

**Migration:** Simply change `model = "LPCM"` to `model = "PCM"` in existing code.

## Major New Features

### Weights Support ✨

* **ALL model families now support weights** (frequency and probability weights)
  - GLMM (Gaussian, Binomial, Poisson)
  - IRT models (Rasch, 2PL, 3PL, GRM, PCM, GPCM, NRM)
  - EIRT (Explanatory IRT)
  - Ordinal models
  - Latent Class Analysis
  - Multinomial models
  - Survival models

* **Usage**: Simply add `weights` parameter to any fitting function
  ```r
  fit <- gllamm(y ~ x + (1 | group), data = data,
                family = binomial(), weights = data$weight)
  ```

* **Automatic validation**: Weights must be non-negative, correct length, and non-missing
* **IRT special handling**: Person-level weights automatically expanded to item-response level
* **New vignette**: `vignette("weights")` with comprehensive examples

### Marginal Predictions ✨

* **Population-averaged predictions via Monte Carlo integration**
  - Available for ALL model families
  - Computes E[Y|X] = ∫ g^(-1)(X'β + Z'u) f(u) du
  - Differs from conditional predictions for nonlinear links (Jensen's inequality)

* **New predict methods for all model types**
  - `predict.gllamm()` - Extended with `type = "marginal"`
  - `predict.gllamm_irt()` - IRT marginal predictions (NEW)
  - `predict.gllamm_eirt()` - EIRT marginal predictions with new item support (NEW)
  - `predict.gllamm_ordinal()` - Ordinal marginal probabilities (NEW)
  - Internal: `predict_multinomial()`, `predict_survival()` (NEW)

* **Key features**:
  - Adjustable precision via `n_sim` parameter (default: 1000)
  - Standard errors available via `se.fit = TRUE`
  - Optimized for Gaussian-identity link (no MC needed)
  - Support for new data predictions

* **Core utilities** (`R/marginal_utils.R`):
  - `mc_integrate_marginal()` - Monte Carlo integration engine
  - `extract_random_vcov()` - Random effects variance extraction
  - `get_inverse_link()` - Inverse link function retrieval

* **New vignette**: `vignette("marginal-predictions")` with mathematical details

### Implementation Details

* **8 TMB templates modified** for weighted likelihood
* **6 new R files** for predict methods (~1200 lines)
* **Comprehensive testing**: Core utilities fully tested, integration tests ready
* **NAMESPACE updated**: New S3 methods exported
* **Documentation complete**: Roxygen docs for all new functions

### Selective Guessing for 3PL ✨

* **NEW: `mc_items` parameter for 3PL models**
  - Specify which items have guessing parameters (multiple choice)
  - Other items use 2PL likelihood (no guessing)
  - Enables realistic modeling of mixed-format tests

* **Usage**:
  ```r
  # 20 items: first 15 are MC, last 5 are open-ended
  fit <- fit_irt(responses, model = "3PL", mc_items = 1:15)
  # OR with logical vector
  fit <- fit_irt(responses, model = "3PL", mc_items = c(rep(TRUE, 15), rep(FALSE, 5)))
  ```

* **Default**: `mc_items = NULL` applies guessing to all items (backward compatible)

### Mixed Item Types in EIRT ✨

* **PCM/GPCM now handle mixed dichotomous and polytomous items**
  - Binary items (K=2) are special case of PCM
  - Can mix 2-category and multi-category items in same assessment
  - Both `difficulty_formula` and `threshold_formula` apply appropriately

* **Example**:
  ```r
  # Assessment: 20 binary + 10 four-category items
  responses <- cbind(
    matrix(sample(0:1, 50*20, TRUE), 50, 20),  # Binary (0/1)
    matrix(sample(1:4, 50*10, TRUE), 50, 10)   # 4-category (1/2/3/4)
  )

  fit <- fit_eirt(responses, item_data,
                  difficulty_formula = ~ word_freq,
                  threshold_formula = ~ abstractness,
                  model = "PCM")  # Handles both!
  # Note: Binary items auto-recoded from 0/1 to 1/2 (message displayed)
  ```

* **Auto-recoding**: Binary items coded 0/1 are automatically recoded to 1/2
  - Polytomous models require 1-based coding (1, 2, ..., K)
  - 0-based items detected and recoded automatically with message

## Additional Enhancements

### IRT Models

* **Polytomous IRT models** now fully supported:
  - GRM (Graded Response Model)
  - PCM (Partial Credit Model)
  - GPCM (Generalized Partial Credit Model)
  - NRM (Nominal Response Model)

* **Explanatory IRT (EIRT)**:
  - Threshold-level predictors for polytomous models
  - Item parameter regression
  - DIF analysis
  - Predictions for new items using item covariates

### Bug Fixes

* **CRITICAL: Fixed threshold_formula being completely ignored** (introduced when removing LPCM)
  - After removing `model = "LPCM"`, code still checked for it
  - Result: `threshold_formula` was never applied to PCM/GPCM
  - Fixed: Now checks for PCM/GPCM models with threshold_formula
* Fixed TMB compilation issues for polytomous EIRT
* Improved convergence for complex random effects structures
* Fixed parameter extraction for EIRT models

---

# GLLAMMR 0.1.0 (Initial Development)

## Phase 1: Foundation & Basic GLMM

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
