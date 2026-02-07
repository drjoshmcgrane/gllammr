# IRT Enhancement Plan - GLLAMM Parity

## Current State
- ✅ Basic dichotomous IRT (Rasch, 2PL, 3PL) via `fit_irt()`
- ❌ Polytomous IRT models
- ❌ DIF analysis capabilities
- ❌ Explanatory IRT (item covariates)
- ❌ Integration with main `gllamm()` interface

## Required Enhancements

### 1. Polytomous IRT Models

#### A. Graded Response Model (GRM)
**Description**: Samejima's (1969) model for ordered categorical responses
**Parameterization**:
- Item discrimination: `a_j` (slope)
- Category thresholds: `b_{jk}` for k=1,...,K-1 categories
- Probability model: `P(Y_ij ≥ k | θ_i) = logit^{-1}(a_j(θ_i - b_{jk}))`

**Implementation**:
- New TMB template: `src/gllamm_irt_grm.hpp`
- Cumulative logit or probit link
- Constraint: `b_{j1} < b_{j2} < ... < b_{j,K-1}` (ordered thresholds)
- R interface: `fit_irt(..., model = "GRM")`

#### B. Partial Credit Model (PCM)
**Description**: Rasch model for polytomous items (Masters, 1982)
**Parameterization**:
- All items have equal discrimination (a_j = 1)
- Step difficulties: `δ_{jk}` for k=1,...,K-1
- Probability model using sequential logits

**Implementation**:
- New TMB template: `src/gllamm_irt_pcm.hpp`
- Constraint: Equal discriminations across items
- R interface: `fit_irt(..., model = "PCM")`

#### C. Generalized Partial Credit Model (GPCM)
**Description**: Extension of PCM with varying discriminations
**Parameterization**:
- Item-specific discriminations: `a_j`
- Step difficulties: `δ_{jk}`

**Implementation**:
- New TMB template: `src/gllamm_irt_gpcm.hpp`
- R interface: `fit_irt(..., model = "GPCM")`

#### D. Nominal Response Model (NRM)
**Description**: For unordered categorical responses
**Parameterization**:
- Category-specific slopes and intercepts
- No ordering constraint

**Implementation**:
- New TMB template: `src/gllamm_irt_nrm.hpp`
- R interface: `fit_irt(..., model = "NRM")`

### 2. Differential Item Functioning (DIF)

#### A. Uniform DIF (Difficulty DIF)
**Test**: Do item difficulties differ between groups?
**Method**:
```r
# Difficulty varies by group
b_j = b_j0 + b_j1 * GROUP
```
**Implementation**:
- Extend TMB templates to include group indicators
- Likelihood ratio test: constrained vs unconstrained model
- Wald test for DIF parameters
- Effect size measures (ETS Delta scale)

#### B. Non-uniform DIF (Discrimination DIF)
**Test**: Do item discriminations differ between groups?
**Method**:
```r
# Discrimination varies by group
a_j = a_j0 + a_j1 * GROUP
```

#### C. DIF Detection Functions
**R Interface**:
```r
dif_test(irt_fit, group_variable, items = NULL, type = c("both", "uniform", "nonuniform"))
# Returns:
# - Chi-square tests for each item
# - Effect sizes
# - Flagged items
# - Plots (item characteristic curves by group)
```

**Statistical Tests**:
- Likelihood ratio test
- Wald test
- Raju's area measures
- Lord's chi-square test

**Implementation**:
- New R file: `R/dif.R`
- Functions: `dif_test()`, `dif_plot()`, `dif_summary()`
- Integrate with existing IRT infrastructure

### 3. Explanatory IRT (EIRT) - Dichotomous Items

#### A. Item-Level Covariates
**Description**: Model item parameters as functions of item characteristics

**Difficulty Regression**:
```r
b_j = γ_0 + γ_1 * ITEM_COV1_j + γ_2 * ITEM_COV2_j + ε_j
```

**Discrimination Regression** (for 2PL/3PL):
```r
log(a_j) = δ_0 + δ_1 * ITEM_COV1_j + δ_2 * ITEM_COV2_j + η_j
```

**Example Use Cases**:
- Test whether item difficulty relates to word frequency
- Model discrimination as function of item type (multiple choice vs. open-ended)
- Explain guessing based on number of response options

**Implementation**:
- New TMB template: `src/gllamm_eirt_dichot.hpp`
- R interface:
```r
fit_eirt(response ~ person_covariates,
         item_difficulty ~ item_covariate1 + item_covariate2,
         item_discrimination ~ item_type,
         data = data,
         item_data = item_characteristics,
         model = "2PL")
```

#### B. Person-Level Covariates
**Description**: Already partially supported via fixed effects in main formula
```r
# Person ability depends on covariates
θ_i = β_0 + β_1 * AGE_i + β_2 * GENDER_i + u_i
```

**Enhancement**: Make this explicit in IRT context with proper identification constraints

### 4. Explanatory IRT - Polytomous Items

#### A. Threshold Regression (for GRM)
**Model**:
```r
b_{jk} = γ_{0k} + γ_{1k} * ITEM_COV_j + ε_{jk}
# with constraint: b_{j1} < b_{j2} < ... < b_{j,K-1}
```

**Implementation**:
- New TMB template: `src/gllamm_eirt_polytomous.hpp`
- Parameterize using cumulative sums to ensure ordering
- R interface:
```r
fit_eirt(response ~ person_covariates,
         item_thresholds ~ item_covariate,
         data = data,
         item_data = item_characteristics,
         model = "GRM")
```

#### B. Step Difficulty Regression (for PCM/GPCM)
**Model**:
```r
δ_{jk} = τ_{0k} + τ_{1k} * ITEM_COV_j + ζ_{jk}
```

### 5. Integration with Main GLLAMM Interface

#### A. Unified Syntax
**Goal**: Use `gllamm()` for IRT models with special syntax

**Proposed Syntax**:
```r
# Basic IRT (current)
fit_irt(responses, model = "2PL")

# IRT via gllamm() (new)
gllamm(response ~ 0 + item + (1 | person),
       family = binomial(),
       irt = list(model = "2PL", discrimination = "item"))

# Polytomous IRT
gllamm(response ~ 0 + item + (1 | person),
       family = ordinal(link = "logit"),
       irt = list(model = "GRM"))

# Explanatory IRT
gllamm(response ~ age + gender + (1 | person),
       family = binomial(),
       irt = list(
         model = "2PL",
         item_formula = list(
           difficulty ~ word_frequency,
           discrimination ~ item_type
         )
       ),
       item_data = item_characteristics)

# DIF analysis
gllamm(response ~ 0 + item*group + (1 | person),
       family = binomial(),
       irt = list(model = "2PL", dif = "group"))
```

#### B. Stata GLLAMM Correspondence
**Reference**: Stata GLLAMM manual section on IRT models

**Stata syntax** (example):
```stata
* 2PL model
gllamm y x1 x2, i(id) link(logit) family(binom) nrf(1) eqs(disc) nip(10)

* Explanatory IRT
gllamm y, i(id) link(logit) family(binom) nrf(1) ///
  expand(item it1 it2) load1(1-10) eqs(disc)
```

**Our equivalent**:
```r
gllamm(y ~ x1 + x2 + (1 | id), family = binomial(),
       irt = list(model = "2PL", nip = 10))
```

### 6. Implementation Priority

**Phase 1: Polytomous IRT (2 weeks)**
- [ ] Implement GRM (src/gllamm_irt_grm.hpp, update R/irt.R)
- [ ] Implement PCM (src/gllamm_irt_pcm.hpp)
- [ ] Implement GPCM (src/gllamm_irt_gpcm.hpp)
- [ ] Add tests (tests/testthat/test-irt-polytomous.R)
- [ ] Update vignettes/irt-models.Rmd

**Phase 2: DIF Analysis (1 week)**
- [ ] Create R/dif.R with DIF detection functions
- [ ] Implement likelihood ratio tests
- [ ] Implement Wald tests and effect sizes
- [ ] Add DIF plotting functions
- [ ] Add tests (tests/testthat/test-dif.R)
- [ ] Add DIF section to vignette

**Phase 3: Explanatory IRT - Dichotomous (1.5 weeks)**
- [ ] Create src/gllamm_eirt_dichot.hpp
- [ ] Implement item-level covariate regression
- [ ] Update R/irt.R with fit_eirt()
- [ ] Add tests (tests/testthat/test-eirt-dichot.R)
- [ ] Update vignette with EIRT examples

**Phase 4: Explanatory IRT - Polytomous (1.5 weeks)**
- [ ] Create src/gllamm_eirt_polytomous.hpp
- [ ] Implement threshold regression for GRM
- [ ] Implement step difficulty regression for PCM/GPCM
- [ ] Add tests (tests/testthat/test-eirt-polytomous.R)
- [ ] Update vignette

**Phase 5: Integration & Validation (2 weeks)**
- [ ] Integrate IRT into main gllamm() interface
- [ ] Match Stata GLLAMM on all IRT examples
- [ ] Match mirt on polytomous models
- [ ] Comprehensive testing and validation
- [ ] Final documentation updates

**Total: 8 weeks**

## Success Criteria

### Functionality
- [ ] All 7 major IRT model types working (Rasch, 2PL, 3PL, GRM, PCM, GPCM, NRM)
- [ ] DIF detection with statistical tests and plots
- [ ] Explanatory IRT for dichotomous items
- [ ] Explanatory IRT for polytomous items
- [ ] Integration with main gllamm() interface

### Validation
- [ ] Match mirt package on all comparable models (tolerance < 2%)
- [ ] Match TAM package on polytomous models
- [ ] Match Stata GLLAMM on IRT examples from gllamm.org
- [ ] Simulation-recovery tests for all model types

### Documentation
- [ ] Updated IRT vignette with all model types
- [ ] DIF analysis tutorial
- [ ] Explanatory IRT examples
- [ ] Stata GLLAMM migration guide updated for IRT

### Testing
- [ ] 40+ new tests covering all IRT functionality
- [ ] Edge cases (boundary parameters, sparse data, missing responses)
- [ ] Performance benchmarks

## File Structure After Enhancements

```
GLLAMMR/
├── src/
│   ├── gllamm_irt.hpp                    # Existing: Rasch/2PL/3PL
│   ├── gllamm_irt_grm.hpp                # NEW: Graded Response Model
│   ├── gllamm_irt_pcm.hpp                # NEW: Partial Credit Model
│   ├── gllamm_irt_gpcm.hpp               # NEW: Gen. Partial Credit
│   ├── gllamm_irt_nrm.hpp                # NEW: Nominal Response
│   ├── gllamm_eirt_dichot.hpp            # NEW: Explanatory IRT (binary)
│   └── gllamm_eirt_polytomous.hpp        # NEW: Explanatory IRT (polytomous)
│
├── R/
│   ├── irt.R                              # ENHANCED: Add polytomous models
│   ├── eirt.R                             # NEW: Explanatory IRT interface
│   └── dif.R                              # NEW: DIF analysis functions
│
├── tests/testthat/
│   ├── test-irt.R                         # EXISTING: Basic IRT tests
│   ├── test-irt-polytomous.R              # NEW: Polytomous IRT
│   ├── test-dif.R                         # NEW: DIF tests
│   ├── test-eirt-dichot.R                 # NEW: EIRT dichotomous
│   └── test-eirt-polytomous.R             # NEW: EIRT polytomous
│
└── vignettes/
    └── irt-models.Rmd                     # ENHANCED: Complete IRT guide
```

## References

### IRT Models
- Samejima, F. (1969). Graded Response Model. *Psychometrika*.
- Masters, G. N. (1982). Partial Credit Model. *Psychometrika*.
- Muraki, E. (1992). Generalized Partial Credit Model. *Applied Psychological Measurement*.
- Bock, R. D. (1972). Nominal Response Model.

### DIF Analysis
- Raju, N. S. (1990). Area measures for DIF.
- Lord, F. M. (1980). Chi-square test for DIF.
- Thissen, D., Steinberg, L., & Wainer, H. (1993). Detection of DIF using IRT.

### Explanatory IRT
- De Boeck, P., & Wilson, M. (2004). *Explanatory Item Response Models*.
- Rijmen, F., et al. (2003). Explaining item and person characteristics.

### Software
- Stata GLLAMM Manual: http://www.gllamm.org/
- mirt package: Chalmers (2012)
- TAM package: Kiefer, Robitzsch, & Wu (2021)
