# GLLAMMR vs Stata GLLAMM: Feature Comparison

**Date:** 7 Feb 2026
**Purpose:** Comprehensive comparison to ensure GLLAMMR has all critical Stata GLLAMM functionality

---

## Research Sources

Based on information from:
- [Stata GLLAMM Manual](https://www.stata.com/manuals14/rgllamm.pdf)
- [GLLAMM Primer by Stas Kolenikov](https://staskolenikov.net/stata/gllamm-demo.html)
- [Bristol CMM GLLAMM Review](https://www.bristol.ac.uk/cmm/learning/mmsoftware/gllamm.html)
- [Stata GLLAMM Features](https://www.stata.com/features/generalized-linear-models/)
- Various Statalist discussions on gllapred, weights, and robust SEs

---

## Core Model Families

| Family | Stata GLLAMM | GLLAMMR | Notes |
|--------|--------------|---------|-------|
| **Gaussian** | ✅ | ✅ | Continuous response |
| **Binomial** | ✅ | ✅ | Binary response |
| **Poisson** | ✅ | ✅ | Count data |
| **Negative Binomial** | ❌ | ❌ | Neither package has this |
| **Multinomial** | ✅ | ✅ | Unordered categorical |
| **Ordinal** | ✅ (logit/probit only) | ✅ (6 link functions!) | **GLLAMMR has more** |
| **Gamma** | ❓ | ❌ | Unclear in Stata docs |
| **Inverse Gaussian** | ❓ | ❌ | Unclear in Stata docs |

### GLLAMMR Extensions Beyond Stata

1. **IRT Models (7 types)**: Rasch, 2PL, 3PL, GRM, PCM, GPCM, NRM - Stata doesn't have these
2. **Explanatory IRT (EIRT)**: Item parameter regression - Stata doesn't have this
3. **Latent Class Analysis**: GLLAMMR has dedicated fit_lca() - Stata requires manual setup
4. **DIF Analysis**: Built-in dif_test() - Stata doesn't have this

---

## Link Functions

| Link Function | Stata GLLAMM | GLLAMMR | Notes |
|---------------|--------------|---------|-------|
| **Identity** | ✅ | ✅ | Linear regression |
| **Logit** | ✅ | ✅ | Standard for binary |
| **Probit** | ✅ | ✅ | Alternative for binary |
| **Complementary log-log (cloglog)** | ✅ | ❌ | **MISSING in GLLAMMR** |
| **Log** | ✅ | ✅ | For Poisson/count |
| **ACL (Adjacent Category)** | ❌ | ✅ | **GLLAMMR exclusive** |
| **CRL (Continuation Ratio)** | ❌ | ✅ | **GLLAMMR exclusive** |
| **PPO (Partial Proportional Odds)** | ❌ | ✅ | **GLLAMMR exclusive** |

### Critical Missing Link: cloglog

The **complementary log-log** link is important for:
- Asymmetric binary outcomes
- Survival analysis with discrete time
- Gompertz/extreme value distributions
- Common in epidemiology and survival analysis

**Recommendation:** ✅ **ADD cloglog link to binomial family**

---

## Numerical Integration

| Method | Stata GLLAMM | GLLAMMR | Notes |
|--------|--------------|---------|-------|
| **Gaussian Quadrature** | ✅ | ✅ (via TMB) | Standard approach |
| **Adaptive Quadrature** | ✅ (recommended) | ⚠️ | **GLLAMMR uses Laplace** |
| **Laplace Approximation** | ❓ | ✅ | GLLAMMR default via TMB |
| **Discrete (Non-parametric)** | ✅ | ❌ | Freely estimated point masses |

### Integration Comparison

**Stata GLLAMM:**
- `nip()` option: number of integration points per dimension
- `adapt` option: enables adaptive quadrature (repositions integration points)
- `ip(f)`: discrete factor model with freely estimated masses
- Computation time scales with product of points across dimensions

**GLLAMMR:**
- Uses TMB's Laplace approximation for random effects
- Automatic differentiation for efficiency
- No user control over integration method currently

**Assessment:**
- Laplace approximation is fast and accurate for many models
- Adaptive quadrature can be more accurate but slower
- Our approach is more modern (TMB automatic differentiation)
- **No critical gap** - different but valid approaches

---

## Model Types and Capabilities

### Multilevel/Hierarchical Models

| Feature | Stata GLLAMM | GLLAMMR | Status |
|---------|--------------|---------|--------|
| Nested random effects | ✅ | ✅ | Both support |
| Random coefficients/slopes | ✅ | ✅ | Both support |
| Multiple levels | ✅ | ✅ | Both support |
| Cross-classified | ✅ | ⚠️ | Unclear in GLLAMMR |

### Latent Variable Models

| Feature | Stata GLLAMM | GLLAMMR | Status |
|---------|--------------|---------|--------|
| Factor models | ✅ | ✅ | Both support |
| Structural equation models | ✅ | ⚠️ | fit_sem() not implemented |
| IRT models | ❌ | ✅ | **GLLAMMR exclusive** |
| Latent class models | ✅ (manual) | ✅ (dedicated) | **GLLAMMR easier** |

### Panel Data Models

| Feature | Stata GLLAMM | GLLAMMR | Status |
|---------|--------------|---------|--------|
| Repeated measurements | ✅ | ✅ | Both support |
| Unbalanced panels | ✅ | ✅ | Both support |

### Multiple Equation Models

| Feature | Stata GLLAMM | GLLAMMR | Status |
|---------|--------------|---------|--------|
| Different families per equation | ✅ | ❌ | **MISSING** |
| Different links per equation | ✅ | ❌ | **MISSING** |

**Stata Approach:** Uses `lv()`, `fv()` options to vary families/links by observation type. Requires manual data restructuring with `expand` and indicator variables.

**Assessment:** This is a complex feature. fit_mixed_response() TMB template exists but not implemented. **Lower priority** unless user needs it.

---

## Estimation Options and Features

### Starting Values

| Feature | Stata GLLAMM | GLLAMMR | Status |
|---------|--------------|---------|--------|
| User-specified starting values | ✅ `from()` | ✅ `start` parameter | Both support |
| Copy from previous model | ✅ `copy` option | ⚠️ | Can manually extract |
| Warm starts | ✅ | ✅ | Both support |

### Weights

| Feature | Stata GLLAMM | GLLAMMR | Status |
|---------|--------------|---------|--------|
| Probability weights (pweights) | ✅ `pweight()` | ❌ | **MISSING** |
| Frequency weights (fweights) | ✅ | ❌ | **MISSING** |
| Analytic weights (aweights) | ✅ | ❌ | **MISSING** |
| Sampling weights | ✅ Multilevel | ❌ | **MISSING** |

**Stata:** `pweight(stubname)` specifies stubname1, stubname2, etc. contain sampling weights for level 1, 2, etc.

**Assessment:** Weights are **important for survey data**. This is a **significant gap**.

**Recommendation:** ✅ **ADD weights support**
- Priority: HIGH for survey research applications
- Implementation: Add `weights` parameter to gllamm()
- TMB supports observation-level weights

### Robust Standard Errors

| Feature | Stata GLLAMM | GLLAMMR | Status |
|---------|--------------|---------|--------|
| Sandwich/Robust SEs | ✅ | ❓ | **Unclear** |
| Cluster-robust SEs | ✅ | ❓ | **Unclear** |

**Stata:** Provides robust SEs based on sandwich estimator (was notable before Stata 11).

**Assessment:** Need to check if TMB sdreport provides robust SE options.

**Recommendation:** ⚠️ **CHECK and potentially ADD**

### Constraints

| Feature | Stata GLLAMM | GLLAMMR | Status |
|---------|--------------|---------|--------|
| Linear constraints | ✅ `constr()` | ❌ | **MISSING** |
| Fixed parameters | ✅ | ⚠️ | Can use offsets |
| Equality constraints | ✅ | ❌ | **MISSING** |

**Stata:** Uses `constr()` option for model identification and parameter restrictions.

**Assessment:** Useful for identification in complex models. **Medium priority**.

**Recommendation:** ⚠️ **CONSIDER adding** (via TMB's MAP parameter feature)

### Variance Modeling

| Feature | Stata GLLAMM | GLLAMMR | Status |
|---------|--------------|---------|--------|
| Heteroscedastic errors | ✅ Scale equations | ❌ | **MISSING** |
| Constrained variances | ✅ | ⚠️ | Unclear |

**Stata:** Allows scale equations to model heteroscedastic error variance.

**Assessment:** Less common feature. **Low priority** unless requested.

---

## Post-Estimation and Prediction

### gllapred Command (Stata)

Stata's `gllapred` provides extensive prediction options:

| Option | Description | GLLAMMR Equivalent | Status |
|--------|-------------|-------------------|--------|
| **xb** | Fixed part of linear predictor | predict(..., type = "link") | ✅ |
| **u** | Posterior means & SDs of random effects | ranef() | ✅ |
| **linpred** | Linear predictor with posterior means | predict(..., type = "link") | ✅ |
| **mu** | Mean response E[g⁻¹(ν)] | predict(..., type = "response") | ✅ |
| **marginal** | Marginal predictions (integrated over REs) | ❌ | **MISSING** |
| **us(varname)** | Conditional predictions for specific RE values | ❌ | **MISSING** |

### gllasim Command (Stata)

Stata's `gllasim` simulates from fitted models:

| Option | Description | GLLAMMR Equivalent | Status |
|--------|-------------|-------------------|--------|
| Default | Simulate responses | simulate.gllamm() | ✅ |
| **u** | Simulate latent variables | ⚠️ | Partial |
| **us(varname)** | Simulate for specified REs | ❌ | **MISSING** |
| **from(matrix)** | Simulate with alternative parameters | ❌ | **MISSING** |

### GLLAMMR Prediction Methods

```r
# What we have
predict(fit, type = "response")  # Fitted values
predict(fit, type = "link")      # Linear predictor
fitted(fit)                      # Fitted values
residuals(fit)                   # Residuals
ranef(fit)                       # Random effects (posterior modes)
fixef(fit)                       # Fixed effects
simulate(fit)                    # Simulate responses
```

### Critical Missing Features

1. **Marginal vs Conditional Predictions**
   - Stata: Can integrate over RE distribution (`marginal`) or condition on specific RE values (`us()`)
   - GLLAMMR: Only conditional predictions (using posterior modes)
   - **Impact:** Important for population-averaged vs subject-specific inference
   - **Recommendation:** ✅ **ADD marginal prediction option**

2. **Posterior Standard Deviations**
   - Stata: Returns both posterior means AND SDs of random effects
   - GLLAMMR: ranef() returns only modes, not uncertainty
   - **Recommendation:** ⚠️ **ENHANCE ranef() to include SEs**

---

## Diagnostic and Model Comparison Features

### What Stata GLLAMM Has

| Feature | Available | Notes |
|---------|-----------|-------|
| Log-likelihood | ✅ | Standard |
| AIC/BIC | ✅ | Standard |
| LR tests | ✅ | Manual via lrtest |
| Score tests | ❓ | Not clearly documented |
| Modification indices | ❓ | Not clearly documented |

### What GLLAMMR Has

| Feature | Available | Notes |
|---------|-----------|-------|
| Log-likelihood | ✅ | logLik() method |
| AIC/BIC | ✅ | Standard |
| LR tests | ✅ | compare_eirt() for EIRT |
| Model-specific fit() | ✅ | **GLLAMMR exclusive** |
| Model-specific plot() | ✅ | **GLLAMMR exclusive** |
| DIF analysis | ✅ | **GLLAMMR exclusive** |
| EIRT utilities | ✅ | **GLLAMMR exclusive** |

**Assessment:** GLLAMMR has **more extensive diagnostics** for specific model types (IRT, EIRT, LCA).

---

## Summary: What's Missing from GLLAMMR

### High Priority Additions ⚠️

1. **Complementary log-log (cloglog) link**
   - Common for asymmetric binary outcomes and survival analysis
   - Easy to add to binomial family
   - **Estimated time:** 2-3 hours

2. **Weights support (pweights, fweights)**
   - Critical for survey data and sampling designs
   - TMB supports observation-level weights
   - **Estimated time:** 4-6 hours (implementation + testing)

3. **Marginal predictions**
   - Population-averaged inference (integrate over REs)
   - Important distinction from conditional predictions
   - **Estimated time:** 4-6 hours

### Medium Priority ⚠️

4. **Robust/Sandwich standard errors**
   - Check if TMB provides this
   - Add option to gllamm() if available
   - **Estimated time:** 2-4 hours (if TMB supports)

5. **Enhanced ranef() with standard errors**
   - Currently only returns posterior modes
   - Add posterior SDs for uncertainty quantification
   - **Estimated time:** 3-4 hours

6. **Parameter constraints**
   - Linear equality constraints
   - Useful for model identification
   - **Estimated time:** 6-8 hours (use TMB MAP feature)

### Low Priority (Optional)

7. **Multiple equation models with different families**
   - Complex feature requiring significant restructuring
   - fit_mixed_response() template exists but not implemented
   - **Estimated time:** 2-3 weeks (major undertaking)

8. **Heteroscedastic error modeling**
   - Scale equations for variance
   - Less commonly used
   - **Estimated time:** 1-2 weeks

9. **Discrete distribution for random effects**
   - Non-parametric ML with point masses
   - Alternative to normality assumption
   - **Estimated time:** 1-2 weeks

### Not Needed (Stata Doesn't Have Either)

- Negative binomial family (neither package has it)
- Many advanced diagnostics (both packages lack)

---

## What GLLAMMR Has That Stata Doesn't

### Major Advantages ✅

1. **Comprehensive IRT Support (7 models)**
   - Rasch, 2PL, 3PL, GRM, PCM, GPCM, NRM
   - Stata requires manual setup or separate commands

2. **Explanatory IRT (EIRT)**
   - Item parameter regression
   - Model comparison utilities
   - Completely unique to GLLAMMR

3. **Enhanced Ordinal Models**
   - 6 link functions (logit, probit, ACL, CRL forward/backward, PPO)
   - Stata only has logit/probit

4. **DIF Analysis**
   - Built-in dif_test(), dif_plot()
   - Stata requires manual implementation

5. **Model-Specific Diagnostics**
   - fit() methods for each model type
   - Rich model-specific statistics (entropy, pseudo-R², reliability, etc.)
   - Stata has only basic output

6. **Comprehensive Plotting**
   - plot() methods for each model type
   - ICC, IIF, TIF, class profiles, etc.
   - Stata requires manual graphing

7. **Modern TMB Backend**
   - Automatic differentiation
   - Fast and accurate
   - Laplace approximation efficient for many models

8. **Dedicated LCA Support**
   - fit_lca() with automatic restarts
   - Entropy, APPA, classification plots
   - Stata requires complex manual setup

---

## Recommended Implementation Plan

### Phase 1: Critical Missing Features (1-2 weeks)

**Goal:** Add features needed for feature parity with Stata GLLAMM

1. **Add cloglog link** (2-3 hours)
   - Modify binomial family to support cloglog
   - Update TMB templates if needed
   - Add tests

2. **Add weights support** (4-6 hours)
   - Add `weights` parameter to gllamm() and fit_* functions
   - Implement in TMB templates (observation-level weights)
   - Document and test with survey data examples

3. **Add marginal predictions** (4-6 hours)
   - Add `marginal = TRUE` option to predict()
   - Integrate over RE distribution (numerical integration)
   - Add examples comparing marginal vs conditional

4. **Check and add robust SEs** (2-4 hours)
   - Investigate TMB sdreport robust SE options
   - Add `robust = TRUE` parameter if available
   - Document when appropriate to use

5. **Enhance ranef() output** (3-4 hours)
   - Add standard errors/posterior SDs
   - Return both modes and uncertainty
   - Update documentation

**Total Phase 1 Time:** ~20-25 hours (2-3 weeks with testing and documentation)

### Phase 2: Enhancement Features (Optional)

6. **Parameter constraints** (6-8 hours)
   - Add `constraints` parameter using TMB MAP
   - Examples for identification

7. **Additional link functions** (2-3 hours each)
   - Consider: log-log, cauchit, others

### Phase 3: Advanced Features (Long-term)

8. **Multiple equation models** - Only if user requests (2-3 weeks)
9. **Heteroscedastic models** - Low demand (1-2 weeks)
10. **Discrete RE distributions** - Specialized use (1-2 weeks)

---

## Conclusion

### Overall Assessment

**GLLAMMR vs Stata GLLAMM:**

| Aspect | Winner | Notes |
|--------|--------|-------|
| **Basic GLMMs** | TIE | Both well-supported |
| **Ordinal models** | GLLAMMR | 6 links vs 2 |
| **IRT models** | GLLAMMR | 7 models vs none |
| **EIRT** | GLLAMMR | Unique feature |
| **LCA** | GLLAMMR | Easier interface |
| **Weights** | STATA | We lack this |
| **cloglog link** | STATA | We lack this |
| **Marginal predictions** | STATA | We lack this |
| **Diagnostics** | GLLAMMR | Much richer |
| **Plotting** | GLLAMMR | Built-in vs manual |
| **Modern backend** | GLLAMMR | TMB vs older methods |

### Critical Gaps to Address

1. ✅ **cloglog link** - Easy addition, commonly needed
2. ✅ **Weights support** - Important for survey data
3. ✅ **Marginal predictions** - Key for pop-averaged inference
4. ⚠️ **Robust SEs** - Check if TMB supports
5. ⚠️ **ranef() SEs** - Enhanced uncertainty quantification

### GLLAMMR's Unique Strengths

- **IRT/EIRT ecosystem** - Far beyond Stata capabilities
- **Rich diagnostics** - Model-specific fit() methods
- **Comprehensive plotting** - All model types
- **Modern inference** - TMB automatic differentiation
- **User-friendly** - Dedicated functions vs manual setup

### Bottom Line

**GLLAMMR is competitive with Stata GLLAMM** and in many ways **superior** for:
- IRT and educational measurement
- Ordinal data with flexible link functions
- Rich model diagnostics and visualization
- Modern computational backend

**Critical additions needed:**
- cloglog link (HIGH)
- Weights support (HIGH)
- Marginal predictions (HIGH)
- Robust SEs (MEDIUM)

With these additions, **GLLAMMR would match or exceed Stata GLLAMM** in nearly all respects.

---

**Report Date:** 7 Feb 2026
**Next Steps:** Implement Phase 1 critical features
