# EIRT Design Inconsistencies and Proposed Fixes

## Issues Identified

### 1. Inconsistency: Rasch vs PCM Model Specification

**Current design (INCONSISTENT):**

**Dichotomous:**
```r
# Just one model, add predictors via formula
fit_eirt(responses, item_data,
         difficulty_formula = ~ word_freq,
         model = "Rasch")  # One model
```

**Polytomous:**
```r
# Two separate models!
fit_eirt(responses, item_data,
         difficulty_formula = ~ word_freq,
         model = "PCM")  # Without threshold predictors

fit_eirt(responses, item_data,
         difficulty_formula = ~ word_freq,
         threshold_formula = ~ abstractness,
         model = "LPCM")  # With threshold predictors - SEPARATE MODEL!
```

**Problem:** Why have separate PCM and LPCM when we don't have separate Rasch and "Linear Rasch"?

---

### 2. Issue: Can't Do Pure LLTM (No Error Term)

**Current implementation ALWAYS includes item residuals:**

```cpp
// From src/gllamm_eirt.hpp line 121
difficulty(j) = difficulty_pred + epsilon_b(j);
```

**This means:**
- b_i = W × γ + ε_b[i]  ← LLTM + error ✅
- b_i = W × γ           ← Pure LLTM ❌ (not possible!)

**Problem:** No way to fit pure LLTM where item covariates perfectly predict difficulty.

---

### 3. Question: Can You Predict Discrimination?

**Answer: YES!** (Already implemented but maybe not clear)

From `src/gllamm_eirt.hpp` lines 123-128:
```cpp
Type log_discrim_pred = 0.0;
for (int p = 0; p < p_disc; p++) {
  log_discrim_pred += delta(p) * W_discrimination(j, p);
}
discrimination(j) = exp(log_discrim_pred + epsilon_a(j));
```

**This works NOW:**
```r
# 2PL with discrimination predictors
fit_eirt(responses, item_data,
         difficulty_formula = ~ word_freq,
         discrimination_formula = ~ item_type,  # ← Works!
         model = "2PL")

# GPCM with discrimination predictors
fit_eirt(responses, item_data,
         difficulty_formula = ~ word_freq,
         discrimination_formula = ~ item_type,  # ← Works!
         model = "GPCM")
```

**Status:** ✅ Already implemented, just needs better documentation

---

### 4. Clarification: "Prediction at Response Level"

**Question:** What does "prediction should happen at the response level" mean?

**Current implementation:**
- Person-item-response is the observation unit
- Model: P(Y_ij = k | θ_j, item parameters)
- This IS at the response level

**If you mean:**
- Item parameters (b_i, a_i, τ_im) are functions of item covariates
- ✅ Yes, we do this

**Or if you mean:**
- Response-level predictors (like person × item interactions)
- ❌ No, our predictors are item-level only (in item_data)

**Need clarification on what you mean here.**

---

## Proposed Design Changes

### Change 1: Merge LPCM into PCM

**Remove "LPCM" as a separate model.**

**New design:**
```r
# PCM without threshold predictors
fit_eirt(responses, item_data,
         difficulty_formula = ~ word_freq,
         model = "PCM")

# PCM WITH threshold predictors (currently called LPCM)
fit_eirt(responses, item_data,
         difficulty_formula = ~ word_freq,
         threshold_formula = ~ abstractness,  # Optional!
         model = "PCM")
```

**Implementation:**
- Keep model = "PCM" (poly_model_type = 2)
- If `threshold_formula` is provided, use LPCM likelihood
- If `threshold_formula = NULL`, use standard PCM likelihood

**Benefits:**
- ✅ Consistent with dichotomous approach
- ✅ One model with optional extensions
- ✅ More intuitive

---

### Change 2: Add Option for Pure LLTM (No Residuals)

**Add parameter to control residuals:**

```r
fit_eirt(responses, item_data,
         difficulty_formula = ~ word_freq,
         discrimination_formula = ~ item_type,
         model = "2PL",
         item_residuals = TRUE)  # Default: LLTM + error
```

**When `item_residuals = FALSE`:**
- b_i = W_diff × γ (no ε_b)
- log(a_i) = W_disc × δ (no ε_a)
- Pure LLTM / pure discrimination model

**Implementation:**
```cpp
// Modified template
if (item_residuals == 1) {
  difficulty(j) = difficulty_pred + epsilon_b(j);
  discrimination(j) = exp(log_discrim_pred + epsilon_a(j));
} else {
  difficulty(j) = difficulty_pred;  // Pure LLTM
  discrimination(j) = exp(log_discrim_pred);
}
```

**Benefits:**
- ✅ Supports pure LLTM
- ✅ Can test whether residuals are needed
- ✅ More flexible

---

### Change 3: Better Documentation of Discrimination Predictors

**Current status:** Already works, just not well documented!

**Add to documentation:**
```r
#' @param discrimination_formula Formula for discrimination regression.
#'   For 2PL: log(a_i) = W_disc %*% delta + epsilon_a
#'   For GPCM: log(a_i) = W_disc %*% delta + epsilon_a
#'   For Rasch/PCM/LPCM: Discrimination fixed at 1 (this parameter ignored)
```

---

## Complete Proposed API

### Dichotomous Models

```r
# Rasch with item predictors + residuals (LLTM + error)
fit_eirt(responses, item_data,
         difficulty_formula = ~ word_freq + length,
         model = "Rasch",
         item_residuals = TRUE)  # Default

# Rasch with pure LLTM (no residuals)
fit_eirt(responses, item_data,
         difficulty_formula = ~ word_freq + length,
         model = "Rasch",
         item_residuals = FALSE)

# 2PL with discrimination predictors
fit_eirt(responses, item_data,
         difficulty_formula = ~ word_freq,
         discrimination_formula = ~ item_type,
         model = "2PL",
         item_residuals = TRUE)
```

### Polytomous Models

```r
# PCM (standard)
fit_eirt(responses, item_data,
         difficulty_formula = ~ word_freq,
         model = "PCM")

# PCM with threshold predictors (currently LPCM)
fit_eirt(responses, item_data,
         difficulty_formula = ~ word_freq,
         threshold_formula = ~ abstractness,  # Makes it "LPCM"
         model = "PCM")

# GPCM with discrimination predictors
fit_eirt(responses, item_data,
         difficulty_formula = ~ word_freq,
         discrimination_formula = ~ item_type,
         model = "GPCM")

# GPCM with discrimination AND threshold predictors
fit_eirt(responses, item_data,
         difficulty_formula = ~ word_freq,
         discrimination_formula = ~ item_type,
         threshold_formula = ~ abstractness,
         model = "GPCM")
```

---

## Model Matrix

**After changes:**

| Model | Difficulty Predictors | Discrimination Predictors | Threshold Predictors | Item Residuals |
|-------|----------------------|---------------------------|---------------------|----------------|
| Rasch | ✅ | ❌ (fixed = 1) | ❌ | ✅ Optional |
| 2PL | ✅ | ✅ | ❌ | ✅ Optional |
| PCM | ✅ | ❌ (fixed = 1) | ✅ Optional | ✅ Optional |
| GPCM | ✅ | ✅ | ✅ Optional | ✅ Optional |
| GRM | ✅ | ✅ | ❌ (ordered) | ✅ Optional |

**Note:** LPCM is removed - it's just PCM with threshold_formula specified!

---

## Implementation Steps

### Step 1: Merge LPCM into PCM

**File: R/eirt.R**

```r
# Current:
model <- match.arg(model, c("Rasch", "2PL", "GRM", "PCM", "GPCM", "LPCM"))

# New:
model <- match.arg(model, c("Rasch", "2PL", "GRM", "PCM", "GPCM"))

# Internal logic:
if (model == "PCM" && !is.null(threshold_formula)) {
  # Use LPCM likelihood (poly_model_type = 4)
  poly_model_type <- 4L
} else if (model == "PCM") {
  # Use standard PCM likelihood (poly_model_type = 2)
  poly_model_type <- 2L
}
```

**File: src/gllamm_eirt.hpp**

No changes needed - poly_model_type = 4 already implements LPCM!

### Step 2: Add item_residuals Parameter

**File: R/eirt.R**

```r
fit_eirt <- function(response_matrix,
                     item_data,
                     difficulty_formula = ~ 1,
                     discrimination_formula = ~ 1,
                     threshold_formula = NULL,
                     weights = NULL,
                     model = c("Rasch", "2PL", "GRM", "PCM", "GPCM"),
                     item_residuals = TRUE,  # NEW!
                     start = NULL,
                     control = list())
```

**File: src/gllamm_eirt.hpp**

```cpp
// Add to DATA section
DATA_INTEGER(item_residuals);  // 1 = include residuals, 0 = pure LLTM

// Modify parameter computation (line 121)
if (item_residuals == 1) {
  difficulty(j) = difficulty_pred + epsilon_b(j);
  discrimination(j) = exp(log_discrim_pred + epsilon_a(j));
} else {
  difficulty(j) = difficulty_pred;
  discrimination(j) = exp(log_discrim_pred);
}

// Modify priors (only include if item_residuals == 1)
if (item_residuals == 1) {
  for (int j = 0; j < n_items; j++) {
    nll -= dnorm(epsilon_b(j), Type(0.0), sigma_epsilon_b, true);
    nll -= dnorm(epsilon_a(j), Type(0.0), sigma_epsilon_a, true);
  }
}
```

### Step 3: Document Discrimination Predictors

**File: R/eirt.R**

Update documentation:
```r
#' @param discrimination_formula Formula for discrimination regression (log scale).
#'   Applies to 2PL and GPCM models.
#'   Model: log(a_i) = W_disc %*% delta + epsilon_a (if item_residuals = TRUE)
#'   For Rasch, PCM: Discrimination fixed at 1 (parameter ignored)
#'   For GRM: Discrimination varies but uses different parameterization
```

---

## Backward Compatibility

**Breaking changes:**
- ❌ `model = "LPCM"` will error (need to use `model = "PCM"` with `threshold_formula`)

**Migration:**
```r
# Old:
fit_lpcm <- fit_eirt(..., model = "LPCM", threshold_formula = ~ x)

# New:
fit_pcm <- fit_eirt(..., model = "PCM", threshold_formula = ~ x)
```

**Non-breaking additions:**
- ✅ `item_residuals` parameter (default = TRUE preserves current behavior)
- ✅ Better documentation

---

## Benefits of Proposed Design

### 1. Consistency
- Dichotomous and polytomous use same design philosophy
- One model with optional extensions via formulas

### 2. Flexibility
- Pure LLTM vs LLTM + error via `item_residuals`
- Optional threshold predictors via `threshold_formula`
- Optional discrimination predictors via `discrimination_formula`

### 3. Clarity
- Model name describes family (Rasch, PCM, etc.)
- Formulas specify what's being predicted
- No need for separate "LPCM" model

### 4. Power
- Can fit all variants:
  - Pure LLTM (no residuals)
  - LLTM + error (with residuals)
  - Threshold predictors (PCM → LPCM)
  - Discrimination predictors (PCM → GPCM → GPCM+LPCM)

---

## Response to Original Questions

### Q1: "Get rid of LPCM and integrate into PCM"
**A:** ✅ Agree! Proposed above.

### Q2: "For 2PL and GPCM, should be able to predict discrimination"
**A:** ✅ Already works! Just needs documentation.

### Q3: "Prediction should happen at response level"
**A:** ⚠️ Need clarification. Already at person-item-response level.

### Q4: "What about LLTM with no error?"
**A:** ❌ Not currently supported. Proposed `item_residuals = FALSE` parameter.

---

## Example Use Cases After Changes

### Pure LLTM (No Residuals)
```r
fit <- fit_eirt(responses, item_data,
                difficulty_formula = ~ word_freq + length,
                model = "Rasch",
                item_residuals = FALSE)

# Test if residuals needed
fit_with_error <- update(fit, item_residuals = TRUE)
anova(fit, fit_with_error)  # LRT
```

### PCM with Everything
```r
fit <- fit_eirt(responses, item_data,
                difficulty_formula = ~ word_freq,    # Item location
                threshold_formula = ~ abstractness,   # Threshold spacing (makes it "LPCM")
                model = "PCM")
```

### GPCM with All Predictors
```r
fit <- fit_eirt(responses, item_data,
                difficulty_formula = ~ word_freq,
                discrimination_formula = ~ item_type,
                threshold_formula = ~ abstractness,
                model = "GPCM")
```

---

## Summary

**Current issues:**
1. ❌ Inconsistency: Separate LPCM model
2. ❌ Missing: Pure LLTM option
3. ⚠️ Unclear: Discrimination predictors work but not documented
4. ❓ Unclear: What "response level" means

**Proposed fixes:**
1. ✅ Merge LPCM into PCM (use threshold_formula to activate)
2. ✅ Add item_residuals parameter for pure LLTM
3. ✅ Document discrimination_formula better
4. ❓ Clarify "response level" requirement

**Result:** More consistent, flexible, and powerful API! 🎉
