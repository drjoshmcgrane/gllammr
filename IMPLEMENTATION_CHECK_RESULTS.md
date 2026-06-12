# Implementation Check Results - Mixed Items & Selective Guessing

**Date:** February 9, 2026

---

## Issue 1: 3PL with Item-Specific Guessing Parameters

### Current Status: ❌ NOT IMPLEMENTED

**Problem:** All items get guessing parameters when using 3PL. No way to specify that only some items are multiple choice.

### Current Implementation

**File:** `R/irt.R` (line 161)
```r
# All items get guessing parameters
guessing_init <- rep(0.1, n_items)
```

**File:** `src/gllamm_irt.hpp` (line 57)
```cpp
// 3PL applies guessing to ALL items
Type c = guessing(item);
prob = c + (Type(1.0) - c) * invlogit(eta);
```

### What Needs to Be Added

**1. R Interface (R/irt.R)**

Add `mc_items` parameter:
```r
fit_irt <- function(response_matrix,
                    model = c("Rasch", "2PL", "3PL", "GRM", "PCM", "GPCM", "NRM"),
                    weights = NULL,
                    mc_items = NULL,  # NEW: Which items have guessing?
                    start = NULL,
                    control = list())
```

Pass to TMB:
```r
tmb_data <- list(
  ...,
  mc_items = as.integer(mc_indicator)  # NEW: 0/1 indicator vector
)
```

**2. TMB Template (src/gllamm_irt.hpp)**

Add data input:
```cpp
DATA_IVECTOR(mc_items);  // 1 if item has guessing, 0 otherwise
```

Modify 3PL likelihood:
```cpp
} else {
  // 3PL model - guessing only for MC items
  Type eta = discrimination(item) * (theta(person) - difficulty(item));

  if (mc_items(item) == 1) {
    // MC item: apply guessing parameter
    Type c = guessing(item);
    prob = c + (Type(1.0) - c) * invlogit(eta);
  } else {
    // Non-MC item: no guessing (same as 2PL)
    prob = invlogit(eta);
  }
}
```

**3. Parameter Initialization**

Only estimate guessing for MC items:
```r
# In R/irt.R
n_mc_items <- sum(mc_indicator)
guessing_init <- rep(0.1, n_mc_items)  # Only MC items
```

### Example Usage (After Implementation)

```r
# 20 items: first 15 are MC, last 5 are open-ended
responses <- matrix(rbinom(1000, 1, 0.6), 50, 20)

fit_3pl <- fit_irt(
  responses,
  model = "3PL",
  mc_items = c(1:15)  # Only these items get guessing parameter
)

# OR with logical vector
fit_3pl <- fit_irt(
  responses,
  model = "3PL",
  mc_items = c(rep(TRUE, 15), rep(FALSE, 5))
)
```

---

## Issue 2: Mixed Dichotomous/Polytomous Items in EIRT

### Current Status: ✅ MOSTLY WORKS, ⚠️ NEEDS BUG FIX

**Finding:** The infrastructure is already there, but there's a bug in the threshold formula handling!

### What Already Works ✅

**1. Per-Item Category Counting**

**File:** `R/eirt.R` (line 144-147)
```r
n_categories_per_item <- apply(response_matrix, 2, function(x) {
  length(unique(x[!is.na(x)]))
})
max_categories <- max(n_categories_per_item)
```

✅ Correctly computes categories per item
✅ Handles mixed items (some with 2 cats, some with 3+)

**2. TMB Template Uses Per-Item Categories**

**File:** `src/gllamm_eirt.hpp` (line 182)
```cpp
int K = n_categories_per_item(item);  // Get K for THIS specific item
```

✅ Each item uses its own K
✅ Loop over categories adapts: `for (int m = 1; m < K; m++)`
✅ Works for K=2, K=3, K=4, etc.

**3. Parameter Matrices Sized for max_categories**

**File:** `src/gllamm_eirt.hpp` (line 57, 60, 63)
```cpp
PARAMETER_MATRIX(step_param);   // [n_items x (max_categories-1)]
PARAMETER_MATRIX(xi);           // [p_thresh x (max_categories-1)]
PARAMETER_MATRIX(e_step);       // [n_items x (max_categories-1)]
```

✅ Sized for maximum
✅ Items with fewer categories just don't use all columns

### The Bug ❌

**File:** `R/eirt.R` (line 162-167)
```r
# Threshold regression design matrix (LPCM only)
if (model == "LPCM" && !is.null(threshold_formula)) {  # ← BUG: "LPCM" no longer exists!
  W_threshold <- model.matrix(threshold_formula, data = item_data)
} else {
  W_threshold <- matrix(0, n_items, 1)
}
```

**Problem:**
- We removed "LPCM" model
- But this code still checks for it
- Now `model == "LPCM"` is ALWAYS FALSE
- So threshold_formula is NEVER used!

### The Fix ✅

**File:** `R/eirt.R` (line 162-167)

**Replace this:**
```r
if (model == "LPCM" && !is.null(threshold_formula)) {
  W_threshold <- model.matrix(threshold_formula, data = item_data)
} else {
  W_threshold <- matrix(0, n_items, 1)
}
```

**With this:**
```r
# Threshold regression design matrix (for PCM/GPCM with threshold_formula)
if (!is.null(threshold_formula) && model %in% c("PCM", "GPCM")) {
  W_threshold <- model.matrix(threshold_formula, data = item_data)
} else {
  W_threshold <- matrix(0, n_items, 1)
}
```

### How Mixed Items Work with PCM

**Dichotomous items (K=2):**
- Category 0 (or 1) vs category 1 (or 2)
- 1 threshold: δ_i1
- PCM with K=2: `cumsum(1) = θ - δ_i1`
- Prob(Y=0) = 1/(1 + exp(θ - δ_i1))
- Prob(Y=1) = exp(θ - δ_i1)/(1 + exp(θ - δ_i1))
- **This is identical to Rasch/1PL!**

**Polytomous items (K>2):**
- Multiple categories
- K-1 thresholds: δ_i1, δ_i2, ..., δ_i(K-1)
- Standard PCM likelihood

**With EIRT:**
- **difficulty_formula**: Predicts b_i for ALL items (both 2-cat and multi-cat)
- **threshold_formula**: Predicts threshold spacing for ALL thresholds
  - 2-category items: 1 threshold
  - 3-category items: 2 thresholds
  - 4-category items: 3 thresholds
  - Each threshold gets predicted from item covariates

### Example: Mixed Items in EIRT

```r
# Vocabulary assessment: 30 items
# - 20 dichotomous (0/1): true/false questions
# - 10 polytomous (1/2/3/4): multiple choice with 4 options

# Response matrix (persons × items)
responses <- cbind(
  matrix(sample(0:1, 50*20, replace=TRUE), 50, 20),  # Dichotomous
  matrix(sample(1:4, 50*10, replace=TRUE), 50, 10)   # Polytomous
)

# Item-level covariates
item_data <- data.frame(
  word_freq = rnorm(30),
  abstractness = rnorm(30),
  n_cats = c(rep(2, 20), rep(4, 10))  # Informational
)

# Fit EIRT with mixed items
fit_mixed <- fit_eirt(
  response_matrix = responses,
  item_data = item_data,

  # Item difficulty: ALL items
  difficulty_formula = ~ word_freq + abstractness,

  # Threshold spacing: ALL thresholds
  # - Items 1-20: 1 threshold each (20 thresholds)
  # - Items 21-30: 3 thresholds each (30 thresholds)
  # - Total: 50 thresholds predicted
  threshold_formula = ~ abstractness,

  model = "PCM"  # PCM handles both!
)

# Results:
# - γ coefficients: Predict item location for all 30 items
# - ξ matrix: [p_thresh × 3] for max_categories-1
#   - ξ[,1]: Effect on 1st threshold (all items have this)
#   - ξ[,2]: Effect on 2nd threshold (only items 21-30 use this)
#   - ξ[,3]: Effect on 3rd threshold (only items 21-30 use this)
```

### Interaction Example

```r
# Same covariate affects both levels
fit_interact <- fit_eirt(
  responses,
  item_data,

  difficulty_formula = ~ word_freq,     # Overall item difficulty
  threshold_formula = ~ word_freq,      # Threshold spacing

  model = "PCM"
)

# Interpretation for a 4-category item:
# b_i = γ₀ + γ₁*word_freq_i + ε_b
# δ_i1 = b_i + (ξ₀1 + ξ₁1*word_freq_i) + ε_threshold
# δ_i2 = b_i + (ξ₀2 + ξ₁2*word_freq_i) + ε_threshold
# δ_i3 = b_i + (ξ₀3 + ξ₁3*word_freq_i) + ε_threshold

# word_freq has:
# - ONE effect on item location: γ₁
# - THREE effects on thresholds: ξ₁1, ξ₁2, ξ₁3 (different per threshold)
```

---

## Summary of Needed Changes

### Priority 1: Fix Threshold Formula Bug 🔴

**File:** `R/eirt.R` (line 162-167)

**Current (broken):**
```r
if (model == "LPCM" && !is.null(threshold_formula)) {
```

**Fixed:**
```r
if (!is.null(threshold_formula) && model %in% c("PCM", "GPCM")) {
```

**Impact:** Without this fix, threshold_formula is COMPLETELY IGNORED (always gets dummy matrix)

---

### Priority 2: Add mc_items Parameter for 3PL 🟡

**Files to modify:**
1. `R/irt.R` - Add mc_items parameter, create indicator vector
2. `src/gllamm_irt.hpp` - Add DATA_IVECTOR(mc_items), conditional guessing
3. Parameter initialization - Only estimate guessing for MC items

**Impact:** Enables realistic 3PL for mixed-format tests

---

### Priority 3: Add Validation & Documentation 🟢

**1. Validation for mixed items:**
```r
# In fit_eirt()
if (any(n_categories_per_item == 2) && any(n_categories_per_item > 2)) {
  message("Note: Data contains mixed dichotomous (", sum(n_categories_per_item == 2),
          ") and polytomous (", sum(n_categories_per_item > 2), ") items. ",
          "Using PCM which handles both as special cases of partial credit model.")
}
```

**2. Documentation:**
- Update fit_eirt() docs to explicitly mention mixed items work
- Add example with mixed items
- Clarify that PCM with K=2 is equivalent to Rasch

---

## Testing Recommendations

### Test 1: Mixed Items in EIRT

```r
# Create mixed data
set.seed(123)
responses_mixed <- cbind(
  matrix(sample(0:1, 100*10, replace=TRUE), 100, 10),  # 10 dichotomous
  matrix(sample(1:4, 100*5, replace=TRUE), 100, 5)     # 5 polytomous (4-cat)
)

item_data <- data.frame(
  item_id = 1:15,
  covar = rnorm(15)
)

# Should work after bug fix
fit <- fit_eirt(
  responses_mixed,
  item_data,
  difficulty_formula = ~ covar,
  threshold_formula = ~ covar,
  model = "PCM"
)

# Check:
# - Estimates should be reasonable
# - n_categories_per_item should be c(2,2,...,2,4,4,4,4,4)
# - No errors or warnings
```

### Test 2: 3PL with Selective Guessing (after implementation)

```r
# 20 items: 15 MC, 5 open-ended
responses_3pl <- matrix(rbinom(50*20, 1, 0.6), 50, 20)

fit <- fit_irt(
  responses_3pl,
  model = "3PL",
  mc_items = 1:15
)

# Check:
# - Only 15 guessing parameters estimated
# - Items 16-20 have no guessing (effectively 2PL)
```

---

*Check completed: February 9, 2026*
