# EIRT Predictor Capabilities - Current Implementation

## Summary

After reviewing the polytomous EIRT implementation, here's what predictor support exists for each model:

---

## Dichotomous Models

### Rasch Model
- ✅ **Item-level predictors**: `difficulty_formula`
- ❌ **Discrimination predictors**: N/A (discrimination = 1, fixed)

### 2PL Model
- ✅ **Item-level predictors**: `difficulty_formula`
- ✅ **Discrimination predictors**: `discrimination_formula`

---

## Polytomous Models

### GRM (Graded Response Model)
- ✅ **Item-level predictors**: `difficulty_formula` (item location)
- ✅ **Discrimination predictors**: `discrimination_formula`
- ❌ **Threshold-level predictors**: N/A (thresholds are free parameters with ordering constraints)

**Code reference:** `src/gllamm_eirt.hpp` lines 191, 248
```cpp
Type a = discrimination(item);  // Uses discrimination from predictors
```

### PCM (Partial Credit Model)
- ✅ **Item-level predictors**: `difficulty_formula` (item location)
- ❌ **Discrimination predictors**: N/A (discrimination = 1, Rasch family)
- ❌ **Threshold-level predictors**: N/A (step deviations are free parameters)

### GPCM (Generalized Partial Credit Model)
- ✅ **Item-level predictors**: `difficulty_formula` (item location)
- ✅ **Discrimination predictors**: `discrimination_formula`
- ❌ **Threshold-level predictors**: N/A (step deviations are free parameters)

**Code reference:** `src/gllamm_eirt.hpp` line 248
```cpp
Type a = discrimination(item);  // Uses discrimination from predictors
```

### LPCM (Linear Partial Credit Model)
- ✅ **Item-level predictors**: `difficulty_formula` (item location b_i)
- ✅ **Threshold-level predictors**: `threshold_formula` (threshold-specific effects)
- ❌ **Discrimination predictors**: N/A (discrimination = 1, Rasch family)

**Code reference:** `src/gllamm_eirt.hpp` lines 287-291
```cpp
Type delta_im = difficulty(item);  // b_i from difficulty_formula
for (int p = 0; p < p_thresh; p++) {
  delta_im += W_threshold(item, p) * xi(p, m - 1);  // threshold predictors
}
delta_im += e_step(item, m - 1);  // threshold residual
```

**Model:** delta_im = b_i + sum_k xi_{k,m} * x_{i,k} + e_{i,m}

---

## Answer to Your Questions

### Question 1: Threshold AND item-level predictors for all polytomous models?

**Answer:** ❌ **Only LPCM has both**

- **LPCM**: ✅ Has both item-level (`difficulty_formula`) and threshold-level (`threshold_formula`) predictors
- **GRM, PCM, GPCM**: ❌ Only have item-level predictors; thresholds/steps are free parameters

This is by design, following the Kim & Wilson (2019) framework where LPCM is specifically designed for threshold-difficulty regression.

### Question 2: Discrimination predictors for models with varying discrimination?

**Answer:** ✅ **Partially - for GRM and GPCM only**

Models with varying discrimination:
- **2PL** (dichotomous): ✅ Has `discrimination_formula`
- **GRM** (polytomous): ✅ Has `discrimination_formula`
- **GPCM** (polytomous): ✅ Has `discrimination_formula`

Models with fixed discrimination (Rasch family):
- **Rasch** (dichotomous): ❌ discrimination = 1 (by definition)
- **PCM** (polytomous): ❌ discrimination = 1 (Rasch family)
- **LPCM** (polytomous): ❌ discrimination = 1 (Rasch family)

---

## Feature Matrix

| Model | Item-Level Difficulty | Threshold-Level | Discrimination | Family |
|-------|----------------------|-----------------|----------------|--------|
| Rasch | ✅ | N/A | ❌ (fixed=1) | Rasch |
| 2PL | ✅ | N/A | ✅ | 2PL |
| GRM | ✅ | ❌ (free params) | ✅ | 2PL-like |
| PCM | ✅ | ❌ (free params) | ❌ (fixed=1) | Rasch |
| GPCM | ✅ | ❌ (free params) | ✅ | 2PL-like |
| LPCM | ✅ | ✅ | ❌ (fixed=1) | Rasch |

---

## What's NOT Currently Implemented

### LGPCM (Hypothetical Extension)

A **Linear Generalized Partial Credit Model** would combine:
- ✅ Item-level difficulty predictors (`difficulty_formula`)
- ✅ Threshold-level predictors (`threshold_formula`)
- ✅ Discrimination predictors (`discrimination_formula`)

**This is NOT currently implemented.**

It would be the discrimination-varying version of LPCM, with model:
```
delta_im = b_i + sum_k xi_{k,m} * x_{i,k} + e_{i,m}
cumsum(m) = cumsum(m-1) + a_i * (theta - delta_im)
```

Where:
- b_i ~ difficulty_formula
- a_i ~ discrimination_formula
- xi_{k,m} ~ threshold_formula

---

## Usage Examples

### LPCM with both predictor types
```r
# Item data with both item-level and threshold-level covariates
item_data <- data.frame(
  abstractness = rnorm(20),      # Item-level covariate
  cognitive_level = rnorm(20)    # Threshold-level covariate
)

# Fit LPCM with both
fit_lpcm <- fit_eirt(
  responses,
  item_data = item_data,
  difficulty_formula = ~ abstractness,        # Item location
  threshold_formula = ~ cognitive_level,      # Threshold-specific effects
  model = "LPCM"
)

# Model: delta_im = b_i(abstractness) + xi_m(cognitive_level) + e_im
```

### GPCM with discrimination predictors
```r
# Item data
item_data <- data.frame(
  difficulty_pred = rnorm(20),
  discrimination_pred = rnorm(20)
)

# Fit GPCM with discrimination predictors
fit_gpcm <- fit_eirt(
  responses,
  item_data = item_data,
  difficulty_formula = ~ difficulty_pred,
  discrimination_formula = ~ discrimination_pred,  # Predicts log(a_i)
  model = "GPCM"
)

# Discrimination: log(a_i) = delta * discrimination_pred + epsilon_a
```

### GRM with discrimination predictors
```r
fit_grm <- fit_eirt(
  responses,
  item_data = item_data,
  difficulty_formula = ~ abstractness,
  discrimination_formula = ~ item_type,  # Predicts log(a_i)
  model = "GRM"
)
```

---

## Recommendation

**For your use case:**

1. **If you need threshold-level predictors**: Use **LPCM**
   - Has both item-level and threshold-level predictors
   - Discrimination fixed at 1 (Rasch family)

2. **If you need discrimination predictors**: Use **GPCM** or **GRM**
   - Have item-level difficulty and discrimination predictors
   - No threshold-level predictors (thresholds are free parameters)

3. **If you need BOTH threshold predictors AND discrimination predictors**:
   - ❌ Not currently implemented
   - Would require implementing LGPCM (Linear Generalized Partial Credit Model)
   - This would be a natural extension of the current framework

---

## Implementation Notes

From `src/gllamm_eirt.hpp`:

**Discrimination is always computed from predictors:**
```cpp
// Line 123-128
Type log_discrim_pred = 0.0;
for (int p = 0; p < p_disc; p++) {
  log_discrim_pred += delta(p) * W_discrimination(j, p);
}
discrimination(j) = exp(log_discrim_pred + epsilon_a(j));
```

**LPCM uses threshold predictors:**
```cpp
// Lines 287-291
Type delta_im = difficulty(item);  // b_i
for (int p = 0; p < p_thresh; p++) {
  delta_im += W_threshold(item, p) * xi(p, m - 1);  // threshold effects
}
delta_im += e_step(item, m - 1);  // residual
```

**But LPCM does NOT use discrimination:**
```cpp
// Line 293 - no "a *" multiplier
cumsum(m) = cumsum(m-1) + (theta(person) - delta_im);
```

Compare to **GPCM which DOES use discrimination:**
```cpp
// Line 265
cumsum(m) = cumsum(m-1) + a * (theta(person) - delta_im);
```

---

## Conclusion

✅ **Yes, LPCM has both item-level and threshold-level predictors**
✅ **Yes, GRM and GPCM have discrimination predictors**
❌ **No model currently combines threshold predictors WITH discrimination predictors**

The current implementation follows the standard IRT literature where:
- Rasch-family models (PCM, LPCM) have fixed discrimination = 1
- 2PL-family models (GPCM, GRM) have varying discrimination but fixed step structures
- LPCM is the only model with threshold-level regression

If you need a model with all three types of predictors (item-level difficulty, threshold-level, AND discrimination), that would be an LGPCM extension not currently implemented.
