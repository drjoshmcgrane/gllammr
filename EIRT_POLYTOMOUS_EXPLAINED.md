# How EIRT Works for Polytomous Data: Complete Explanation

**Date:** 7 Feb 2026

---

## Overview

Polytomous Explanatory IRT combines:
1. **Polytomous IRT models** (items with 3+ ordered categories)
2. **Explanatory modeling** (item parameters as functions of item characteristics)
3. **Multi-level structure** (item-level AND threshold-level predictors)

This creates a very flexible framework for understanding what makes polytomous items harder/easier and how response scales function.

---

## Part 1: Standard Polytomous IRT (No Covariates)

### Graded Response Model (GRM)

For an item with K categories (e.g., 1, 2, 3, 4 for a 4-point Likert scale):

**Cumulative probabilities:**
```
P(Y ≥ k | θ) = 1 / (1 + exp(-a(θ - b_k)))
```

Where:
- θ = person ability
- a = item discrimination (slope)
- b_k = threshold k (K-1 thresholds total)
- Constraint: b_1 < b_2 < ... < b_{K-1} (ordered)

**Category probabilities:**
```
P(Y = 1 | θ) = 1 - P(Y ≥ 2 | θ)
P(Y = k | θ) = P(Y ≥ k | θ) - P(Y ≥ k+1 | θ)  for k = 2, ..., K-1
P(Y = K | θ) = P(Y ≥ K | θ)
```

**Example:** 4-point scale (Strongly Disagree to Strongly Agree)
```
Category 1: P(Y=1) = P(Y < threshold_1)
Category 2: P(Y=2) = P(threshold_1 ≤ Y < threshold_2)
Category 3: P(Y=3) = P(threshold_2 ≤ Y < threshold_3)
Category 4: P(Y=4) = P(Y ≥ threshold_3)
```

### Parameters per Item

For each item j in standard GRM:
- a_j: discrimination (1 parameter)
- b_{j,1}, b_{j,2}, ..., b_{j,K-1}: thresholds (K-1 parameters)

**Total:** K parameters per item

For 20 items with 4 categories each: 20 × 4 = 80 item parameters!

---

## Part 2: Explanatory IRT - Item-Level Covariates Only

### The Problem with Standard IRT

With many items, you have many parameters:
- Hard to interpret patterns across items
- Need large samples
- Can't predict parameters for new items

### The EIRT Solution

**Model item parameters as functions of item characteristics:**

Instead of:
```
a_1 = 1.2, a_2 = 0.8, a_3 = 1.5, ...
b_{1,1} = -0.5, b_{2,1} = 0.3, b_{3,1} = -0.2, ...
```

We model:
```
log(a_j) = δ_0 + δ_1 × item_type_j + ε_{a,j}
b_{j,k} = γ_0 + γ_1 × word_frequency_j + ε_{b,j} + offset_k
```

### Mathematical Formulation

**Discrimination (item-level):**
```
log(a_j) = W_{disc,j}' δ + ε_{a,j}

where:
  W_{disc,j} = design matrix row for item j (item characteristics)
  δ = discrimination regression coefficients
  ε_{a,j} ~ N(0, σ²_a) = item-specific residual
```

**Difficulty (item-level):**
```
b_j = W_{diff,j}' γ + ε_{b,j}

where:
  W_{diff,j} = design matrix row for item j
  γ = difficulty regression coefficients
  ε_{b,j} ~ N(0, σ²_b) = item-specific residual
```

**Thresholds (basic):**
```
b_{j,1} = b_j + τ_1
b_{j,2} = b_j + τ_1 + exp(τ_2)
b_{j,3} = b_j + τ_1 + exp(τ_2) + exp(τ_3)
...
```

The `exp()` transformation ensures thresholds are ordered: b_{j,1} < b_{j,2} < ...

### What This Achieves

**Advantages:**
1. Reduces parameters: Instead of K × n_items, we have p_diff + p_disc coefficients
2. Interpretable: "Word frequency decreases difficulty by 0.3 units"
3. Can predict for new items: Given word_frequency, predict difficulty
4. Tests hypotheses: "Do abstract items have higher discrimination?"

**Example:**
```r
# Item characteristics
item_data <- data.frame(
  word_frequency = c(2.3, 1.1, 3.5, ...), # log frequency
  length = c(5, 8, 4, ...),                # word length
  abstractness = c(0, 1, 0, ...)           # 0=concrete, 1=abstract
)

# Fit EIRT
fit <- fit_eirt(
  responses,
  item_data,
  difficulty_formula = ~ word_frequency + length + abstractness,
  discrimination_formula = ~ abstractness,
  model = "GRM"
)

# Results
fit$regression_coefficients$difficulty
# γ_0 (Intercept):     0.52
# γ_1 (word_frequency): -0.31  # Common words easier
# γ_2 (length):         0.15   # Longer words harder
# γ_3 (abstractness):   0.68   # Abstract words harder

fit$regression_coefficients$discrimination
# δ_0 (Intercept):    0.82  # log(a) for concrete items
# δ_1 (abstractness): 0.25  # Abstract items have higher discrimination
```

**Interpretation:**
- A common word (high frequency) is easier (lower b)
- Longer words are harder (higher b)
- Abstract words are harder and discriminate better

---

## Part 3: Threshold-Level Covariates (NEW Enhancement)

### The Additional Problem

In the item-level EIRT above, ALL thresholds shift together:
- If an item is "hard", all its thresholds shift up
- But threshold SPACING is still fixed

**Real-world phenomenon:** Item characteristics might affect not just overall difficulty, but HOW people use the response scale.

### Examples Where Threshold Covariates Matter

**Example 1: Abstract Items**
```
Concrete item: "I like apples"
  Category 1-2 threshold: -1.5
  Category 2-3 threshold:  0.0
  Category 3-4 threshold:  1.5
  Spacing: 1.5 units between each

Abstract item: "I believe in justice"
  Category 1-2 threshold: -0.8
  Category 2-3 threshold:  0.0
  Category 3-4 threshold:  0.8
  Spacing: 0.8 units (COMPRESSED)
```

People might avoid extreme responses for abstract items (compressed scale), but use the full range for concrete items.

**Example 2: Positively vs Negatively Worded Items**
```
Positive item: "I am happy"
  - Might have expanded upper thresholds (easy to strongly agree)

Negative item: "I am sad"
  - Might have compressed upper thresholds (people avoid strong agreement)
```

### Mathematical Formulation with Threshold Covariates

**Without threshold covariates (item-level only):**
```
b_{j,1} = b_j + offset_1 + ε_{thr,j,1}
b_{j,2} = b_{j,1} + exp(offset_2 + ε_{thr,j,2})
b_{j,3} = b_{j,2} + exp(offset_3 + ε_{thr,j,3})

where:
  b_j = W_{diff,j}' γ + ε_{b,j}  (item-level difficulty)
  offset_k = fixed threshold positions
  ε_{thr,j,k} ~ N(0, σ²_thr) = threshold residuals
```

**With threshold covariates (NEW):**
```
b_{j,1} = b_j + W_{thr,j}' τ + ε_{thr,j,1}
b_{j,2} = b_{j,1} + exp(W_{thr,j}' τ + ε_{thr,j,2})
b_{j,3} = b_{j,2} + exp(W_{thr,j}' τ + ε_{thr,j,3})

where:
  W_{thr,j} = threshold covariate design matrix for item j
  τ = threshold regression coefficients
```

**Key insight:** W_{thr,j}' τ models threshold SPACING/POSITIONING as a function of item characteristics.

### Complete Model Structure

For item j with person i, category k:

**Person ability:**
```
θ_i ~ N(0, σ²_θ)
```

**Item discrimination:**
```
log(a_j) = W_{disc,j}' δ + ε_{a,j}
ε_{a,j} ~ N(0, σ²_a)
```

**Item difficulty (baseline for all thresholds):**
```
b_j = W_{diff,j}' γ + ε_{b,j}
ε_{b,j} ~ N(0, σ²_b)
```

**Thresholds (with covariates):**
```
b_{j,1} = b_j + W_{thr,j}' τ + ε_{thr,j,1}
b_{j,2} = b_{j,1} + exp(W_{thr,j}' τ + ε_{thr,j,2})
b_{j,3} = b_{j,2} + exp(W_{thr,j}' τ + ε_{thr,j,3})
...
ε_{thr,j,k} ~ N(0, σ²_thr)
```

**Response probability:**
```
P(Y_{ij} = k | θ_i) = P(Y_{ij} ≥ k | θ_i) - P(Y_{ij} ≥ k+1 | θ_i)

where:
P(Y_{ij} ≥ k | θ_i) = 1 / (1 + exp(-a_j(θ_i - b_{j,k})))
```

### Three-Level Decomposition

For each threshold b_{j,k}, we have:

```
b_{j,k} = [Item-level]  +  [Threshold-level]  +  [Residuals]
        = (W_{diff,j}' γ) + (W_{thr,j}' τ)      + (ε_{b,j} + ε_{thr,j,k})
          └─ Overall difficulty                  └─ Unexplained variation
                            └─ Threshold spacing/position
```

**Example:**
```r
item_data <- data.frame(
  word_frequency = 2.5,
  abstractness = 1      # Abstract item
)

# Item-level difficulty (affects ALL thresholds):
b_j = γ_0 + γ_1 × 2.5 = 0.5 + (-0.3) × 2.5 = -0.25

# Threshold-level spacing (affects how thresholds are positioned):
threshold_shift = τ_0 + τ_1 × 1 = 0.2 + (-0.4) × 1 = -0.2

# Final thresholds:
b_{j,1} = -0.25 + (-0.2) + ε_{thr,j,1} = -0.45 + residual
b_{j,2} = b_{j,1} + exp(-0.2 + ε_{thr,j,2})
b_{j,3} = b_{j,2} + exp(-0.2 + ε_{thr,j,3})
```

The negative τ_1 = -0.4 for abstractness means abstract items have **compressed thresholds** (smaller spacing).

---

## Part 4: Implementation in TMB

### Data Structures

**Input to TMB:**
```cpp
DATA_VECTOR(y);                    // Responses: [n_obs]
DATA_IVECTOR(person_id);           // Person index: [n_obs]
DATA_IVECTOR(item_id);             // Item index: [n_obs]
DATA_MATRIX(W_difficulty);         // Item covariates for difficulty: [n_items × p_diff]
DATA_MATRIX(W_discrimination);     // Item covariates for discrimination: [n_items × p_disc]
DATA_MATRIX(W_threshold);          // Item covariates for thresholds: [n_items × p_thr]
DATA_IVECTOR(n_categories_per_item); // Number of categories per item: [n_items]
DATA_INTEGER(threshold_covariate_model); // 0=no threshold covariates, 1=yes
```

**Parameters:**
```cpp
PARAMETER_VECTOR(theta);           // Person abilities: [n_persons]
PARAMETER_VECTOR(gamma);           // Difficulty coefficients: [p_diff]
PARAMETER_VECTOR(delta);           // Discrimination coefficients: [p_disc]
PARAMETER_VECTOR(tau);             // Threshold coefficients: [p_thr]
PARAMETER_VECTOR(epsilon_b);       // Item difficulty residuals: [n_items]
PARAMETER_VECTOR(epsilon_a);       // Item discrimination residuals: [n_items]
PARAMETER_MATRIX(threshold_resid); // Threshold residuals: [n_items × max_K-1]
PARAMETER(log_sigma_theta);        // Person ability SD
PARAMETER(log_sigma_epsilon_b);    // Difficulty residual SD
PARAMETER(log_sigma_epsilon_a);    // Discrimination residual SD
PARAMETER(log_sigma_threshold);    // Threshold residual SD
```

### Computation Flow

**Step 1: Compute item parameters from covariates**
```cpp
for (int j = 0; j < n_items; j++) {
  // Difficulty
  Type difficulty_pred = 0.0;
  for (int p = 0; p < p_diff; p++) {
    difficulty_pred += gamma(p) * W_difficulty(j, p);
  }
  difficulty(j) = difficulty_pred + epsilon_b(j);

  // Discrimination
  Type log_discrim_pred = 0.0;
  for (int p = 0; p < p_disc; p++) {
    log_discrim_pred += delta(p) * W_discrimination(j, p);
  }
  discrimination(j) = exp(log_discrim_pred + epsilon_a(j));
}
```

**Step 2: For each observation, compute thresholds**
```cpp
for (int i = 0; i < n_obs; i++) {
  int person = person_id(i);
  int item = item_id(i);
  int K = n_categories_per_item(item);

  // Build thresholds
  vector<Type> ordered_threshold(K - 1);

  if (threshold_covariate_model == 0) {
    // No threshold covariates
    ordered_threshold(0) = difficulty(item) + threshold_resid(item, 0);
    for (int k = 1; k < K - 1; k++) {
      ordered_threshold(k) = ordered_threshold(k-1) + exp(threshold_resid(item, k));
    }
  } else {
    // WITH threshold covariates
    for (int k = 0; k < K - 1; k++) {
      // Predicted threshold shift from covariates
      Type threshold_pred = 0.0;
      for (int p = 0; p < p_thr; p++) {
        threshold_pred += tau(p) * W_threshold(item, p);
      }

      if (k == 0) {
        // First threshold
        ordered_threshold(k) = difficulty(item) + threshold_pred + threshold_resid(item, k);
      } else {
        // Subsequent thresholds (maintain ordering)
        ordered_threshold(k) = ordered_threshold(k-1) +
                               exp(threshold_pred + threshold_resid(item, k));
      }
    }
  }
  // ... continue with likelihood computation
}
```

**Step 3: Compute category probabilities using GRM**
```cpp
int obs_cat = y(i) - 1;  // Convert to 0-indexed

Type prob_cat;
if (obs_cat == 0) {
  // Lowest category
  prob_cat = invlogit(discrimination(item) * (theta(person) - ordered_threshold(0)));
} else if (obs_cat == K - 1) {
  // Highest category
  prob_cat = 1.0 - invlogit(discrimination(item) *
                            (theta(person) - ordered_threshold(K - 2)));
} else {
  // Middle categories
  Type p_le_k = invlogit(discrimination(item) *
                         (theta(person) - ordered_threshold(obs_cat)));
  Type p_le_k_minus_1 = invlogit(discrimination(item) *
                                 (theta(person) - ordered_threshold(obs_cat - 1)));
  prob_cat = p_le_k - p_le_k_minus_1;
}
```

**Step 4: Add to negative log-likelihood**
```cpp
nll -= log(prob_cat + 1e-10);  // Avoid log(0)
```

### Random Effects

Which parameters are random (integrated out via Laplace approximation)?
```r
random = c("theta", "epsilon_b", "epsilon_a", "threshold_resid")
```

- `theta`: Person abilities (n_persons)
- `epsilon_b`: Item difficulty residuals (n_items)
- `epsilon_a`: Item discrimination residuals (n_items)
- `threshold_resid`: Item-threshold residuals (n_items × max_K-1)

**Total random effects:** n_persons + n_items + n_items + n_items × (max_K-1)

For 200 persons, 20 items, 4 categories:
- 200 + 20 + 20 + 20 × 3 = 300 random effects

These are integrated out using Laplace approximation (very fast in TMB).

---

## Part 5: Practical Example

### Setup

```r
# 15 items, 4-point Likert scale
item_data <- data.frame(
  item_id = 1:15,
  word_frequency = c(2.5, 1.8, 3.2, ...), # log frequency
  abstractness = c(0, 0, 1, 0, 1, ...)    # 0=concrete, 1=abstract
)

# 200 persons × 15 items
responses <- matrix(..., 200, 15)  # Values: 1, 2, 3, 4
```

### Fit Model

```r
fit <- fit_eirt(
  response_matrix = responses,
  item_data = item_data,

  # Item-level: Overall difficulty depends on word frequency
  difficulty_formula = ~ word_frequency,

  # Threshold-level: Abstract items have different threshold spacing
  threshold_formula = ~ abstractness,

  model = "GRM"
)
```

### Estimated Parameters

**Item-level difficulty:**
```
γ_0 (Intercept):      0.52   # Baseline difficulty for log(freq)=0
γ_1 (word_frequency): -0.31  # Each unit increase in log(freq) decreases difficulty
```

**Threshold-level:**
```
τ_0 (Intercept):     0.15   # Baseline threshold spacing for concrete items
τ_1 (abstractness):  -0.40  # Abstract items have compressed thresholds
```

**Variance components:**
```
σ_θ:         1.02  # Person ability variation
σ_epsilon_b: 0.25  # Item difficulty residual variation
σ_epsilon_a: 0.18  # Item discrimination residual variation
σ_threshold: 0.31  # Threshold residual variation
```

### Interpretation

**Item 3 (concrete, word_frequency = 3.2):**
```
b_3 = γ_0 + γ_1 × 3.2 + ε_{b,3}
    = 0.52 + (-0.31) × 3.2 + 0.12
    = -0.35

threshold_shift = τ_0 + τ_1 × 0 = 0.15

Thresholds:
  b_{3,1} = -0.35 + 0.15 = -0.20
  b_{3,2} = -0.20 + exp(0.15) = 0.96
  b_{3,3} = 0.96 + exp(0.15) = 2.12

Threshold spacing: 1.16, 1.16 (evenly spaced)
```

**Item 5 (abstract, word_frequency = 1.8):**
```
b_5 = γ_0 + γ_1 × 1.8 + ε_{b,5}
    = 0.52 + (-0.31) × 1.8 + (-0.08)
    = -0.12

threshold_shift = τ_0 + τ_1 × 1 = 0.15 + (-0.40) = -0.25

Thresholds:
  b_{5,1} = -0.12 + (-0.25) = -0.37
  b_{5,2} = -0.37 + exp(-0.25) = 0.41
  b_{5,3} = 0.41 + exp(-0.25) = 1.19

Threshold spacing: 0.78, 0.78 (COMPRESSED due to abstractness)
```

**Conclusion:** Abstract items have compressed threshold spacing, making middle categories more likely and extreme categories harder to endorse.

---

## Part 6: Comparison with Standard IRT

### Standard GRM (No Covariates)

**Parameters estimated:** 20 items × 4 params = 80 parameters
```
Item 1: a_1, b_{1,1}, b_{1,2}, b_{1,3}
Item 2: a_2, b_{2,1}, b_{2,2}, b_{2,3}
...
Item 20: a_20, b_{20,1}, b_{20,2}, b_{20,3}
```

**Limitations:**
- 80 parameters hard to interpret
- Can't test hypotheses about item characteristics
- Can't predict for new items
- Need large sample sizes

### EIRT with Item-Level Covariates

**Parameters estimated:** p_diff + p_disc + n_items × 2 residuals
```
Difficulty: γ_0, γ_1 (2 params)
Discrimination: δ_0 (1 param)
Residuals: ε_{b,1}, ..., ε_{b,20}, ε_{a,1}, ..., ε_{a,20} (40 random effects)
```

**Total:** 3 fixed + 40 random = 43 effective parameters (vs 80)

**Advantages:**
- Interpretable coefficients
- Can test hypotheses
- Can predict for new items
- More parsimonious

### EIRT with Item + Threshold Covariates

**Parameters estimated:** p_diff + p_disc + p_thr + residuals
```
Difficulty: γ_0, γ_1 (2 params)
Discrimination: δ_0 (1 param)
Threshold: τ_0, τ_1 (2 params)
Residuals: ε_b, ε_a, ε_thr (20 + 20 + 20×3 = 100 random effects)
```

**Total:** 5 fixed + 100 random = 105 effective parameters

**Additional advantages:**
- Models threshold structure
- Tests hypotheses about response scale usage
- Even more flexible

---

## Summary

### Three-Level Hierarchy

GLLAMMR's polytomous EIRT models item parameters at THREE levels:

1. **Item-level (difficulty, discrimination)**
   - `difficulty_formula = ~ word_frequency + complexity`
   - Affects overall item difficulty/discrimination
   - All thresholds shift together

2. **Threshold-level (NEW)**
   - `threshold_formula = ~ abstractness`
   - Affects threshold spacing and positioning
   - Models how response scales function differently

3. **Residual-level**
   - Item-specific deviations: ε_b, ε_a
   - Threshold-specific deviations: ε_thr
   - Captures unexplained variation

### Mathematical Model Summary

```
For item j, person i, category k:

P(Y_{ij} = k | θ_i) = GRM(θ_i, a_j, b_{j,1}, ..., b_{j,K-1})

where:
  θ_i ~ N(0, σ²_θ)

  log(a_j) = W_{disc,j}' δ + ε_{a,j}

  b_j = W_{diff,j}' γ + ε_{b,j}

  b_{j,1} = b_j + W_{thr,j}' τ + ε_{thr,j,1}
  b_{j,k} = b_{j,k-1} + exp(W_{thr,j}' τ + ε_{thr,j,k})  for k > 1
```

### Why This Is Powerful

1. **Reduces parameters** dramatically (80 → 5 fixed effects in example)
2. **Interpretable** - understand what makes items hard/easy
3. **Testable** - hypothesis tests about item characteristics
4. **Predictive** - predict parameters for new items
5. **Flexible** - models both item-level and threshold-level effects
6. **Unique** - no other software offers this flexibility!

---

**Created:** 7 Feb 2026
**For:** Understanding polytomous EIRT implementation in GLLAMMR
