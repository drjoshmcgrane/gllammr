# De Boeck/Wilson vs GLLAMM: Relationship and Differences

## De Boeck & Wilson lme4 Approach

### 1. Rasch Model
```r
response ~ 0 + item + (1 | person)
```

**Equivalent to:**
- **Fixed effects:** Item difficulties (one per item)
- **Random effects:** Person abilities θ ~ N(0, σ²_θ)

**Model:** logit(P(Y_ij = 1)) = θ_j - b_i

### 2. LLTM (Linear Logistic Test Model)
```r
response ~ 0 + item_predictor + (1 | person)
```

**Equivalent to:**
- **Fixed effects:** Item difficulty predicted from covariates (W × γ)
- **Random effects:** Person abilities θ ~ N(0, σ²_θ)

**Model:**
- b_i = W_i' × γ (no residual!)
- logit(P(Y_ij = 1)) = θ_j - b_i

**Problem:** Assumes item covariates perfectly predict difficulty (no residual variance)

### 3. LLTM + Error
```r
response ~ 0 + item_predictor + (1 | person) + (1 | item)
```

**Equivalent to:**
- **Fixed effects:** Item difficulty predicted from covariates
- **Random effects:**
  - Person abilities: θ ~ N(0, σ²_θ)
  - Item residuals: ε_b ~ N(0, σ²_ε)

**Model:**
- b_i = W_i' × γ + ε_b[i]
- logit(P(Y_ij = 1)) = θ_j - b_i

**This is Explanatory IRT!**

---

## Our EIRT Implementation: Dichotomous

### Model Structure
```cpp
// From src/gllamm_eirt.hpp
difficulty(j) = difficulty_pred + epsilon_b(j);
// where difficulty_pred = W_difficulty * gamma

// Likelihood:
prob = invlogit(theta(person) - difficulty(item));
```

**In statistical notation:**
- b_i = W_i' × γ + ε_b[i], where ε_b ~ N(0, σ²_ε_b)
- logit(P(Y_ij = 1)) = θ_j - b_i
- θ ~ N(0, σ²_θ)

**This is EXACTLY De Boeck's LLTM + error!**

### Usage Comparison

**De Boeck (lme4):**
```r
# Must reshape data to long format first
data_long <- expand.grid(person = 1:100, item = 1:20)
data_long$response <- ...
data_long$word_freq <- item_data$word_freq[data_long$item]

# Fit
library(lme4)
fit <- glmer(response ~ 0 + word_freq + (1 | person) + (1 | item),
             data = data_long,
             family = binomial)
```

**GLLAMM (our package):**
```r
# Data in wide format (person × item matrix)
responses <- matrix(...)  # 100 × 20

# Item covariates
item_data <- data.frame(word_freq = ...)

# Fit
fit <- fit_eirt(responses,
                item_data = item_data,
                difficulty_formula = ~ word_freq,
                model = "Rasch")  # or "2PL"
```

**Same model, cleaner syntax!**

---

## Polytomous Models: The Complexity

### De Boeck & Wilson Approach

For polytomous models, De Boeck & Wilson resort to **data expansion**:

#### Data Expansion Example

Original data (4 categories: 0, 1, 2, 3):
```
Person 1, Item 1: Response = 2
```

Expanded data (3 thresholds):
```
Person 1, Item 1, Threshold 1: "Passed" = 1  (response ≥ 1)
Person 1, Item 1, Threshold 2: "Passed" = 1  (response ≥ 2)
Person 1, Item 1, Threshold 3: "Passed" = 0  (response ≥ 3)
```

**lme4 syntax:**
```r
# Expanded data
expanded_data <- data.frame(
  person = c(1, 1, 1, ...),
  item = c(1, 1, 1, ...),
  threshold = c(1, 2, 3, ...),
  passed = c(1, 1, 0, ...)
)

# Fit (example for GRM-like model)
fit <- glmer(passed ~ 0 + threshold + (1 | person) + (1 | item),
             data = expanded_data,
             family = binomial)
```

**Problem:** The 3 observations from same person-item are NOT independent!
- Response = 2 means: passed threshold 1 AND 2 BUT NOT 3
- These are correlated observations

### Rabe-Hesketh "Exploded Likelihood"

To handle the dependency, Rabe-Hesketh et al. use **exploded likelihood**:

Instead of treating each threshold observation as independent Bernoulli, they:
1. Expand the data to person-item-threshold format
2. Model threshold crossing probabilities
3. Use the **full categorical likelihood** that respects dependency

**Exploded likelihood:** Treats the expanded data as if from independent binary responses, but then corrects the likelihood to account for the constraint that all thresholds for same person-item are dependent.

---

## Our GLLAMM Implementation: Polytomous

### Key Difference: We DON'T Use Exploded Likelihood

**We use direct categorical likelihood in TMB.**

#### For GRM (lines 181-206 in gllamm_eirt.hpp):

```cpp
// We compute P(Y = k) directly for each category
if (obs_cat == 0) {
  prob_cat = Type(1.0) - invlogit(a * (theta(person) - ordered_threshold(0)));
} else if (obs_cat == K - 1) {
  prob_cat = invlogit(a * (theta(person) - ordered_threshold(K - 2)));
} else {
  Type p_ge_k = invlogit(a * (theta(person) - ordered_threshold(obs_cat - 1)));
  Type p_ge_k_plus_1 = invlogit(a * (theta(person) - ordered_threshold(obs_cat)));
  prob_cat = p_ge_k - p_ge_k_plus_1;
}

// Categorical log-likelihood
nll -= w_i * log(prob_cat + Type(1e-10));
```

**This is the proper categorical likelihood:**
- P(Y_ij = 0) = 1 - F(θ - τ₁)
- P(Y_ij = k) = F(θ - τ_k) - F(θ - τ_{k+1}) for k = 1, ..., K-2
- P(Y_ij = K-1) = F(θ - τ_{K-1})

**No data expansion required!**

#### For PCM/GPCM (lines 215-273):

```cpp
// Adjacent-categories formulation with direct categorical likelihood
vector<Type> cumsum(K);
cumsum(0) = Type(0.0);

for (int m = 1; m < K; m++) {
  Type delta_im = difficulty(item) + s_im;
  cumsum(m) = cumsum(m-1) + a * (theta(person) - delta_im);
}

Type denom = Type(0.0);
for (int m = 0; m < K; m++) {
  denom += exp(cumsum(m));
}
prob_cat = exp(cumsum(obs_cat)) / denom;
```

**This is the multinomial logit likelihood (softmax):**
- P(Y_ij = k) = exp(sum_{m=0}^k (θ - δ_im)) / sum_{m'=0}^{K-1} exp(sum_{m=0}^{m'} (θ - δ_im))

**No data expansion required!**

---

## LPCM: Our Implementation

### Mathematical Model

```
delta_im = b_i + sum_k xi_{k,m} * x_{i,k} + e_{i,m}
```

Where:
- **b_i**: Item location (from difficulty_formula)
- **xi_{k,m}**: Threshold regression coefficients (from threshold_formula)
- **e_{i,m}**: Threshold residuals ~ N(0, σ²_e)

### Code (lines 281-301):

```cpp
for (int m = 1; m < K; m++) {
  // Threshold difficulty = item location + threshold-specific regression
  Type delta_im = difficulty(item);  // b_i
  for (int p = 0; p < p_thresh; p++) {
    delta_im += W_threshold(item, p) * xi(p, m - 1);  // threshold covariates
  }
  delta_im += e_step(item, m - 1);  // threshold residual

  cumsum(m) = cumsum(m-1) + (theta(person) - delta_im);
}
```

### De Boeck Equivalent?

**There is no clean lme4 equivalent for LPCM!**

You would need to:
1. Expand data to person-item-threshold format
2. Model threshold-specific effects
3. Handle dependency via exploded likelihood

**Our approach is much cleaner:**
```r
fit_lpcm <- fit_eirt(
  responses,  # Wide format: person × item
  item_data = item_data,
  difficulty_formula = ~ abstractness,     # Item-level
  threshold_formula = ~ cognitive_level,   # Threshold-level
  model = "LPCM"
)
```

---

## Summary Table

| Aspect | De Boeck (lme4) | Rabe-Hesketh GLLAMM | Our GLLAMM |
|--------|----------------|---------------------|------------|
| **Dichotomous** | LLTM + error via glmer | Direct binary likelihood | Same as Rabe-Hesketh |
| **Data format** | Long (person-item rows) | Long | **Wide (matrix)** ✨ |
| **Polytomous approach** | Data expansion | Exploded likelihood | **Direct categorical** ✨ |
| **Complexity** | High (manual expansion) | Medium | **Low (automatic)** ✨ |
| **LPCM support** | ❌ No clean solution | ✅ Yes | ✅ Yes |
| **Threshold predictors** | ❌ Not practical | ✅ Yes | ✅ Yes |
| **Software** | lme4 + manual coding | Stata gllamm | **R + TMB** ✨ |

---

## Relationship Summary

### Dichotomous EIRT

**De Boeck LLTM + error:**
```r
response ~ 0 + item_predictor + (1 | person) + (1 | item)
```

**Our EIRT:**
```r
fit_eirt(responses, item_data,
         difficulty_formula = ~ item_predictor,
         model = "Rasch")
```

**✅ Same model, different syntax**

### Polytomous EIRT

**De Boeck approach:**
- Expand data → person-item-threshold format
- Fit binary model to expanded data
- Issues with dependency

**Rabe-Hesketh approach:**
- Expand data
- Use exploded likelihood to handle dependency
- Complex implementation

**Our approach:**
- **No data expansion** ✨
- Direct categorical likelihood in TMB
- Clean matrix interface
- Supports all polytomous models (GRM, PCM, GPCM, LPCM)

**✅ Same goals, better implementation**

---

## Advantages of Our Approach

### 1. No Manual Data Expansion
```r
# De Boeck: Manual expansion required
data_long <- expand.grid(person = 1:n, item = 1:J, threshold = 1:K)
# ... complex reshaping ...

# GLLAMM: Automatic
fit <- fit_eirt(responses, ...)  # Just pass matrix
```

### 2. Proper Categorical Likelihood
- De Boeck: Treat thresholds as independent → need exploded likelihood correction
- **GLLAMM: Direct categorical likelihood → no correction needed**

### 3. Cleaner Syntax for LPCM
```r
# GLLAMM
fit_lpcm <- fit_eirt(
  responses,
  item_data,
  difficulty_formula = ~ item_chars,
  threshold_formula = ~ threshold_chars,
  model = "LPCM"
)

# lme4: Would require complex manual setup with data expansion
```

### 4. Unified Framework
All models (dichotomous and polytomous) use same interface:
```r
fit_eirt(responses, item_data, difficulty_formula, model = "...")
```

---

## What We're Doing Differently

1. **Data structure:** Wide format (person × item matrix) instead of long format
2. **Likelihood:** Direct categorical likelihood instead of exploded likelihood
3. **Implementation:** TMB C++ templates instead of lme4 R code
4. **Interface:** Specialized IRT functions instead of general GLMM syntax

**Result:** Cleaner, faster, more intuitive for IRT modeling!

---

## Conclusion

**Conceptually:** Our EIRT models are equivalent to De Boeck's LLTM + error approach and consistent with Rabe-Hesketh's GLLAMM framework.

**Technically:** We implement them more efficiently:
- No manual data expansion
- Direct categorical likelihood (no exploded likelihood needed)
- Matrix interface (more natural for IRT)
- TMB automatic differentiation (faster optimization)

**For polytomous models:** We go beyond what's practical in lme4, providing clean implementations of GRM, PCM, GPCM, and especially LPCM with threshold-level predictors.

**Bottom line:** We're implementing the same statistical models as De Boeck and Rabe-Hesketh, but with a modern, efficient, IRT-focused interface! 🎉
