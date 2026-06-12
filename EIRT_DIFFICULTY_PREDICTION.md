# How Item Difficulty is Predicted in EIRT

## The Model

In our EIRT implementation, item difficulty follows this structure:

```
b_i = (W_difficulty[i,] %*% gamma) + epsilon_b[i]
```

Where:
- **W_difficulty**: Design matrix from `difficulty_formula` applied to `item_data`
- **gamma**: Regression coefficients (estimated)
- **epsilon_b[i]**: Item-specific residual ~ N(0, sigma_epsilon_b²)

**Code reference:** `src/gllamm_eirt.hpp` lines 117-121

---

## Key Point: You Use ITEM COVARIATES, Not Item Dummies

### ❌ **NOT** like this (item as categorical predictor):
```r
# WRONG INTERPRETATION
difficulty_formula = ~ factor(item_id)  # This is NOT how it works!
```

### ✅ **Actually** like this (item characteristics):
```r
# CORRECT - Use item-level covariates
item_data <- data.frame(
  word_frequency = c(5.2, 3.1, 7.8, ...),  # Log frequency
  word_length = c(4, 7, 5, ...),            # Number of letters
  abstractness = c(0.3, 0.8, 0.2, ...)      # Abstractness rating
)

fit <- fit_eirt(responses,
                item_data = item_data,
                difficulty_formula = ~ word_frequency + word_length,
                model = "PCM")
```

**Model:** b_i = γ₀ + γ₁×word_frequency[i] + γ₂×word_length[i] + ε_b[i]

---

## Examples

### Example 1: Predict Difficulty from Word Frequency

```r
item_data <- data.frame(
  item_id = 1:20,
  word_freq = c(5.2, 3.1, 7.8, ...)  # Log frequency from corpus
)

fit_pcm <- fit_eirt(
  responses,
  item_data = item_data,
  difficulty_formula = ~ word_freq,  # More frequent = easier
  model = "PCM"
)

# Resulting model:
# b_i = gamma_0 + gamma_1 * word_freq[i] + epsilon_b[i]
```

**Interpretation:**
- `gamma_0`: Average item difficulty
- `gamma_1`: Effect of word frequency on difficulty (likely negative - frequent words are easier)
- `epsilon_b[i]`: Item-specific deviations not explained by frequency

### Example 2: Multiple Predictors

```r
item_data <- data.frame(
  item_id = 1:20,
  word_freq = ...,
  word_length = ...,
  abstractness = ...
)

fit_pcm <- fit_eirt(
  responses,
  item_data = item_data,
  difficulty_formula = ~ word_freq + word_length + abstractness,
  model = "PCM"
)

# Model: b_i = gamma_0 + gamma_1*freq + gamma_2*length + gamma_3*abstract + epsilon_b[i]
```

### Example 3: Intercept-Only (Closest to "Standard" IRT)

```r
fit_pcm <- fit_eirt(
  responses,
  item_data = item_data,
  difficulty_formula = ~ 1,  # Just intercept
  model = "PCM"
)

# Model: b_i = gamma_0 + epsilon_b[i]
```

**This is closest to "standard" PCM** where each item has its own difficulty, but:
- `gamma_0` is the grand mean difficulty
- `epsilon_b[i]` is the item-specific deviation

**Effectively:** Each item gets its own difficulty through `epsilon_b[i]`, which is like a random effect.

---

## Standard IRT vs EIRT

### Standard PCM (in packages like mirt, TAM)
```
b_i = free parameter (one per item)
```
- 20 items = 20 difficulty parameters
- No explanation of WHY items differ in difficulty

### PCM in EIRT (our implementation)
```
b_i = X_i' * gamma + epsilon_b[i]
```
- Regression coefficients (gamma) explain difficulty
- Residuals (epsilon_b) capture unexplained variation
- Can test hypotheses about what makes items difficult

---

## What Gets Estimated?

### With `difficulty_formula = ~ word_freq + length`

**Fixed effects (gamma):**
- γ₀: Intercept (baseline difficulty)
- γ₁: Effect of word frequency
- γ₂: Effect of word length

**Random effects (epsilon_b):**
- 20 values, one per item (deviations from prediction)

**Variance component:**
- σ²_epsilon_b: Residual variance in difficulty not explained by covariates

**Total difficulty parameters:**
- 3 regression coefficients (γ₀, γ₁, γ₂)
- 1 residual variance (σ²_epsilon_b)
- = 4 parameters (instead of 20 free difficulties)

---

## Comparison Table

| Specification | Model | Parameters | Interpretation |
|--------------|-------|-----------|----------------|
| `~ 1` | b_i = γ₀ + ε_b[i] | 2 (intercept + SD) | "Free" difficulties via random effects |
| `~ word_freq` | b_i = γ₀ + γ₁×freq + ε_b[i] | 3 | Frequency explains difficulty |
| `~ freq + length` | b_i = γ₀ + γ₁×freq + γ₂×len + ε_b[i] | 4 | Multiple predictors |

In all cases, epsilon_b[i] gives each item its own difficulty, but covariates **explain** some of that variation.

---

## How This Differs from "Item as Predictor"

### If You Did This (hypothetically):
```r
item_data$item_id <- factor(1:20)
difficulty_formula = ~ item_id  # DON'T DO THIS in EIRT!
```

**What would happen:**
- W_difficulty would be a 20×20 identity matrix (dummy coding)
- gamma would be 20 coefficients (one per item)
- This is equivalent to free parameters

**But this defeats the purpose of EIRT!** The whole point is to explain item difficulties with **item characteristics**, not just estimate separate parameters for each item.

---

## Key Insight

In EIRT, you're doing **item-level regression**, not person-level regression:

```
Person-level regression (standard GLMM):
  Y_ij ~ person characteristics (age, gender, etc.)

Item-level regression (EIRT):
  b_i ~ item characteristics (frequency, length, abstractness, etc.)
```

---

## Summary

**Your question:** "Is the prediction of item difficulty for PCM just in terms of 'item' being in the model?"

**Answer:** ❌ **No!**

In our EIRT implementation:
- Item difficulty is predicted from **item covariates** (word frequency, length, etc.)
- NOT from item ID as a categorical variable
- You provide item-level data with substantive predictors
- The model estimates regression coefficients showing how those predictors affect difficulty

**Example:**
```r
# You provide:
item_data <- data.frame(
  word = c("cat", "dog", "philosophy", ...),
  frequency = c(8.2, 7.9, 4.3, ...),      # Real item characteristics
  length = c(3, 3, 10, ...)
)

# Model estimates:
# b_i = gamma_0 + gamma_1 * frequency + gamma_2 * length + epsilon_b[i]

# NOT:
# b_i = gamma_1 * I(item==1) + gamma_2 * I(item==2) + ...  # This is NOT what we do
```

The **explanatory** in "Explanatory IRT" means we **explain item parameters with item covariates**, not just estimate separate parameters per item!
