# PCM vs LPCM: What's the Difference?

## The Key Distinction

Both PCM and LPCM have **item-level difficulty predictors** in our EIRT framework. The difference is in how they model the **step/threshold parameters**:

---

## PCM (Partial Credit Model)

### Model Structure
```
delta_im = b_i + s_im
```

Where:
- **b_i** = item location (explained by `difficulty_formula`)
- **s_im** = step deviations (FREE PARAMETERS, one set per item)

### Example: 3 Items, 4 Categories (3 thresholds each)

**Item 1:** s_{1,1} = -0.8, s_{1,2} = 0.3, s_{1,3} = 0.5
**Item 2:** s_{2,1} = -1.2, s_{2,2} = 0.8, s_{2,3} = 0.4
**Item 3:** s_{3,1} = -0.5, s_{3,2} = -0.2, s_{3,3} = 0.7

Each item has **completely independent** step parameters (9 parameters total).

### What's Being Estimated
- **Item location (b_i)**: Predicted from item covariates via `difficulty_formula`
- **Step deviations (s_im)**: Free parameters, **unique to each item**

**Total parameters:** n_items × (K-2) step parameters + parameters in difficulty_formula

---

## LPCM (Linear Partial Credit Model)

### Model Structure
```
delta_im = b_i + sum_k xi_{k,m} * x_{i,k} + e_{i,m}
```

Where:
- **b_i** = item location (explained by `difficulty_formula`)
- **xi_{k,m}** = threshold regression coefficients (explained by `threshold_formula`)
- **x_{i,k}** = item covariates
- **e_{i,m}** = threshold-specific residuals

### Example: Same 3 Items with Cognitive Level Covariate

Suppose `cognitive_level` is an item covariate:
- Item 1: cognitive_level = 0.5
- Item 2: cognitive_level = 1.2
- Item 3: cognitive_level = -0.3

**Threshold regression coefficients (shared across items):**
- xi_{cognitive_level, threshold1} = 0.6
- xi_{cognitive_level, threshold2} = 0.3
- xi_{cognitive_level, threshold3} = -0.2

**Resulting thresholds:**
- **Item 1:** δ_{1,1} = b_1 + 0.6×0.5 + e_{1,1}
- **Item 2:** δ_{2,1} = b_2 + 0.6×1.2 + e_{2,1}
- **Item 3:** δ_{3,1} = b_3 + 0.6×(-0.3) + e_{3,1}

Notice: All items use the **same coefficient** (0.6) for the effect of cognitive_level on threshold 1.

### What's Being Estimated
- **Item location (b_i)**: Predicted from item covariates via `difficulty_formula`
- **Threshold effects (xi_{k,m})**: Predicted from item covariates via `threshold_formula`, **shared across items**
- **Residuals (e_{i,m})**: Unexplained threshold variation

**Total parameters:** (K-1) × p_threshold + parameters in difficulty_formula + 1 residual SD

---

## The Fundamental Difference

### PCM: Unstructured Steps
- Each item has its own unique set of step parameters
- No assumption about WHY items have different steps
- More parameters, less structure

**Analogy:** Fixed effects for each item-step combination

### LPCM: Structured Steps
- Step parameters are **functions of item characteristics**
- Assumes item features explain threshold differences
- Fewer parameters, more structure
- Can make predictions for new items

**Analogy:** Random effects with covariates explaining variance

---

## Concrete Example

### Research Question: Do abstract words have more spread-out thresholds?

**Item data:**
```
Item 1: "cat"   - abstractness = 0.2
Item 2: "dog"   - abstractness = 0.3
Item 3: "love"  - abstractness = 0.9
Item 4: "truth" - abstractness = 0.95
```

### PCM Approach
```r
fit_pcm <- fit_eirt(
  responses,
  item_data = item_data,
  difficulty_formula = ~ abstractness,  # Item location varies with abstractness
  model = "PCM"
)
```

**What it estimates:**
- Effect of abstractness on item **location** (overall difficulty)
- Separate free step parameters for each item (unrelated to abstractness)

**Result:** You know abstract words are harder/easier, but no explanation for why their thresholds are spaced differently.

### LPCM Approach
```r
fit_lpcm <- fit_eirt(
  responses,
  item_data = item_data,
  difficulty_formula = ~ 1,              # Item location (no predictors)
  threshold_formula = ~ abstractness,    # Threshold spacing varies with abstractness
  model = "LPCM"
)
```

**What it estimates:**
- Effect of abstractness on threshold **spacing** (how spread out the categories are)
- Common regression coefficients applied to all items

**Result:** You can test whether abstract words have more/less compressed response scales.

### Best Approach (Use Both!)
```r
fit_lpcm_full <- fit_eirt(
  responses,
  item_data = item_data,
  difficulty_formula = ~ abstractness,   # Abstract words harder overall
  threshold_formula = ~ abstractness,     # Abstract words have different threshold spacing
  model = "LPCM"
)
```

**Result:** Separate effects of abstractness on overall difficulty AND threshold structure.

---

## When to Use Which?

### Use PCM when:
- ✅ You have no theory about threshold structure
- ✅ You want to estimate unique steps for each item
- ✅ You're doing exploratory analysis
- ✅ Sample size is large (can afford many parameters)

### Use LPCM when:
- ✅ You have item covariates that might explain threshold spacing
- ✅ You want to test hypotheses about threshold structure
- ✅ You want to predict thresholds for new items
- ✅ You want a more parsimonious model
- ✅ Sample size is limited

---

## Parameter Count Comparison

For **50 items, 5 categories (4 thresholds each)**:

### PCM
- Item locations: Predicted by difficulty_formula (e.g., 3 coefficients)
- Step deviations: 50 items × 3 free parameters = **150 parameters**
- **Total:** ~150 threshold-related parameters

### LPCM
- Item locations: Predicted by difficulty_formula (e.g., 3 coefficients)
- Threshold regression: 4 thresholds × p_threshold (e.g., 3) = **12 parameters**
- Residual variance: **1 parameter**
- **Total:** ~13 threshold-related parameters

**LPCM is MUCH more parsimonious** when you have many items!

---

## Can They Be Integrated?

### Conceptually: Yes!

You could think of PCM as LPCM with:
```r
# PCM is LPCM with dummy variables for each item-threshold:
threshold_formula = ~ item1:threshold1 + item1:threshold2 + ... + itemN:thresholdK
```

This would give each item-threshold its own parameter, equivalent to PCM.

### Practically: No!

They're kept separate because:
1. **Different purposes**: PCM for exploration, LPCM for explanation
2. **Different estimation**: PCM estimates free parameters, LPCM does regression
3. **Different interpretation**: PCM is descriptive, LPCM is explanatory
4. **Historical**: They come from different traditions (LPCM from Wilson's work on construct maps)

---

## Summary

| Aspect | PCM | LPCM |
|--------|-----|------|
| Item location | ✅ difficulty_formula | ✅ difficulty_formula |
| Threshold structure | ❌ Free parameters per item | ✅ threshold_formula (shared effects) |
| Parameters | Many (n_items × K) | Few (K × p_covariates) |
| Interpretation | Descriptive | Explanatory |
| New item prediction | ❌ No | ✅ Yes |
| Best for | Exploration | Hypothesis testing |

**Bottom line:** Both models have item-level difficulty predictors. The difference is that **LPCM adds structure to the threshold parameters** by explaining them with item covariates, while **PCM leaves them as free parameters**.

---

## Your Question Answered

> "If PCM already has item difficulty predictors, why do we need LPCM?"

**Answer:** PCM explains item **location** (overall difficulty) with predictors, but treats threshold **spacing** as free parameters unique to each item. LPCM explains BOTH location AND threshold spacing with predictors, providing a more structured and explanatory model.

Think of it like:
- **PCM**: Some items are harder than others (explained), but they all have different response patterns (unexplained)
- **LPCM**: Some items are harder than others (explained), AND their response patterns differ in systematic ways (also explained)

The key innovation of LPCM is explaining **within-item threshold structure** using item characteristics!
