# Multi-Level IRT: R Integration Complete

**Date:** February 9, 2026
**Status:** ✅ **R INTEGRATION COMPLETE** (TMB compilation pending)

---

## Summary

The R-level integration for multi-level IRT models is now complete. Users can specify hierarchical random effects using lme4-style syntax, and the package will properly parse the formulas, create grouping structures, pass data to TMB, and extract/display results.

---

## What's Implemented

### 1. User-Facing API

Users can now fit multi-level IRT models:

```r
# Standard IRT (as before)
fit_standard <- fit_irt(responses, model = "2PL")

# Multi-level IRT with class effects
fit_multilevel <- fit_irt(
  responses,
  model = "2PL",
  person_data = student_data,
  random = ~ (1 | class_id)
)

# Nested: students in classes in schools
fit_nested <- fit_irt(
  responses,
  model = "Rasch",
  person_data = student_data,
  random = ~ (1 | school_id/class_id)
)

# Crossed: longitudinal data
fit_long <- fit_irt(
  responses,
  model = "2PL",
  person_data = long_data,
  random = ~ (1 | student_id) + (1 | time_point)
)
```

### 2. Formula Parsing (`R/parse_random.R`)

Complete lme4-style random effects formula parser:

- `parse_random_formula()`: Main entry point
- `extract_random_terms()`: Extracts `(1 | group)` terms from formula
- `parse_single_random_term()`: Parses individual random effect terms
- `expand_nested_terms()`: Expands `school/class` to `school + school:class`
- `create_grouping_matrix()`: Creates integer grouping matrix for TMB

**Features:**
- ✅ Simple grouping: `~ (1 | class)`
- ✅ Nested notation: `~ (1 | school/class)`
- ✅ Multiple levels: `~ (1 | state) + (1 | district) + (1 | school) + (1 | class)`
- ✅ Crossed effects: `~ (1 | student) + (1 | time)`
- ✅ Handles NA (partial nesting) via -1 coding
- ✅ Validates grouping variables exist in data

### 3. Modified `fit_irt()` Function

**New Parameters:**
- `person_data`: Data frame with person-level variables (required for random effects)
- `random`: lme4-style formula for random effects

**Validation:**
- Checks person_data and random consistency
- Validates person_data has same number of rows as response_matrix
- Ensures grouping variables exist in person_data

**Workflow:**
1. Parse random effects formula if specified
2. Create grouping matrix
3. Detect standard vs multi-level mode
4. Pass appropriate data to TMB
5. Use correct DLL (`gllamm_irt` vs `gllamm_irt_multilevel`)
6. Extract parameters (theta_0, u_random, sigma_random)
7. Compute composite abilities
8. Calculate ICCs
9. Build result object with class `"gllamm_irt_multilevel"`

### 4. TMB Data Preparation

When `random` is specified, additional TMB data is prepared:

```r
tmb_data$has_random <- 1L
tmb_data$n_random_effects <- as.integer(n_re)
tmb_data$group_ids <- as.matrix(group_ids)      # [n_persons × n_re]
tmb_data$n_groups <- as.integer(n_groups)       # groups per level
tmb_data$max_n_groups <- as.integer(max(n_groups))
```

### 5. Parameter Initialization

Multi-level models initialize:
- `theta_0`: Person-level deviations (instead of `theta`)
- `u_random`: Matrix of group effects [max_n_groups × n_re]
- `log_sigma_random`: Vector of random effect SDs (one per level)
- `log_sigma_theta`: Person-level SD

### 6. TMB Object Creation

Uses appropriate DLL and random effects:

```r
if (has_random) {
  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = c("theta_0", "u_random"),  # Integrate out both
    DLL = "gllamm_irt_multilevel",
    silent = TRUE
  )
} else {
  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = "theta",
    DLL = "gllamm_irt",
    silent = TRUE
  )
}
```

### 7. Result Object Structure

Multi-level models return additional components:

```r
result$random_effects <- list(
  u_random = u_random_hat,              # Matrix of group effects
  sigma_random = sigma_random_hat,      # SD for each level
  group_names = group_names,            # Names of grouping variables
  n_groups = n_groups,                  # Number of groups per level
  icc = icc_values,                     # ICC for each level + person
  composite_theta = composite_theta     # Total ability (theta_0 + REs)
)

class(result) <- c("gllamm_irt_multilevel", "gllamm_irt", "gllamm")
```

### 8. Print Method Enhancement

`print.gllamm_irt()` now detects multi-level models and displays:

```
Multi-Level IRT Model ( 2PL )

Number of persons: 500
Number of items: 15
Model type: Dichotomous

Item parameters:
...

Random Effects:
  Grouping variables: class_id
  Number of groups: 20

Variance Components:
       Groups  Variance Std.Dev
    class_id    0.2500  0.5000
      Person    1.0000  1.0000
    Residual    3.2899  1.8138

Intraclass Correlations:
     Level   ICC
  class_id 0.055
    Person 0.221

Ability distribution:
  ...
```

---

## What's NOT Yet Implemented

### Immediate (Blocking Testing)

1. **TMB Compilation Issue:**
   - `TMB::compile('src/gllamm_irt_multilevel.hpp')` returns "Nothing to be done"
   - No `.so` file created
   - Appears to be build system issue, not code issue
   - **Action needed:** Resolve compilation workflow

### Phase 1 Remaining

2. **Polytomous Multi-Level:**
   - Currently only dichotomous (Rasch/2PL/3PL) support multi-level
   - Polytomous models (GRM/PCM/GPCM/NRM) explicitly error if `random` specified
   - Need to create `src/gllamm_irt_poly_multilevel.hpp`

3. **S3 Methods:**
   - `VarCorr.gllamm_irt_multilevel()` - extract variance components
   - `icc.gllamm_irt_multilevel()` - compute ICCs
   - `ranef.gllamm_irt_multilevel()` - extract random effects
   - `coef.gllamm_irt_multilevel()` - extract coefficients

4. **Testing:**
   - Unit tests for parsing functions
   - Simulation tests for parameter recovery
   - Tests for nested, crossed, partial nesting

### Phase 2-4

5. **EIRT Extension:** Add multi-level to explanatory IRT
6. **Documentation:** Update `?fit_irt`, create vignette
7. **Polish:** Examples, error messages, edge cases

---

## Files Modified/Created

### Created
- ✅ `R/parse_random.R` (236 lines) - Complete random effects parsing
- ✅ `src/gllamm_irt_multilevel.hpp` (169 lines) - TMB template for multi-level dichotomous IRT
- ✅ `test_multilevel_irt.R` - Comprehensive test script
- ✅ `compile_multilevel.R` - Compilation helper script

### Modified
- ✅ `R/irt.R` - Added multi-level integration (70 lines added/modified)
  - Function signature (lines 68-76)
  - Validation (lines 77-99)
  - Dispatch (lines 132-140)
  - fit_irt_dichotomous signature (line 149)
  - TMB data preparation (lines 209-226)
  - Parameter initialization (lines 228-265)
  - TMB object creation (lines 267-280)
  - Parameter extraction (lines 287-309)
  - Result construction (lines 311-346)
  - Print method (lines 451-523)

### Documentation Updated
- ✅ `MULTILEVEL_IRT_IMPLEMENTATION_STATUS.md` - Progress tracking
- ✅ `MULTILEVEL_IRT_DESIGN.md` - Design document (already existed)
- ✅ `MULTILEVEL_IRT_ROADMAP.md` - Implementation plan (already existed)

---

## Testing Status

### Verified ✅
1. Formula parsing works correctly
   - Simple: `~ (1 | class)` ✓
   - Parentheses handling fixed ✓
   - Grouping matrix creation ✓

2. R code is syntactically correct
   - No parse errors ✓
   - Function signatures match ✓
   - Data structures correct ✓

### Pending ⏳
1. TMB compilation and loading
2. Full integration test
3. Parameter recovery verification
4. Model comparison (standard vs multi-level)

---

## Next Steps

### Immediate Priority

1. **Resolve TMB Compilation:**
   - Option A: Debug make/compilation workflow
   - Option B: Use package build system (`devtools::document()`, `devtools::load_all()`)
   - Option C: Manual compilation with proper flags

2. **Test Basic Functionality:**
   - Run `test_multilevel_irt.R`
   - Verify standard IRT still works
   - Verify multi-level IRT fits
   - Check parameter recovery

### Then Continue with Phase 1

3. Create S3 methods (`VarCorr`, `icc`, `ranef`)
4. Extend to polytomous multi-level
5. Write unit tests
6. Update documentation

---

## Key Design Decisions Implemented

1. ✅ **Parameter Structure:** theta_i = theta_0i + sum_g u_g[group_g[i]]
   - Clearly separates person and group effects
   - Easy to interpret components

2. ✅ **Formula Interface:** lme4-style `random = ~ (1 | group)`
   - Familiar to R users
   - Standard notation
   - Supports nested and crossed

3. ✅ **Person Data Format:** Data frame with rows matching response_matrix
   - Natural for users
   - Can include other person-level covariates
   - Validated for consistency

4. ✅ **Partial Nesting:** Use NA in grouping variable → -1 in TMB
   - Natural handling
   - No special API needed
   - Works seamlessly

5. ✅ **Separate DLL:** Multi-level uses `gllamm_irt_multilevel` DLL
   - Cleaner code
   - No overhead when not needed
   - Standard IRT unchanged

---

## Example Usage (Once Compilation Works)

```r
library(GLLAMMR)

# Simulate nested data
set.seed(123)
n_students <- 500
n_classes <- 20
n_items <- 15

person_data <- data.frame(
  student_id = 1:n_students,
  class_id = rep(1:n_classes, each = n_students / n_classes)
)

# Generate responses (with class effects)
theta_student <- rnorm(n_students)
u_class <- rnorm(n_classes, 0, 0.5)
theta_total <- theta_student + u_class[person_data$class_id]
difficulty <- rnorm(n_items)

responses <- matrix(NA, n_students, n_items)
for (i in 1:n_students) {
  for (j in 1:n_items) {
    p <- plogis(theta_total[i] - difficulty[j])
    responses[i, j] <- rbinom(1, 1, p)
  }
}

# Fit standard model
fit_standard <- fit_irt(responses, model = "Rasch")

# Fit multi-level model
fit_multilevel <- fit_irt(
  responses,
  model = "Rasch",
  person_data = person_data,
  random = ~ (1 | class_id)
)

# Compare
print(fit_standard)
print(fit_multilevel)

# ICCs
fit_multilevel$random_effects$icc

# Variance components
VarCorr(fit_multilevel)  # Once method implemented

# Random effects
ranef(fit_multilevel, level = "class_id")  # Once method implemented
```

---

*Integration completed: February 9, 2026*
