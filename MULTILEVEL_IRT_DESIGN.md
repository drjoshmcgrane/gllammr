# Multi-Level IRT Design

**Date:** February 9, 2026
**Status:** 🔄 **DESIGN PHASE**

---

## Overview

Add capability for multi-level random effects in IRT models to handle hierarchical/clustered data structures common in educational assessment.

---

## Use Cases

### Educational Assessment Examples

**1. Fully Nested Structure:**
```
Students > Classes > Schools > Districts > States
```

**2. Partially Nested:**
- Some students in classrooms, others not (e.g., home-schooled)
- Some schools in districts, others independent

**3. Crossed Random Effects:**
- Students × Time points (repeated measures)
- Students × Raters (multiple scorers)
- Items × Schools (DIF by school)

---

## Current Limitation

**Current IRT models:**
```r
fit_irt(response_matrix, model = "2PL")
```

**Random effects:**
- Only person abilities: `theta_i ~ N(0, sigma_theta²)`
- No clustering structure
- Assumes all persons independent

**Problem:**
- Ignores clustering (students in classes)
- Underestimates standard errors
- Can't model between-group variation
- Can't estimate ICC at different levels

---

## Proposed Design

### API Extension

**Standard IRT:**
```r
fit_irt(
  response_matrix,
  person_data = NULL,        # NEW: Data frame with person-level variables
  random = NULL,             # NEW: Random effects formula
  model = "2PL",
  ...
)
```

**EIRT:**
```r
fit_eirt(
  response_matrix,
  item_data,
  person_data = NULL,        # NEW: Person-level data
  difficulty_formula = ~ word_freq,
  random = NULL,             # NEW: Random effects formula
  model = "2PL",
  ...
)
```

---

## Formula Syntax

### Nested Structures

**Example 1: Students in classes**
```r
fit_irt(
  responses,
  person_data = student_data,
  random = ~ (1 | class),
  model = "2PL"
)
# theta_i = theta_0i + u_class[class[i]]
# theta_0i ~ N(0, sigma_theta²)
# u_class ~ N(0, sigma_class²)
```

**Example 2: Fully nested (students > classes > schools)**
```r
fit_irt(
  responses,
  person_data = student_data,
  random = ~ (1 | school/class),  # Nested notation
  model = "2PL"
)
# Equivalent to: ~ (1 | school) + (1 | school:class)
# theta_i = theta_0i + u_school[school[i]] + u_class[class[i]]
```

**Example 3: Explicit nesting**
```r
fit_irt(
  responses,
  person_data = student_data,
  random = ~ (1 | school) + (1 | class:school),
  model = "2PL"
)
```

**Example 4: Deep nesting**
```r
fit_irt(
  responses,
  person_data = student_data,
  random = ~ (1 | state/district/school/class),
  model = "2PL"
)
# theta_i = theta_0i + u_state + u_district + u_school + u_class
```

### Crossed Structures

**Example 5: Crossed random effects (students × time)**
```r
# Longitudinal IRT: same students, multiple time points
fit_irt(
  responses,
  person_data = long_data,  # Each row = person × time
  random = ~ (1 | student) + (1 | time),
  model = "2PL"
)
# theta_it = u_student[i] + u_time[t]
```

**Example 6: Rater effects**
```r
# Multiple raters scoring same students
fit_irt(
  responses,
  person_data = rating_data,
  random = ~ (1 | student) + (1 | rater),
  model = "2PL"
)
# theta_ir = u_student[i] + u_rater[r]
```

### Partially Nested

**Example 7: Some students in classes, others not**
```r
# person_data$class has NA for students not in classes
fit_irt(
  responses,
  person_data = student_data,
  random = ~ (1 | class),
  model = "2PL"
)
# theta_i = theta_0i + u_class[class[i]]  if class[i] not NA
# theta_i = theta_0i                       if class[i] is NA
```

---

## Mathematical Model

### Standard IRT (Current)

**2PL Model:**
```
P(Y_ij = 1 | theta_i) = logit^(-1)(a_j * (theta_i - b_j))
theta_i ~ N(0, sigma_theta²)
```

### Multi-Level IRT (Proposed)

**Nested Structure (Students in Classes in Schools):**
```
P(Y_ijk = 1 | theta_ijk) = logit^(-1)(a_j * (theta_ijk - b_j))
theta_ijk = theta_0ijk + u_school[school[k]] + u_class[class[k]]
theta_0ijk ~ N(0, sigma_theta²)
u_school ~ N(0, sigma_school²)
u_class ~ N(0, sigma_class²)
```

**Crossed Structure (Students × Time):**
```
P(Y_ijt = 1 | theta_ijt) = logit^(-1)(a_j * (theta_ijt - b_j))
theta_ijt = u_student[i] + u_time[t]
u_student ~ N(0, sigma_student²)
u_time ~ N(0, sigma_time²)
```

**General Formulation:**
```
theta_i = sum_g u_g[group_g[i]]
u_g ~ N(0, sigma_g²)
```

---

## Data Structure

### Person Data Format

**Nested Example:**
```r
person_data <- data.frame(
  person_id = 1:1000,
  class_id = rep(1:40, each = 25),      # 40 classes, 25 students each
  school_id = rep(1:10, each = 100),    # 10 schools, 4 classes each
  district_id = rep(1:2, each = 500)    # 2 districts, 5 schools each
)

# Response matrix: 1000 persons × 50 items
responses <- matrix(rbinom(50000, 1, 0.6), 1000, 50)

fit <- fit_irt(
  responses,
  person_data = person_data,
  random = ~ (1 | district_id/school_id/class_id),
  model = "2PL"
)
```

**Crossed Example:**
```r
# Longitudinal: 200 students × 3 time points = 600 observations
person_data <- expand.grid(
  student_id = 1:200,
  time_point = 1:3
)
person_data$person_id <- 1:600  # Unique identifier for each observation

# Response matrix: 600 observations × 30 items
responses <- matrix(rbinom(18000, 1, 0.6), 600, 30)

fit <- fit_irt(
  responses,
  person_data = person_data,
  random = ~ (1 | student_id) + (1 | time_point),
  model = "2PL"
)
```

**Partially Nested Example:**
```r
person_data <- data.frame(
  person_id = 1:1000,
  class_id = c(rep(1:30, each = 30), rep(NA, 100))  # 100 students not in classes
)

fit <- fit_irt(
  responses,
  person_data = person_data,
  random = ~ (1 | class_id),  # Only applies to 900 students
  model = "2PL"
)
```

---

## Implementation Plan

### Phase 1: Single Random Effect (Beyond Person)

**Goal:** Add one additional random effect level

**Example:**
```r
fit_irt(responses, person_data, random = ~ (1 | class))
```

**Changes:**
1. Add `person_data` and `random` parameters to `fit_irt()`
2. Parse random effects formula
3. Create grouping factor mapping
4. Modify TMB template to handle 2-level structure
5. Estimate two variance components

### Phase 2: Multiple Nested Random Effects

**Goal:** Handle arbitrary nesting depth

**Example:**
```r
fit_irt(responses, person_data, random = ~ (1 | school/class))
```

**Changes:**
1. Parse nested formula notation
2. Handle multiple random effects simultaneously
3. Estimate multiple variance components
4. Compute ICCs at each level

### Phase 3: Crossed Random Effects

**Goal:** Handle non-nested structures

**Example:**
```r
fit_irt(responses, person_data, random = ~ (1 | student) + (1 | time))
```

**Changes:**
1. Detect crossed vs nested structure
2. Modify TMB to handle crossed effects
3. Different parameterization

### Phase 4: EIRT Extension

**Goal:** Multi-level + explanatory item parameters

**Example:**
```r
fit_eirt(responses, item_data, person_data,
         difficulty_formula = ~ word_freq,
         random = ~ (1 | school/class))
```

---

## TMB Template Changes

### Current Structure (Person-Level Only)

```cpp
// Current: gllamm_irt.hpp
PARAMETER_VECTOR(theta);     // Person abilities
PARAMETER(log_sigma_theta);  // SD of abilities

// Prior
for (int p = 0; p < n_persons; p++) {
  nll -= dnorm(theta(p), Type(0.0), sigma_theta, true);
}

// Likelihood
Type eta = discrimination(item) * (theta(person) - difficulty(item));
```

### Proposed: Single Additional Level

```cpp
// New: gllamm_irt_multilevel.hpp
DATA_IVECTOR(group_id);         // Group identifier for each person
DATA_INTEGER(n_groups);         // Number of groups

PARAMETER_VECTOR(theta_0);      // Person-level deviations
PARAMETER_VECTOR(u_group);      // Group-level random effects
PARAMETER(log_sigma_theta);     // Person-level SD
PARAMETER(log_sigma_group);     // Group-level SD

// Priors
for (int p = 0; p < n_persons; p++) {
  nll -= dnorm(theta_0(p), Type(0.0), sigma_theta, true);
}
for (int g = 0; g < n_groups; g++) {
  nll -= dnorm(u_group(g), Type(0.0), sigma_group, true);
}

// Compose ability
for (int i = 0; i < n_obs; i++) {
  int person = person_id(i);
  int group = group_id(person);

  Type theta_total = theta_0(person) + u_group(group);

  Type eta = discrimination(item) * (theta_total - difficulty(item));
  // ... likelihood
}
```

### Proposed: Multiple Nested Levels

```cpp
DATA_INTEGER(n_random_effects);           // Number of RE levels
DATA_IMATRIX(group_ids);                  // [n_persons × n_random_effects]
DATA_IVECTOR(n_groups_per_level);         // Number of groups per level

PARAMETER_VECTOR(theta_0);                // Person-level deviations
PARAMETER_VECTOR(u_random_combined);      // All random effects concatenated
PARAMETER_VECTOR(log_sigma_random);       // SD for each RE level

// Build theta with multiple levels
for (int i = 0; i < n_obs; i++) {
  int person = person_id(i);

  Type theta_total = theta_0(person);

  // Add each random effect level
  for (int re = 0; re < n_random_effects; re++) {
    int group = group_ids(person, re);
    if (group >= 0) {  // -1 indicates NA (partial nesting)
      theta_total += u_random_combined(group);
    }
  }

  // Likelihood using theta_total
}
```

---

## Output and Methods

### Summary Output

```r
> summary(fit)
Multi-Level IRT Model (2PL)

Random Effects:
                Variance    SD      ICC
(Person)        0.850      0.922   0.425
(Class)         0.450      0.671   0.225
(School)        0.700      0.837   0.350
Total           2.000      1.414   1.000

Number of groups:
  Classes:  40
  Schools:  10

Fixed Effects (Item Parameters):
  ... (as before)
```

### Methods

**VarCorr() - Variance Components:**
```r
VarCorr(fit)
#   Groups   Name  Variance Std.Dev.
#   School   (Int) 0.700    0.837
#   Class    (Int) 0.450    0.671
#   Person   (Int) 0.850    0.922
```

**ICC() - Intraclass Correlations:**
```r
icc(fit, level = "school")  # School-level ICC
icc(fit, level = "class")   # Class-level ICC
icc(fit)                    # All levels
```

**ranef() - Random Effects:**
```r
ranef(fit, level = "school")  # School effects
ranef(fit, level = "class")   # Class effects
ranef(fit, level = "person")  # Person effects (theta)
```

---

## Validation and Testing

### Test Cases

**Test 1: Verify single-level reproduces current results**
```r
# Without clustering
fit1 <- fit_irt(responses, model = "2PL")

# With clustering but no variance
fit2 <- fit_irt(responses, person_data,
                random = ~ (1 | class), model = "2PL")

# If sigma_class ≈ 0, estimates should match
```

**Test 2: Simulated nested data**
```r
# Generate data with known structure
sigma_school <- 0.8
sigma_class <- 0.5
sigma_person <- 1.0

# Simulate and recover parameters
```

**Test 3: Crossed effects**
```r
# Longitudinal data
# Check that student and time effects are separated
```

**Test 4: Partial nesting**
```r
# Some students without class assignment
# Check that NA handling works
```

---

## Benefits

### 1. Realistic Modeling
- Accounts for clustering in data
- Correct standard errors
- Valid inference

### 2. Substantive Questions
- How much variation is between vs within schools?
- Are there school-level effects on ability?
- ICCs at different levels

### 3. Efficiency
- Shares information across groups
- Better estimates for small groups
- Handles unbalanced designs

### 4. Flexibility
- Fully nested
- Crossed
- Partially nested
- Any combination

---

## Related Features

### Future Extensions

**1. Random Slopes (not just intercepts):**
```r
random = ~ (1 + time | student)  # Growth models
```

**2. Group-Level Predictors:**
```r
random = ~ (1 | school),
school_formula = ~ mean_ses + urbanicity
```

**3. Item-Level Random Effects (Random DIF):**
```r
random = ~ (1 | person) + (1 | school:item)  # Item × school interaction
```

**4. Spatial Random Effects:**
```r
random = ~ (1 | school) + spatial(lat, lon)
```

---

## Priority

**High Priority:**
- Single additional random effect (Phase 1)
- Nested structures (Phase 2)
- Standard IRT models first

**Medium Priority:**
- Crossed random effects (Phase 3)
- EIRT extension (Phase 4)

**Low Priority:**
- Random slopes
- Group-level predictors
- Complex spatial structures

---

## Next Steps

1. **User feedback** on API design
2. **Prototype** single-level implementation
3. **Test** with simulated data
4. **Extend** to multiple levels
5. **Document** with examples
6. **Vignette** on multi-level IRT

---

*Design document created: February 9, 2026*
