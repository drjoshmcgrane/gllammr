# Multi-Level IRT Implementation Roadmap

**Date:** February 9, 2026

---

## Implementation Phases

### Phase 1: Single Additional Random Effect (PRIORITY 1)
**Duration:** 5-7 days
**Complexity:** Medium

**Goal:** Add one clustering level beyond person

**Example:**
```r
fit_irt(responses, person_data = data,
        random = ~ (1 | class), model = "2PL")
```

**Deliverables:**
- R interface with `person_data` and `random` parameters
- Random effects formula parser
- New TMB template for 2-level IRT
- Variance component estimation
- Basic tests and documentation

---

### Phase 2: Multiple Nested Levels (PRIORITY 1)
**Duration:** 3-5 days
**Complexity:** Medium-High

**Goal:** Handle arbitrary nesting depth

**Example:**
```r
fit_irt(responses, person_data = data,
        random = ~ (1 | state/district/school/class),
        model = "2PL")
```

**Deliverables:**
- Nested formula parsing
- Multiple variance components
- ICC calculations at each level
- Extended tests

---

### Phase 3: Crossed Random Effects (PRIORITY 2)
**Duration:** 4-6 days
**Complexity:** High

**Goal:** Handle non-nested structures

**Example:**
```r
fit_irt(responses, person_data = data,
        random = ~ (1 | student) + (1 | time),
        model = "Longitudinal 2PL")
```

**Deliverables:**
- Crossed effects detection
- Modified TMB template
- Longitudinal data handling
- Advanced tests

---

### Phase 4: EIRT Extension (PRIORITY 2)
**Duration:** 3-4 days
**Complexity:** Medium

**Goal:** Multi-level + explanatory items

**Example:**
```r
fit_eirt(responses, item_data, person_data,
         difficulty_formula = ~ word_freq,
         random = ~ (1 | school/class))
```

**Deliverables:**
- Extend EIRT with person_data
- Integration with existing EIRT
- Tests and examples

---

## Critical Design Decisions

### Decision 1: Parameter Structure

**Option A: Separate theta + group effects**
```cpp
theta_i = theta_0i + u_group[group[i]]
```
- Pro: Clear interpretation
- Pro: Easy to extract components
- Con: More parameters

**Option B: Composite theta only**
```cpp
theta_i = u_person[i] + u_group[group[i]]
```
- Pro: Fewer parameters
- Pro: More efficient
- Con: Less clear what's person vs group

**RECOMMENDATION:** Option A for interpretability

---

### Decision 2: Formula Interface

**Option A: lme4-style (RECOMMENDED)**
```r
random = ~ (1 | class) + (1 | school)
```
- Pro: Familiar to users
- Pro: Standard notation
- Con: Need to parse

**Option B: List-based**
```r
random = list(class = "class_id", school = "school_id")
```
- Pro: Simple
- Con: Non-standard
- Con: Limited expressiveness

**RECOMMENDATION:** Option A (lme4-style)

---

### Decision 3: Person Data Format

**Option A: person_data as data frame (RECOMMENDED)**
```r
person_data <- data.frame(
  person_id = 1:1000,
  class_id = ...,
  school_id = ...
)
fit_irt(responses, person_data = person_data,
        random = ~ (1 | class_id))
```
- Pro: Natural
- Pro: Can include other person-level covariates
- Con: Need to match rows to response_matrix

**Option B: Separate grouping vectors**
```r
fit_irt(responses, groups = list(class = class_vec, school = school_vec))
```
- Pro: Simple
- Con: Limited
- Con: No other person data

**RECOMMENDATION:** Option A with automatic row matching

---

### Decision 4: Missing Groups (Partial Nesting)

**Approach:** Use NA in grouping variable
```r
person_data$class_id <- c(1, 2, 2, NA, NA, 3, ...)
# Persons with NA have no class-level effect
```

**Implementation:**
- In TMB, check for NA (-1 in integer coding)
- Only add group effect if group[i] >= 0
- Natural handling of partial nesting

---

## Implementation Details

### Step 1: R Interface Changes

**File:** R/irt.R

**Add parameters:**
```r
fit_irt <- function(response_matrix,
                    model = c("Rasch", "2PL", "3PL", "GRM", "PCM", "GPCM", "NRM"),
                    person_data = NULL,  # NEW
                    random = NULL,       # NEW
                    weights = NULL,
                    mc_items = NULL,
                    start = NULL,
                    control = list())
```

**Process random effects:**
```r
if (!is.null(random)) {
  # Parse formula
  re_terms <- parse_random_formula(random, person_data)

  # Extract grouping factors
  group_factors <- lapply(re_terms, function(term) {
    factor(person_data[[term$group_var]])
  })

  # Pass to TMB
  tmb_data$has_random <- TRUE
  tmb_data$n_random_effects <- length(re_terms)
  tmb_data$group_ids <- do.call(cbind, lapply(group_factors, as.integer)) - 1L
  tmb_data$n_groups <- sapply(group_factors, nlevels)
} else {
  tmb_data$has_random <- FALSE
}
```

---

### Step 2: Formula Parsing

**Create:** R/parse_random.R

```r
#' Parse random effects formula
#' @keywords internal
parse_random_formula <- function(formula, data) {
  # Extract terms
  terms <- lme4::findbars(formula)

  if (is.null(terms)) {
    stop("No random effects found in formula")
  }

  # Parse each term
  result <- lapply(terms, function(term) {
    # term is like: (1 | class)
    # Extract grouping variable
    group_var <- as.character(term[[3]])

    # Check if in data
    if (!group_var %in% names(data)) {
      stop("Grouping variable '", group_var, "' not found in person_data")
    }

    list(
      type = "intercept",  # Only intercepts for now
      group_var = group_var,
      nesting = NULL       # Handle later
    )
  })

  # Check for nesting (e.g., school/class)
  # Expand nested notation
  result <- expand_nested_terms(result)

  return(result)
}

expand_nested_terms <- function(terms) {
  # Handle notation like (1 | school/class)
  # Convert to (1 | school) + (1 | school:class)
  # Implementation needed
}
```

---

### Step 3: TMB Template

**Create:** src/gllamm_irt_multilevel.hpp

```cpp
// Multi-Level IRT: Person-level + additional random effects
// Supports nested and crossed structures

template<class Type>
Type objective_function<Type>::operator() ()
{
  // ============================================================================
  // DATA INPUTS
  // ============================================================================

  // Standard IRT data
  DATA_VECTOR(y);
  DATA_IVECTOR(person_id);
  DATA_IVECTOR(item_id);
  DATA_INTEGER(n_persons);
  DATA_INTEGER(n_items);
  DATA_INTEGER(n_obs);
  DATA_VECTOR(weights);
  DATA_INTEGER(model_type);        // 1=Rasch, 2=2PL, 3=3PL
  DATA_IVECTOR(mc_items);

  // Multi-level structure
  DATA_INTEGER(has_random);        // 0=standard IRT, 1=multi-level
  DATA_INTEGER(n_random_effects);  // Number of RE levels
  DATA_IMATRIX(group_ids);         // [n_persons × n_random_effects]
  DATA_IVECTOR(n_groups);          // Number of groups per level

  // ============================================================================
  // PARAMETERS
  // ============================================================================

  // Item parameters (as before)
  PARAMETER_VECTOR(difficulty);
  PARAMETER_VECTOR(discrimination);
  PARAMETER_VECTOR(guessing);

  // Person-level deviations
  PARAMETER_VECTOR(theta_0);
  PARAMETER(log_sigma_theta);

  // Random effects (if has_random == 1)
  PARAMETER_MATRIX(u_random);      // [max_n_groups × n_random_effects]
  PARAMETER_VECTOR(log_sigma_random); // SD for each RE level

  // ============================================================================
  // INITIALIZE
  // ============================================================================

  Type nll = 0.0;
  Type sigma_theta = exp(log_sigma_theta);

  // ============================================================================
  // PRIORS
  // ============================================================================

  // Person-level deviations
  for (int p = 0; p < n_persons; p++) {
    nll -= dnorm(theta_0(p), Type(0.0), sigma_theta, true);
  }

  // Random effects priors
  if (has_random == 1) {
    for (int re = 0; re < n_random_effects; re++) {
      Type sigma_re = exp(log_sigma_random(re));
      int n_groups_re = n_groups(re);

      for (int g = 0; g < n_groups_re; g++) {
        nll -= dnorm(u_random(g, re), Type(0.0), sigma_re, true);
      }
    }
  }

  // ============================================================================
  // LIKELIHOOD
  // ============================================================================

  for (int i = 0; i < n_obs; i++) {
    int person = person_id(i);
    int item = item_id(i);

    // Compose ability: person + random effects
    Type theta = theta_0(person);

    if (has_random == 1) {
      for (int re = 0; re < n_random_effects; re++) {
        int group = group_ids(person, re);
        if (group >= 0) {  // -1 indicates NA
          theta += u_random(group, re);
        }
      }
    }

    // IRT likelihood (as before)
    Type prob;

    if (model_type == 1) {
      // Rasch
      Type eta = theta - difficulty(item);
      prob = invlogit(eta);

    } else if (model_type == 2) {
      // 2PL
      Type eta = discrimination(item) * (theta - difficulty(item));
      prob = invlogit(eta);

    } else {
      // 3PL with selective guessing
      Type eta = discrimination(item) * (theta - difficulty(item));

      if (mc_items(item) == 1) {
        Type c = guessing(item);
        prob = c + (Type(1.0) - c) * invlogit(eta);
      } else {
        prob = invlogit(eta);
      }
    }

    // Log-likelihood
    Type w_i = weights(i);
    nll -= w_i * (y(i) * log(prob + Type(1e-10)) +
                   (Type(1.0) - y(i)) * log(Type(1.0) - prob + Type(1e-10)));
  }

  // ============================================================================
  // REPORT
  // ============================================================================

  ADREPORT(theta_0);
  ADREPORT(difficulty);

  if (model_type >= 2) {
    ADREPORT(discrimination);
  }

  if (model_type == 3) {
    ADREPORT(guessing);
  }

  ADREPORT(sigma_theta);

  if (has_random == 1) {
    vector<Type> sigma_random(n_random_effects);
    for (int re = 0; re < n_random_effects; re++) {
      sigma_random(re) = exp(log_sigma_random(re));
    }
    ADREPORT(sigma_random);
    ADREPORT(u_random);
  }

  return nll;
}
```

---

### Step 4: Result Object

**Extend result structure:**
```r
result <- list(
  model = model,
  item_parameters = ...,
  person_abilities = theta_hat,
  ability_sd = sigma_theta_hat,

  # NEW: Multi-level components
  random_effects = if (has_random) {
    list(
      u_random = u_random_hat,      # Matrix of group effects
      sigma_random = sigma_random_hat, # Vector of SDs
      group_names = names(re_terms),
      n_groups = n_groups,
      icc = compute_icc(...)
    )
  } else NULL,

  ...
)

class(result) <- c("gllamm_irt_multilevel", "gllamm_irt", "gllamm")
```

---

### Step 5: S3 Methods

**VarCorr() method:**
```r
#' @export
VarCorr.gllamm_irt_multilevel <- function(x, ...) {
  if (is.null(x$random_effects)) {
    return(NULL)
  }

  result <- data.frame(
    Groups = c(x$random_effects$group_names, "Person", "Residual"),
    Variance = c(x$random_effects$sigma_random^2,
                 x$ability_sd^2,
                 1),  # Logistic variance
    Std.Dev. = c(x$random_effects$sigma_random,
                 x$ability_sd,
                 1)
  )

  class(result) <- c("VarCorr.gllamm", "data.frame")
  return(result)
}
```

**ICC() method:**
```r
#' @export
icc.gllamm_irt_multilevel <- function(x, level = NULL, ...) {
  if (is.null(x$random_effects)) {
    return(NULL)
  }

  # Total variance
  var_total <- sum(x$random_effects$sigma_random^2) +
               x$ability_sd^2 + pi^2/3  # Logistic

  if (is.null(level)) {
    # All levels
    iccs <- c(x$random_effects$sigma_random^2,
              x$ability_sd^2) / var_total
    names(iccs) <- c(x$random_effects$group_names, "Person")
    return(iccs)
  } else {
    # Specific level
    idx <- which(x$random_effects$group_names == level)
    if (length(idx) == 0) {
      stop("Level '", level, "' not found")
    }
    return(x$random_effects$sigma_random[idx]^2 / var_total)
  }
}
```

---

## Testing Strategy

### Unit Tests

**Test 1: Single level matches no clustering when sigma ≈ 0**
```r
test_that("Multi-level reduces to standard when sigma_group = 0", {
  # Fit both models
  # Check equivalence
})
```

**Test 2: Variance component recovery**
```r
test_that("Recovers known variance components", {
  # Simulate data with known sigmas
  # Fit and check recovery
})
```

**Test 3: ICC calculation**
```r
test_that("ICC computed correctly", {
  # Check ICC formula
  # Sum to 1
})
```

### Integration Tests

**Test 4: Real nested data**
```r
# Use educational dataset
# Schools > Classes > Students
# Check reasonable results
```

---

## Timeline

### Week 1 (Phase 1 - Part 1)
- Design finalization
- R interface implementation
- Formula parsing

### Week 2 (Phase 1 - Part 2)
- TMB template creation
- Testing and debugging
- Documentation

### Week 3 (Phase 2)
- Multiple nested levels
- Extended tests
- Vignette

### Week 4 (Phase 3-4)
- Crossed effects (if needed)
- EIRT extension
- Polish and release

---

## Questions for User

1. **Priority:** Is Phase 1 (single level) sufficient initially, or do you need multiple nested levels immediately?

2. **Interface:** Does the lme4-style `random = ~ (1 | group)` syntax work for you?

3. **Data format:** Is `person_data` as a data frame acceptable, or do you prefer a different format?

4. **Polytomous:** Should we start with dichotomous IRT only, or include polytomous from the start?

5. **EIRT:** Is multi-level EIRT needed immediately, or can it wait?

6. **Use cases:** What are your primary use cases? (helps prioritize features)

---

*Roadmap created: February 9, 2026*
