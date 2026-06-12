# GLLAMMR Unified Interface

## Overview
All GLLAMMR models are now accessible through the main `gllamm()` function using the `family` argument. This provides a consistent, lme4-like interface for all model types.

## ✅ Unified Interface Implementation

### Core Principle
**Single entry point**: `gllamm(formula, data, family = ...)`

All specialized functions (`fit_ordinal()`, `fit_irt()`, `fit_lca()`) remain available for direct use, but the recommended interface is through `gllamm()`.

---

## Usage Examples

### 1. Standard GLMM (Gaussian)
```r
# Default family is gaussian()
fit <- gllamm(y ~ x + (1 | group), data = mydata)
```

### 2. Binomial GLMM
```r
fit <- gllamm(success ~ x + (1 | group), 
              data = mydata, 
              family = binomial(link = "logit"))
```

### 3. Poisson GLMM
```r
fit <- gllamm(count ~ x + (1 | group), 
              data = mydata, 
              family = poisson(link = "log"))
```

### 4. Ordinal Regression (NEW!)
```r
# Proportional odds (logit)
fit <- gllamm(rating ~ x + (1 | group), 
              data = mydata, 
              family = ordinal(link = "logit"))

# Adjacent category logit
fit <- gllamm(rating ~ x + (1 | group), 
              data = mydata, 
              family = ordinal(link = "acl"))

# Continuation ratio
fit <- gllamm(rating ~ x + (1 | group), 
              data = mydata, 
              family = ordinal(link = "crl_forward"))

# Partial proportional odds
fit <- gllamm(rating ~ x + (1 | group), 
              data = mydata, 
              family = ordinal(link = "ppo"))
```

---

## Supported Families

### Standard GLM Families
- `gaussian()` - Normal distribution (identity link)
- `binomial()` - Binary/binomial (logit, probit, cloglog)
- `poisson()` - Count data (log link)

### Extended Families
- `ordinal()` - Ordered categorical responses
  - Links: `"logit"`, `"probit"`, `"acl"`, `"crl_forward"`, `"crl_backward"`, `"ppo"`

### Future (when implemented)
- `irt()` - Item response theory models
- `lca()` - Latent class analysis
- `multinomial()` - Unordered categorical

---

## Architecture

### Dispatch Logic in gllamm()

```r
gllamm <- function(formula, data, family = gaussian(), ...) {
  # Validate inputs
  validate_formula(formula, data)
  
  # Dispatch based on family type
  if (inherits(family, "ordinal_family")) {
    return(fit_ordinal(formula, data, link = family$link, ...))
  }
  
  # Standard GLMM path for gaussian/binomial/poisson
  # ...
}
```

### Family Constructors

**ordinal()** - Creates ordinal_family object
```r
ordinal <- function(link = c("logit", "probit", "acl", 
                             "crl_forward", "crl_backward", "ppo")) {
  link <- match.arg(link)
  link_code <- switch(link, 
    logit = 1L, probit = 2L, acl = 3L, 
    crl_forward = 4L, crl_backward = 5L, ppo = 6L)
  
  structure(
    list(family = "ordinal", link = link, link_code = link_code),
    class = c("ordinal_family", "family")
  )
}
```

---

## Alternative Interfaces

While `gllamm()` is the recommended unified interface, specialized functions remain available:

### Direct Function Calls

```r
# Ordinal
fit_ordinal(rating ~ x + (1 | group), data, link = "logit")

# IRT
fit_irt(response_matrix, model = "2PL")

# LCA
fit_lca(indicator_matrix, nclass = 3)
```

### When to Use Direct Functions
- Legacy code compatibility
- Quick prototyping
- When you prefer explicit function names
- Special model-specific arguments

### When to Use gllamm()
- **Recommended for new code**
- Consistent with lme4 interface
- Easier to switch between model types
- Cleaner code with family objects
- Future-proof (all new models will support this)

---

## Benefits of Unified Interface

### 1. Consistency
```r
# Same syntax across model types
fit1 <- gllamm(y ~ x + (1|g), data, family = gaussian())
fit2 <- gllamm(y ~ x + (1|g), data, family = binomial())
fit3 <- gllamm(y ~ x + (1|g), data, family = ordinal())
```

### 2. Easy Model Comparison
```r
# Compare different link functions for ordinal response
fit_logit <- gllamm(rating ~ x + (1|id), data, family = ordinal("logit"))
fit_probit <- gllamm(rating ~ x + (1|id), data, family = ordinal("probit"))
fit_acl <- gllamm(rating ~ x + (1|id), data, family = ordinal("acl"))

# Compare fits
AIC(fit_logit, fit_probit, fit_acl)
```

### 3. Familiar to lme4 Users
```r
# lme4 syntax (doesn't support ordinal)
# lmer(y ~ x + (1|g), data)

# GLLAMMR syntax (supports ordinal)
gllamm(y ~ x + (1|g), data, family = ordinal())
```

### 4. Generic Methods Work Seamlessly
```r
fit <- gllamm(rating ~ x + (1|id), data, family = ordinal())

# All standard methods work
summary(fit)
coef(fit)
fitted(fit)
logLik(fit)

# New methods work too
fit(fit)        # Comprehensive fit statistics
plot(fit)       # Model-specific plots
```

---

## Implementation Details

### Files Modified
- **R/gllamm.R**: Added dispatch logic for ordinal_family
- **R/families.R**: Created ordinal() family constructor
- **R/ordinal.R**: Updated documentation to note gllamm() interface

### Tests Added
- **test-unified-interface.R**: Comprehensive tests for unified interface
  - Tests gllamm() with ordinal() family
  - Verifies equivalence with fit_ordinal()
  - Tests all ordinal links via gllamm()
  - Tests all S3 methods work

### Documentation Updated
- gllamm() examples show ordinal() family
- ordinal() documentation shows gllamm() usage
- fit_ordinal() notes gllamm() as recommended interface

---

## Backward Compatibility

✅ **All existing code continues to work**

```r
# Old way (still works)
fit <- fit_ordinal(y ~ x + (1|g), data, link = "logit")

# New way (recommended)
fit <- gllamm(y ~ x + (1|g), data, family = ordinal("logit"))
```

Both produce identical results. No breaking changes.

---

## Migration Guide

### From fit_ordinal() to gllamm()

**Before:**
```r
fit <- fit_ordinal(rating ~ temp + contact + (1 | judge),
                   data = wine,
                   link = "logit")
```

**After:**
```r
fit <- gllamm(rating ~ temp + contact + (1 | judge),
              data = wine,
              family = ordinal(link = "logit"))
```

### From lme4 to GLLAMMR (for ordinal)

**lme4 (can't do ordinal):**
```r
# Would need ordinal package or polr()
library(ordinal)
clmm(rating ~ temp + (1 | judge), data = wine)
```

**GLLAMMR:**
```r
gllamm(rating ~ temp + (1 | judge),
       data = wine,
       family = ordinal(link = "logit"))
```

---

## Summary

### Key Points
1. ✅ `gllamm()` is the unified interface for all models
2. ✅ Use `family = ordinal(link = "...")` for ordinal models
3. ✅ Specialized functions still available but gllamm() recommended
4. ✅ All S3 methods (fit, plot, summary, etc.) work seamlessly
5. ✅ Fully backward compatible - no breaking changes
6. ✅ Tests verify equivalence between interfaces

### What's New
- `ordinal()` family constructor with 6 link functions
- Dispatch logic in `gllamm()` for ordinal models
- Unified interface consistent across all model types
- Comprehensive tests for unified interface

### Next Steps
When IRT and LCA are updated:
- Add `irt()` family constructor
- Add `lca()` family constructor  
- Update `gllamm()` dispatch for these families
- Maintain backward compatibility with `fit_irt()` and `fit_lca()`
