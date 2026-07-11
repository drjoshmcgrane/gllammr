# Multilevel GLMMs with gllammr

## Introduction

Note: the code in this vignette is illustrative and is not evaluated
when the vignette is built.

This vignette covers multilevel generalized linear mixed models (GLMMs)
for hierarchical data.

## Two-Level Models

### Random Intercept

``` r

library(gllammr)

# Students in schools
fit1 <- gllamm(math_score ~ ses + (1 | school),
               data = school_data,
               family = gaussian())
summary(fit1)
```

### Random Slopes

``` r

# Random slope for SES
fit2 <- gllamm(math_score ~ ses + (ses | school),
               data = school_data)
summary(fit2)
```

## Three-Level Models

``` r

# Students in classes in schools
fit3 <- gllamm(achievement ~ time + (1 | school/class),
               data = multilevel_data)
```

## Non-Gaussian Responses

### Binomial

``` r

fit4 <- gllamm(passed ~ hours_studied + (1 | school),
               family = binomial())
```

### Poisson

``` r

fit5 <- gllamm(num_visits ~ treatment + (1 | clinic),
               family = poisson())
```

### Ordinal

``` r

fit6 <- fit_ordinal(satisfaction ~ service_quality + (1 | store),
                    data = customer_data)
```

## Model Diagnostics

``` r

# Standard diagnostic plots
plot(fit1)

# Goodness of fit
gof.gllamm(fit1)

# ICC
icc(fit1)
```

## Model Comparison

``` r

AIC(fit1, fit2)
BIC(fit1, fit2)
```
