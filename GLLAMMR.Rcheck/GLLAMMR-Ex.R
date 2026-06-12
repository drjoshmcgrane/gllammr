pkgname <- "GLLAMMR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "GLLAMMR-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('GLLAMMR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("VarCorr.gllamm_irt_multilevel")
### * VarCorr.gllamm_irt_multilevel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: VarCorr.gllamm_irt_multilevel
### Title: Extract Variance Components from Multi-Level IRT Models
### Aliases: VarCorr.gllamm_irt_multilevel

### ** Examples

## Not run: 
##D # Fit multi-level model
##D fit <- fit_irt(responses, model = "2PL",
##D                person_data = data, random = ~ (1 | class))
##D 
##D # Extract variance components
##D VarCorr(fit)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("VarCorr.gllamm_irt_multilevel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("abilities")
### * abilities

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: abilities
### Title: Extract Person Abilities from IRT Models
### Aliases: abilities

### ** Examples

## Not run: 
##D # Person-level deviations
##D abilities(fit)
##D 
##D # Total abilities (including class effects)
##D abilities(fit, composite = TRUE)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("abilities", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("binomial")
### * binomial

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: binomial
### Title: Binomial Family for Binary and Binomial Outcomes
### Aliases: binomial

### ** Examples

## Not run: 
##D # Logistic regression (default)
##D family1 <- binomial()
##D family2 <- binomial(link = "logit")
##D 
##D # Probit regression
##D family3 <- binomial(link = "probit")
##D 
##D # Complementary log-log for rare events
##D family4 <- binomial(link = "cloglog")
##D 
##D # Use with gllamm() - recommended interface
##D fit <- gllamm(outcome ~ age + treatment + (1 | clinic),
##D               data = mydata,
##D               family = binomial(link = "logit"))
##D 
##D # Rare event with cloglog
##D fit_rare <- gllamm(rare_disease ~ exposure + (1 | region),
##D                    data = epi_data,
##D                    family = binomial(link = "cloglog"))
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("binomial", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("coef.gllamm_irt")
### * coef.gllamm_irt

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: coef.gllamm_irt
### Title: Extract Coefficients from IRT Models
### Aliases: coef.gllamm_irt

### ** Examples

## Not run: 
##D # Item parameters
##D coef(fit, type = "item")
##D 
##D # Person abilities
##D coef(fit, type = "person")
##D 
##D # Random effects (multi-level only)
##D coef(fit, type = "random")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("coef.gllamm_irt", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compare_eirt")
### * compare_eirt

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compare_eirt
### Title: Compare Explanatory IRT Models
### Aliases: compare_eirt

### ** Examples

## Not run: 
##D # Fit model with predictor
##D fit1 <- fit_eirt(responses, item_data,
##D                  difficulty_formula = ~ word_freq,
##D                  model = "Rasch")
##D 
##D # Fit model without predictor
##D fit0 <- fit_eirt(responses, item_data,
##D                  difficulty_formula = ~ 1,
##D                  model = "Rasch")
##D 
##D # Compare models
##D compare_eirt(fit0, fit1)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compare_eirt", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dif_test")
### * dif_test

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dif_test
### Title: Test for Differential Item Functioning (DIF)
### Aliases: dif_test

### ** Examples

## Not run: 
##D # Fit IRT model
##D fit <- fit_irt(responses, model = "2PL")
##D 
##D # Test for DIF across gender
##D dif_result <- dif_test(fit, group = gender)
##D print(dif_result)
##D 
##D # Plot ICCs for flagged items
##D dif_plot(dif_result, item = 5)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dif_test", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dif_test_with_data")
### * dif_test_with_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dif_test_with_data
### Title: Test for DIF with explicit response data
### Aliases: dif_test_with_data

### ** Examples

## Not run: 
##D dif_result <- dif_test_with_data(responses, group = gender, model = "2PL")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dif_test_with_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit")
### * fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit
### Title: Generic Fit Statistics Function
### Aliases: fit

### ** Examples

## Not run: 
##D # GLMM
##D fit1 <- gllamm(y ~ x + (1 | group), data = data)
##D fit(fit1)
##D 
##D # IRT
##D fit2 <- fit_irt(responses, model = "2PL")
##D fit(fit2, compute_item_fit = TRUE)
##D 
##D # LCA
##D fit3 <- fit_lca(indicators, nclass = 3)
##D fit(fit3)
##D 
##D # Ordinal
##D fit4 <- fit_ordinal(rating ~ x + (1 | id), data = data)
##D fit(fit4, test_po = TRUE)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_binomial")
### * fit_binomial

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_binomial
### Title: Fit Binomial Regression Models with Random Effects
### Aliases: fit_binomial

### ** Examples

## Not run: 
##D # Simulate binary data
##D set.seed(123)
##D n_groups <- 20
##D n_per_group <- 10
##D data <- data.frame(
##D   group = rep(1:n_groups, each = n_per_group),
##D   x = rnorm(n_groups * n_per_group),
##D   y = rbinom(n_groups * n_per_group, 1, 0.5)
##D )
##D 
##D # Recommended: Use gllamm() with binomial() family
##D fit1 <- gllamm(y ~ x + (1 | group),
##D                data = data,
##D                family = binomial(link = "logit"))
##D summary(fit1)
##D 
##D # Probit link
##D fit2 <- gllamm(y ~ x + (1 | group),
##D                data = data,
##D                family = binomial(link = "probit"))
##D 
##D # Complementary log-log for rare events
##D data$rare_event <- rbinom(nrow(data), 1, 0.05)
##D fit3 <- gllamm(rare_event ~ x + (1 | group),
##D                data = data,
##D                family = binomial(link = "cloglog"))
##D summary(fit3)
##D 
##D # Alternative: Call fit_binomial() directly
##D fit4 <- fit_binomial(y ~ x + (1 | group),
##D                      data = data,
##D                      link = "logit")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_binomial", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_eirt")
### * fit_eirt

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_eirt
### Title: Fit Explanatory Item Response Theory Models
### Aliases: fit_eirt

### ** Examples

## Not run: 
##D # Dichotomous EIRT with item and discrimination predictors
##D item_chars <- data.frame(
##D   word_frequency = rnorm(20),
##D   item_length = rpois(20, 5),
##D   item_type = factor(sample(c("concrete", "abstract"), 20, replace = TRUE))
##D )
##D fit_2pl <- fit_eirt(responses, item_data = item_chars,
##D                     difficulty_formula = ~ word_frequency + item_length,
##D                     discrimination_formula = ~ item_type,
##D                     model = "2PL")
##D 
##D # Pure LLTM (no residuals)
##D fit_lltm <- fit_eirt(responses, item_data = item_chars,
##D                      difficulty_formula = ~ word_frequency,
##D                      model = "Rasch",
##D                      item_residuals = FALSE)
##D 
##D # Polytomous PCM (adjacent-categories logit, Rasch family)
##D fit_pcm <- fit_eirt(poly_responses, item_data = item_chars,
##D                     difficulty_formula = ~ abstractness,
##D                     model = "PCM")
##D 
##D # PCM with threshold predictors (LPCM framework)
##D fit_pcm_thresh <- fit_eirt(poly_responses, item_data = item_chars,
##D                            difficulty_formula = ~ abstractness,
##D                            threshold_formula = ~ cognitive_level,
##D                            model = "PCM")
##D 
##D # GPCM with all predictors
##D fit_gpcm <- fit_eirt(poly_responses, item_data = item_chars,
##D                      difficulty_formula = ~ word_frequency,
##D                      discrimination_formula = ~ item_type,
##D                      threshold_formula = ~ cognitive_level,
##D                      model = "GPCM")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_eirt", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_irt")
### * fit_irt

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_irt
### Title: Fit Item Response Theory Models
### Aliases: fit_irt

### ** Examples

## Not run: 
##D # Dichotomous example (Rasch)
##D set.seed(123)
##D n_persons <- 500
##D n_items <- 20
##D theta <- rnorm(n_persons, 0, 1)
##D difficulty <- rnorm(n_items, 0, 1)
##D 
##D # Generate binary responses
##D responses <- matrix(NA, n_persons, n_items)
##D for (i in 1:n_persons) {
##D   for (j in 1:n_items) {
##D     p <- plogis(theta[i] - difficulty[j])
##D     responses[i, j] <- rbinom(1, 1, p)
##D   }
##D }
##D 
##D # Fit Rasch model
##D fit_rasch <- fit_irt(responses, model = "Rasch")
##D summary(fit_rasch)
##D 
##D # Polytomous example (GRM)
##D # Generate 5-category responses
##D responses_poly <- matrix(NA, n_persons, n_items)
##D thresholds <- matrix(seq(-2, 2, length.out = 4), n_items, 4, byrow = TRUE)
##D for (i in 1:n_persons) {
##D   for (j in 1:n_items) {
##D     probs <- c(plogis(theta[i] - thresholds[j, 1]),
##D                diff(plogis(theta[i] - thresholds[j, ])),
##D                1 - plogis(theta[i] - thresholds[j, 4]))
##D     responses_poly[i, j] <- sample(1:5, 1, prob = probs)
##D   }
##D }
##D 
##D # Fit GRM model
##D fit_grm <- fit_irt(responses_poly, model = "GRM")
##D summary(fit_grm)
##D 
##D # 3PL with selective guessing (mixed MC and non-MC items)
##D # Assessment: 20 items, first 15 are MC, last 5 are open-ended
##D fit_3pl <- fit_irt(responses, model = "3PL", mc_items = 1:15)
##D # Only items 1-15 get guessing parameters
##D # Items 16-20 use 2PL likelihood (no guessing)
##D 
##D # Multi-level IRT: students nested in classes
##D person_data <- data.frame(
##D   person_id = 1:n_persons,
##D   class_id = rep(1:10, each = 50)
##D )
##D fit_multilevel <- fit_irt(responses, model = "2PL",
##D                            person_data = person_data,
##D                            random = ~ (1 | class_id))
##D # theta_i = theta_0i + u_class[class[i]]
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_irt", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_lca")
### * fit_lca

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_lca
### Title: Fit Latent Class Analysis Models
### Aliases: fit_lca

### ** Examples

## Not run: 
##D # Simulate 2-class data
##D set.seed(123)
##D n <- 500
##D 
##D # Class 1: high probability of yes
##D class1_probs <- c(0.8, 0.7, 0.9, 0.75)
##D # Class 2: low probability of yes
##D class2_probs <- c(0.2, 0.3, 0.1, 0.25)
##D 
##D # Generate data
##D true_class <- sample(1:2, n, replace = TRUE, prob = c(0.6, 0.4))
##D data <- matrix(NA, n, 4)
##D for (i in 1:n) {
##D   probs <- if (true_class[i] == 1) class1_probs else class2_probs
##D   data[i, ] <- rbinom(4, 1, probs)
##D }
##D colnames(data) <- paste0("Item", 1:4)
##D 
##D # Fit 2-class model
##D fit <- fit_lca(data, nclass = 2)
##D summary(fit)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_lca", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_mixed")
### * fit_mixed

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_mixed
### Title: Fit Joint Models for Mixed Response Types
### Aliases: fit_mixed

### ** Examples

## Not run: 
##D fit <- fit_mixed(
##D   formulas = list(gaussian = biomarker ~ age,
##D                   binomial = event ~ age + treatment),
##D   random = ~ (1 | patient),
##D   data = d)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_mixed", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_multinomial")
### * fit_multinomial

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_multinomial
### Title: Fit Multinomial Regression Models with Random Effects
### Aliases: fit_multinomial

### ** Examples

## Not run: 
##D # Simulate multinomial data
##D data$choice <- factor(sample(c("A", "B", "C"), 100, replace = TRUE))
##D 
##D # Fit multinomial model
##D fit <- fit_multinomial(choice ~ price + quality + (1 | person),
##D                        data = data)
##D summary(fit)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_multinomial", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_npml")
### * fit_npml

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_npml
### Title: Fit Two-Level GLMMs by Nonparametric Maximum Likelihood (NPML)
### Aliases: fit_npml

### ** Examples

## Not run: 
##D fit <- fit_npml(y ~ x + (1 | group), data = d, k = 3,
##D                 family = stats::binomial())
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_npml", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_ordinal")
### * fit_ordinal

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_ordinal
### Title: Fit Ordinal Regression Models with Random Effects
### Aliases: fit_ordinal

### ** Examples

## Not run: 
##D # Simulate ordinal data
##D data$satisfaction <- factor(sample(1:5, 100, replace = TRUE),
##D                             ordered = TRUE,
##D                             levels = 1:5)
##D 
##D # Recommended: Use gllamm() with ordinal() family
##D fit1 <- gllamm(satisfaction ~ age + (1 | clinic),
##D                data = data,
##D                family = ordinal(link = "logit"))
##D summary(fit1)
##D 
##D # Alternative: Call fit_ordinal() directly
##D fit2 <- fit_ordinal(satisfaction ~ age + (1 | clinic),
##D                     data = data,
##D                     link = "logit")
##D summary(fit2)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_ordinal", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_rank")
### * fit_rank

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_rank
### Title: Fit Rank-Ordered Logit Models with Taste Heterogeneity
### Aliases: fit_rank

### ** Examples

## Not run: 
##D fit <- fit_rank(rank ~ price + quality, case = ~ subject,
##D                 random = ~ (0 + price | region), data = d)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_rank", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_sem")
### * fit_sem

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_sem
### Title: Fit Structural Equation Models with Latent Variables
### Aliases: fit_sem

### ** Examples

## Not run: 
##D fit <- fit_sem(
##D   measurement = list(f1 = ~ x1 + x2 + x3, f2 = ~ y1 + y2 + y3),
##D   structural = list(f2 ~ f1),
##D   data = d)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_sem", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fit_survival")
### * fit_survival

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fit_survival
### Title: Fit Parametric Survival Models with Random Effects (Frailty)
### Aliases: fit_survival

### ** Examples

## Not run: 
##D fit <- fit_survival(Surv(time, status) ~ age + (1 | center),
##D                     data = d, distribution = "weibull")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fit_survival", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gllamm")
### * gllamm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gllamm
### Title: Fit Generalized Linear Latent and Mixed Models
### Aliases: gllamm

### ** Examples

## Not run: 
##D # Basic random intercept model
##D data(sleepstudy, package = "lme4")
##D fit1 <- gllamm(Reaction ~ Days + (1 | Subject),
##D                data = sleepstudy)
##D summary(fit1)
##D 
##D # Random intercept and slope
##D fit2 <- gllamm(Reaction ~ Days + (Days | Subject),
##D                data = sleepstudy)
##D summary(fit2)
##D 
##D # Extract components
##D fixef(fit2)        # Fixed effects
##D ranef(fit2)        # Random effects
##D VarCorr(fit2)      # Variance components
##D fitted(fit2)       # Fitted values
##D residuals(fit2)    # Residuals
##D 
##D # Ordinal regression (proportional odds)
##D data$satisfaction <- ordered(sample(1:5, nrow(data), replace = TRUE))
##D fit3 <- gllamm(satisfaction ~ x + (1 | group),
##D                data = data,
##D                family = ordinal(link = "logit"))
##D 
##D # Ordinal with adjacent category logit
##D fit4 <- gllamm(satisfaction ~ x + (1 | group),
##D                data = data,
##D                family = ordinal(link = "acl"))
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gllamm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gllammr_validate")
### * gllammr_validate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gllammr_validate
### Title: Cross-package validation of GLLAMMR estimates
### Aliases: gllammr_validate

### ** Examples

## No test: 
if (requireNamespace("lme4", quietly = TRUE)) {
  gllammr_validate(cases = "gaussian_sleepstudy")
}
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gllammr_validate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("icc")
### * icc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: icc
### Title: Compute Intraclass Correlation Coefficients
### Aliases: icc

### ** Examples

## Not run: 
##D # All ICCs
##D icc(fit)
##D 
##D # Specific level
##D icc(fit, level = "class")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("icc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("irt")
### * irt

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: irt
### Title: IRT Family for Item Response Theory Models
### Aliases: irt

### ** Examples

## Not run: 
##D fit <- gllamm(response_matrix, family = irt("2PL"))
##D # Multi-level IRT: persons nested in classes
##D fit_ml <- gllamm(response_matrix, data = person_data,
##D                  family = irt("Rasch"), random = ~ (1 | class))
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("irt", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lca")
### * lca

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lca
### Title: Latent Class Family for Finite Mixture Models
### Aliases: lca

### ** Examples

## Not run: 
##D fit <- gllamm(indicator_matrix, family = lca(nclass = 3))
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lca", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("multinomial")
### * multinomial

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: multinomial
### Title: Multinomial Family for Unordered Categorical Outcomes
### Aliases: multinomial

### ** Examples

## Not run: 
##D fit <- gllamm(choice ~ x + (1 | region), data = d, family = multinomial())
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("multinomial", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ordinal")
### * ordinal

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ordinal
### Title: Ordinal Family for Proportional and Non-Proportional Odds Models
### Aliases: ordinal

### ** Examples

## Not run: 
##D # Proportional odds model (default)
##D family1 <- ordinal()
##D family2 <- ordinal(link = "logit")
##D 
##D # Adjacent category logit
##D family3 <- ordinal(link = "acl")
##D 
##D # Partial proportional odds
##D family4 <- ordinal(link = "ppo")
##D 
##D # Use with gllamm() - recommended interface
##D fit <- gllamm(rating ~ temp + (1 | judge),
##D               data = wine,
##D               family = ordinal(link = "logit"))
##D 
##D # Or use fit_ordinal() directly
##D fit2 <- fit_ordinal(rating ~ temp + (1 | judge),
##D                     data = wine,
##D                     link = "acl")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ordinal", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.gllamm_irt")
### * plot.gllamm_irt

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.gllamm_irt
### Title: Plot IRT Model Diagnostics
### Aliases: plot.gllamm_irt

### ** Examples

## Not run: 
##D # Fit 2PL model
##D fit <- fit_irt(responses, model = "2PL")
##D 
##D # Plot all diagnostics for items 1-3
##D plot(fit, which = 1:4, items = 1:3)
##D 
##D # Plot only ICCs
##D plot(fit, which = 1, items = 1:5)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.gllamm_irt", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.gllamm_lca")
### * plot.gllamm_lca

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.gllamm_lca
### Title: Plot Latent Class Analysis Results
### Aliases: plot.gllamm_lca

### ** Examples

## Not run: 
##D # Fit LCA model
##D fit <- fit_lca(indicators, nclass = 3)
##D 
##D # Plot all diagnostics
##D plot(fit, which = 1:3)
##D 
##D # Plot only class profiles
##D plot(fit, which = 1)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.gllamm_lca", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.gllamm_ordinal")
### * plot.gllamm_ordinal

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.gllamm_ordinal
### Title: Plot Ordinal Regression Model Diagnostics
### Aliases: plot.gllamm_ordinal

### ** Examples

## Not run: 
##D # Fit ordinal model
##D fit <- fit_ordinal(rating ~ temp + contact + (1 | judge),
##D                    data = wine, link = "logit")
##D 
##D # Plot all diagnostics for 'temp' covariate
##D plot(fit, which = 1:4, covariate = "temp")
##D 
##D # Plot only cumulative probabilities
##D plot(fit, which = 1, covariate = "contact")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.gllamm_ordinal", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_classification_uncertainty")
### * plot_classification_uncertainty

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_classification_uncertainty
### Title: Plot Individual Classification Uncertainty
### Aliases: plot_classification_uncertainty

### ** Examples

## Not run: 
##D fit <- fit_lca(data, nclass = 3)
##D plot_classification_uncertainty(fit, cases = 1:30)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_classification_uncertainty", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_item_covariates")
### * plot_item_covariates

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_item_covariates
### Title: Plot Item Covariate Effects
### Aliases: plot_item_covariates

### ** Examples

## Not run: 
##D fit <- fit_eirt(responses, item_data,
##D                 difficulty_formula = ~ word_freq)
##D 
##D plot_item_covariates(fit, covariate = "word_freq")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_item_covariates", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_ordinal_effects")
### * plot_ordinal_effects

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_ordinal_effects
### Title: Plot Ordinal Model Effects for Multiple Covariates
### Aliases: plot_ordinal_effects

### ** Examples

## Not run: 
##D fit <- fit_ordinal(rating ~ temp + contact + (1 | judge), data = wine)
##D plot_ordinal_effects(fit)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_ordinal_effects", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict.gllamm")
### * predict.gllamm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict.gllamm
### Title: Predict method for GLLAMM models
### Aliases: predict.gllamm

### ** Examples

## Not run: 
##D fit <- gllamm(y ~ x + (1 | group), data = mydata)
##D 
##D # Fitted values (default - conditional on random effects)
##D pred1 <- predict(fit)
##D 
##D # Population-level predictions (fixed effects only, u=0)
##D pred2 <- predict(fit, re.form = NA)
##D 
##D # Marginal predictions (population-averaged, integrating over u)
##D pred3 <- predict(fit, type = "marginal")
##D 
##D # Marginal predictions with standard errors
##D pred4 <- predict(fit, type = "marginal", se.fit = TRUE)
##D 
##D # Marginal predictions for new data
##D pred5 <- predict(fit, newdata = newdata, type = "marginal")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict.gllamm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict.gllamm_irt")
### * predict.gllamm_irt

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict.gllamm_irt
### Title: Predict method for IRT models
### Aliases: predict.gllamm_irt

### ** Examples

## Not run: 
##D # Fit 2PL model
##D responses <- matrix(rbinom(1000, 1, 0.6), 100, 10)
##D fit <- fit_irt(responses, model = "2PL")
##D 
##D # Person abilities
##D abilities <- predict(fit, type = "ability")
##D 
##D # Item response probabilities for each person
##D probs <- predict(fit, type = "probability")
##D 
##D # Marginal item response probabilities (population-level)
##D marg_probs <- predict(fit, type = "marginal")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict.gllamm_irt", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict.gllamm_ordinal")
### * predict.gllamm_ordinal

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict.gllamm_ordinal
### Title: Predict method for ordinal models
### Aliases: predict.gllamm_ordinal

### ** Examples

## Not run: 
##D # Fit ordinal model
##D fit <- gllamm(rating ~ temp + (1 | judge),
##D               data = wine,
##D               family = ordinal(link = "logit"))
##D 
##D # Predicted classes
##D pred_class <- predict(fit, type = "class")
##D 
##D # Conditional probabilities
##D pred_probs <- predict(fit, type = "probs")
##D 
##D # Marginal probabilities (population-averaged)
##D pred_marg <- predict(fit, type = "marginal")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict.gllamm_ordinal", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict_difficulty")
### * predict_difficulty

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict_difficulty
### Title: Predict Item Difficulties from Covariates
### Aliases: predict_difficulty

### ** Examples

## Not run: 
##D fit <- fit_eirt(responses, item_data,
##D                 difficulty_formula = ~ word_freq)
##D 
##D # Predicted difficulties for fitted data
##D pred_diff <- predict_difficulty(fit)
##D 
##D # Predicted difficulties for new items
##D new_items <- data.frame(word_freq = c(-1, 0, 1))
##D pred_new <- predict_difficulty(fit, newdata = new_items)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict_difficulty", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ranef.gllamm_irt_multilevel")
### * ranef.gllamm_irt_multilevel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ranef.gllamm_irt_multilevel
### Title: Extract Random Effects from Multi-Level IRT Models
### Aliases: ranef.gllamm_irt_multilevel

### ** Examples

## Not run: 
##D # Extract class effects
##D ranef(fit, level = "class")
##D 
##D # Extract all random effects
##D ranef(fit)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ranef.gllamm_irt_multilevel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulate.gllamm")
### * simulate.gllamm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulate.gllamm
### Title: Simulate from a GLLAMM model
### Aliases: simulate.gllamm

### ** Examples

## Not run: 
##D fit <- gllamm(y ~ x + (1 | group), data = mydata)
##D 
##D # Single simulation
##D sim1 <- simulate(fit)
##D 
##D # Multiple simulations
##D sim10 <- simulate(fit, nsim = 10)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulate.gllamm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("test_item_covariates")
### * test_item_covariates

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: test_item_covariates
### Title: Test Item Covariate Effects
### Aliases: test_item_covariates

### ** Examples

## Not run: 
##D item_data <- data.frame(
##D   word_freq = rnorm(20),
##D   length = rpois(20, 5)
##D )
##D 
##D # Test if word frequency matters
##D result <- test_item_covariates(
##D   responses,
##D   item_data,
##D   difficulty_formula = ~ word_freq + length,
##D   model = "Rasch"
##D )
##D 
##D print(result$comparison)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("test_item_covariates", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("test_proportional_odds")
### * test_proportional_odds

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: test_proportional_odds
### Title: Test Proportional Odds Assumption
### Aliases: test_proportional_odds

### ** Examples

## Not run: 
##D # Fit proportional odds model
##D fit_po <- fit_ordinal(rating ~ temp + (1 | judge),
##D                       data = wine, link = "logit")
##D 
##D # Test proportional odds assumption
##D po_test <- test_proportional_odds(fit_po)
##D print(po_test)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("test_proportional_odds", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
