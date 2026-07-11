# Package index

## Model fitting

Top-level fitting entry points (the unified
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
interface plus the standalone eirt/DIF/NPML/LCA fitters) and the
`fit_*()` engines that
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
dispatches to internally but which can also be called directly.

- [`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
  : Fit Generalized Linear Latent and Mixed Models
- [`eirt()`](https://drjoshmcgrane.github.io/gllammr/reference/eirt.md)
  : Explanatory IRT Family
- [`dif_irt()`](https://drjoshmcgrane.github.io/gllammr/reference/dif_irt.md)
  : Confirmatory model-based DIF (IRT likelihood-ratio tests)
- [`npml()`](https://drjoshmcgrane.github.io/gllammr/reference/npml.md)
  : Nonparametric maximum likelihood integration
- [`fit()`](https://drjoshmcgrane.github.io/gllammr/reference/fit.md) :
  Generic Fit Statistics Function
- [`fit_binomial()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_binomial.md)
  : Fit Binomial Regression Models with Random Effects
- [`fit_cdm()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_cdm.md)
  : Fit Cognitive Diagnosis Models
- [`fit_eirt()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_eirt.md)
  : Fit Explanatory Item Response Theory Models
- [`fit_irt()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_irt.md)
  : Fit Item Response Theory Models
- [`fit_lca()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_lca.md)
  : Fit Latent Class Analysis Models
- [`fit_mixed()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_mixed.md)
  : Fit Joint Models for Mixed Response Types
- [`fit_multinomial()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_multinomial.md)
  : Fit Multinomial Regression Models with Random Effects
- [`fit_npml()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_npml.md)
  : Fit Two-Level GLMMs by Nonparametric Maximum Likelihood (NPML)
- [`fit_ordinal()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_ordinal.md)
  : Fit Ordinal Regression Models with Random Effects
- [`fit_rank()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_rank.md)
  : Fit Rank-Ordered Logit Models with Taste Heterogeneity
- [`fit_sem()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_sem.md)
  : Fit Structural Equation Models with Latent Variables
- [`fit_survival()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_survival.md)
  : Fit Parametric Survival Models with Random Effects (Frailty)

## Family and spec constructors

Family objects passed as `gllamm(family = ...)`, plus the adaptive
quadrature integration specification used as
`gllamm(integration = ...)`.

- [`ordinal()`](https://drjoshmcgrane.github.io/gllammr/reference/ordinal.md)
  : Ordinal Family for Proportional and Non-Proportional Odds Models
- [`binomial()`](https://drjoshmcgrane.github.io/gllammr/reference/binomial.md)
  : Binomial Family for Binary and Binomial Outcomes
- [`multinomial()`](https://drjoshmcgrane.github.io/gllammr/reference/multinomial.md)
  : Multinomial Family for Unordered Categorical Outcomes
- [`irt()`](https://drjoshmcgrane.github.io/gllammr/reference/irt.md) :
  IRT Family for Item Response Theory Models
- [`lca()`](https://drjoshmcgrane.github.io/gllammr/reference/lca.md) :
  Latent Class Family for Finite Mixture Models
- [`cdm()`](https://drjoshmcgrane.github.io/gllammr/reference/cdm.md) :
  Cognitive Diagnosis Family for Q-Matrix Models
- [`sem()`](https://drjoshmcgrane.github.io/gllammr/reference/sem.md) :
  SEM Family for Structural Equation Models
- [`mixed_response()`](https://drjoshmcgrane.github.io/gllammr/reference/mixed_response.md)
  : Mixed-Response Family for Joint Outcome Models
- [`ranking()`](https://drjoshmcgrane.github.io/gllammr/reference/ranking.md)
  : Rank-Ordered Logit Family
- [`survival_family()`](https://drjoshmcgrane.github.io/gllammr/reference/survival_family.md)
  : Parametric Frailty Survival Family
- [`aghq()`](https://drjoshmcgrane.github.io/gllammr/reference/aghq.md)
  : Adaptive quadrature integration specification

## Methods and diagnostics

Extractors, plots, model comparison, and diagnostic/testing helpers for
fitted `gllammr` models.

- [`VarCorr()`](https://drjoshmcgrane.github.io/gllammr/reference/VarCorr.md)
  : Extract Variance Components from Multi-Level Models
- [`abilities()`](https://drjoshmcgrane.github.io/gllammr/reference/abilities.md)
  : Extract Person Abilities from IRT Models
- [`fixef()`](https://drjoshmcgrane.github.io/gllammr/reference/fixef.md)
  : Generic fixef
- [`ranef()`](https://drjoshmcgrane.github.io/gllammr/reference/ranef.md)
  : Generic ranef
- [`icc()`](https://drjoshmcgrane.github.io/gllammr/reference/icc.md) :
  Compute Intraclass Correlation Coefficients
- [`gof.gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gof.gllamm.md)
  : Goodness-of-fit tests for GLLAMM
- [`predict_difficulty()`](https://drjoshmcgrane.github.io/gllammr/reference/predict_difficulty.md)
  : Predict Item Difficulties from Covariates
- [`find_outliers()`](https://drjoshmcgrane.github.io/gllammr/reference/find_outliers.md)
  : Outlier detection for GLLAMM
- [`plot_classification_uncertainty()`](https://drjoshmcgrane.github.io/gllammr/reference/plot_classification_uncertainty.md)
  : Plot Individual Classification Uncertainty
- [`plot_item_covariates()`](https://drjoshmcgrane.github.io/gllammr/reference/plot_item_covariates.md)
  : Plot Item Covariate Effects
- [`plot_ordinal_effects()`](https://drjoshmcgrane.github.io/gllammr/reference/plot_ordinal_effects.md)
  : Plot Ordinal Model Effects for Multiple Covariates
- [`dif_plot()`](https://drjoshmcgrane.github.io/gllammr/reference/dif_plot.md)
  : Plot item response curves by DIF group
- [`dif_test()`](https://drjoshmcgrane.github.io/gllammr/reference/dif_test.md)
  : Test for Differential Item Functioning (DIF)
- [`dif_test_with_data()`](https://drjoshmcgrane.github.io/gllammr/reference/dif_test_with_data.md)
  : Test for DIF with explicit response data (deprecated)
- [`compare_eirt()`](https://drjoshmcgrane.github.io/gllammr/reference/compare_eirt.md)
  : Compare Explanatory IRT Models
- [`compare_models()`](https://drjoshmcgrane.github.io/gllammr/reference/compare_models.md)
  : Compare fitted gllammr models
- [`latent_structure_comparison()`](https://drjoshmcgrane.github.io/gllammr/reference/latent_structure_comparison.md)
  : Latent structure comparison: categorization, ordering or
  quantification
- [`eirt_r_squared()`](https://drjoshmcgrane.github.io/gllammr/reference/eirt_r_squared.md)
  : Compute R-squared for Item Parameter Regression
- [`test_item_covariates()`](https://drjoshmcgrane.github.io/gllammr/reference/test_item_covariates.md)
  : Test Item Covariate Effects
- [`test_proportional_odds()`](https://drjoshmcgrane.github.io/gllammr/reference/test_proportional_odds.md)
  : Test Proportional Odds Assumption

## Validation and utilities

Cross-package validation harness.

- [`gllammr_validate()`](https://drjoshmcgrane.github.io/gllammr/reference/gllammr_validate.md)
  : Cross-package validation of gllammr estimates

## Internal

S3 methods (print/plot/predict/simulate/coef/VarCorr/… implementations
for each model class) and unexported helper functions. Documented for
maintainer reference but not part of the public modeling API described
in the sections above.

- [`.dif_item_fullfit()`](https://drjoshmcgrane.github.io/gllammr/reference/dot-dif_item_fullfit.md)
  : Full-model (M2) fit for one item, kept for plotting
- [`.dif_item_tests()`](https://drjoshmcgrane.github.io/gllammr/reference/dot-dif_item_tests.md)
  : Nested-model LR tests for one item
- [`.gllamm_re_parts()`](https://drjoshmcgrane.github.io/gllammr/reference/dot-gllamm_re_parts.md)
  : Per-term random-effects pieces of a fitted GLMM
- [`.isotonic_poset()`](https://drjoshmcgrane.github.io/gllammr/reference/dot-isotonic_poset.md)
  : Weighted isotonic regression over a partial order
- [`.mvn_em_saturated()`](https://drjoshmcgrane.github.io/gllammr/reference/dot-mvn_em_saturated.md)
  : EM for the saturated multivariate-normal model under missingness
- [`.ordinal_category_probs()`](https://drjoshmcgrane.github.io/gllammr/reference/dot-ordinal_category_probs.md)
  : Category probabilities for every ordinal link
- [`.ordinal_re_parts()`](https://drjoshmcgrane.github.io/gllammr/reference/dot-ordinal_re_parts.md)
  : Per-term random-effects pieces of a fitted ordinal model
- [`.pava_weighted()`](https://drjoshmcgrane.github.io/gllammr/reference/dot-pava_weighted.md)
  : Weighted isotonic regression (pool-adjacent-violators)
- [`.rasch_mstep()`](https://drjoshmcgrane.github.io/gllammr/reference/dot-rasch_mstep.md)
  : Weighted logistic M-step for the located latent class (Rasch) model
- [`.rmsea_ci()`](https://drjoshmcgrane.github.io/gllammr/reference/dot-rmsea_ci.md)
  : RMSEA 90 percent confidence interval (noncentral chi-square
  inversion)
- [`.topological_order()`](https://drjoshmcgrane.github.io/gllammr/reference/dot-topological_order.md)
  : Topological order of classes under a partial order (Kahn's
  algorithm)
- [`.val_row()`](https://drjoshmcgrane.github.io/gllammr/reference/dot-val_row.md)
  : Build one validation result row
- [`abilities(`*`<gllamm_irt_multilevel>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/abilities.gllamm_irt_multilevel.md)
  : Total abilities for multi-level IRT fits
- [`align_weights()`](https://drjoshmcgrane.github.io/gllammr/reference/align_weights.md)
  : Align observation weights with listwise-deleted model data
- [`coef(`*`<gllamm_eirt>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/coef.gllamm_eirt.md)
  : Extract Item Parameters from EIRT Model
- [`coef(`*`<gllamm_irt>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/coef.gllamm_irt.md)
  : Extract Coefficients from IRT Models
- [`compute_irt_reliability()`](https://drjoshmcgrane.github.io/gllammr/reference/compute_irt_reliability.md)
  : Compute IRT reliability
- [`compute_item_fit_sx2()`](https://drjoshmcgrane.github.io/gllammr/reference/compute_item_fit_sx2.md)
  : Compute S-X^2 item fit statistic
- [`compute_multinomial_probs()`](https://drjoshmcgrane.github.io/gllammr/reference/compute_multinomial_probs.md)
  : Compute multinomial probabilities
- [`compute_person_fit_outfit_infit()`](https://drjoshmcgrane.github.io/gllammr/reference/compute_person_fit_outfit_infit.md)
  : Compute outfit/infit person fit statistics
- [`compute_test_information_summary()`](https://drjoshmcgrane.github.io/gllammr/reference/compute_test_information_summary.md)
  : Compute test information summary
- [`construct_Z_matrix()`](https://drjoshmcgrane.github.io/gllammr/reference/construct_Z_matrix.md)
  : Construct random effects design matrix for newdata
- [`cooks.distance(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/cooks.distance.gllamm.md)
  : Group-level Cook's distance for GLLAMM models
- [`create_grouping_matrix()`](https://drjoshmcgrane.github.io/gllammr/reference/create_grouping_matrix.md)
  : Create grouping factor matrix
- [`diagnostics`](https://drjoshmcgrane.github.io/gllammr/reference/diagnostics.md)
  : Diagnostic Methods for GLLAMM Models
- [`drop_intercept_column()`](https://drjoshmcgrane.github.io/gllammr/reference/drop_intercept_column.md)
  : Drop the intercept column from a fixed-effects design matrix
- [`eirt_item_thresholds()`](https://drjoshmcgrane.github.io/gllammr/reference/eirt_item_thresholds.md)
  : Reconstruct per-item threshold parameters from a fitted polytomous
  EIRT model
- [`eirt_poly_model_name()`](https://drjoshmcgrane.github.io/gllammr/reference/eirt_poly_model_name.md)
  : Map an EIRT poly_model_type code to the shared probability helper's
  model name
- [`expand_nested_random_term()`](https://drjoshmcgrane.github.io/gllammr/reference/expand_nested_random_term.md)
  : Expand a nested random-effects term into one term per level
- [`expand_nested_terms()`](https://drjoshmcgrane.github.io/gllammr/reference/expand_nested_terms.md)
  : Expand nested terms
- [`extract_random_terms()`](https://drjoshmcgrane.github.io/gllammr/reference/extract_random_terms.md)
  : Extract random effects terms from formula
- [`extract_random_vcov()`](https://drjoshmcgrane.github.io/gllammr/reference/extract_random_vcov.md)
  : Extract random effects variance-covariance matrix from fitted model
- [`families`](https://drjoshmcgrane.github.io/gllammr/reference/families.md)
  : Family Objects for GLLAMM Models
- [`fit_cdm_em()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_cdm_em.md)
  : EM estimation for cognitive diagnosis models
- [`fit_irt_dichotomous()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_irt_dichotomous.md)
  : Internal function for dichotomous IRT models
- [`fit_irt_em()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_irt_em.md)
  : MML-EM estimation for IRT models (Bock-Aitkin)
- [`fit_irt_polytomous()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_irt_polytomous.md)
  : Internal function for polytomous IRT models
- [`fit_lca_em()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_lca_em.md)
  : EM estimation for latent class models
- [`fit_multinomial_multi()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_multinomial_multi.md)
  : Multinomial model with multiple random-effects terms
- [`fit_ordinal_multi()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_ordinal_multi.md)
  : Ordinal model with multiple random-effects terms
- [`fit_sem_ml()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_sem_ml.md)
  : Covariance-based ML estimation for SEM (Wishart / FIML likelihood)
- [`fit_statistics`](https://drjoshmcgrane.github.io/gllammr/reference/fit_statistics.md)
  : Model Fit Statistics
- [`fit_tmb_gllamm_aghq()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_tmb_gllamm_aghq.md)
  : Fit a two-level GLMM by adaptive Gauss-Hermite quadrature
- [`fit_tmb_gllamm_multi()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_tmb_gllamm_multi.md)
  : GLMM engine for multiple random-effects terms (crossed/nested)
- [`fit_tmb_gllamm_v2()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_tmb_gllamm_v2.md)
  : Enhanced interface to TMB for gllammr models
- [`fit_tmb_objective_only()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_tmb_objective_only.md)
  : Build a TMB objective (no optimization) for one cluster
- [`fit(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/fit.gllamm.md)
  : Fit Statistics for Standard GLLAMM Models
- [`fit(`*`<gllamm_eirt>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/fit.gllamm_eirt.md)
  : Fit Statistics for Explanatory IRT Models
- [`fit(`*`<gllamm_irt>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/fit.gllamm_irt.md)
  : Fit Statistics for IRT Models
- [`fit(`*`<gllamm_lca>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/fit.gllamm_lca.md)
  : Fit Statistics for Latent Class Analysis
- [`fit(`*`<gllamm_multinomial>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/fit.gllamm_multinomial.md)
  : Fit Statistics for Multinomial Regression Models
- [`fit(`*`<gllamm_ordinal>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/fit.gllamm_ordinal.md)
  : Fit Statistics for Ordinal Regression Models
- [`fixef(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/fixef.gllamm.md)
  : Extract fixed effects
- [`gauss_hermite()`](https://drjoshmcgrane.github.io/gllammr/reference/gauss_hermite.md)
  : Gauss-Hermite nodes and weights (Golub-Welsch)
- [`get_inverse_link()`](https://drjoshmcgrane.github.io/gllammr/reference/get_inverse_link.md)
  : Get inverse link function for a family
- [`print(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm-class.md)
  [`summary(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm-class.md)
  [`coef(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm-class.md)
  [`vcov(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm-class.md)
  [`logLik(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm-class.md)
  [`fitted(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm-class.md)
  [`residuals(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm-class.md)
  : GLLAMM model object class
- [`gllammr`](https://drjoshmcgrane.github.io/gllammr/reference/gllammr-package.md)
  [`gllammr-package`](https://drjoshmcgrane.github.io/gllammr/reference/gllammr-package.md)
  : gllammr: Generalized Linear Latent and Mixed Models
- [`icc(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/icc.gllamm.md)
  : Variance decomposition (ICC)
- [`influence(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/influence.gllamm.md)
  : Influence diagnostics for GLLAMM
- [`irt_category_probs()`](https://drjoshmcgrane.github.io/gllammr/reference/irt_category_probs.md)
  : Category response probabilities for polytomous IRT models
- [`make_model_matrices()`](https://drjoshmcgrane.github.io/gllammr/reference/make_model_matrices.md)
  : Extract model matrices
- [`marginal_utils`](https://drjoshmcgrane.github.io/gllammr/reference/marginal_utils.md)
  : Utility Functions for Marginal Predictions
- [`mc_integrate_fixed_samples()`](https://drjoshmcgrane.github.io/gllammr/reference/mc_integrate_fixed_samples.md)
  : Monte Carlo integration for marginal predictions
- [`mc_integrate_marginal()`](https://drjoshmcgrane.github.io/gllammr/reference/mc_integrate_marginal.md)
  : Monte Carlo integration for marginal predictions (one sample at a
  time)
- [`parse_formula()`](https://drjoshmcgrane.github.io/gllammr/reference/parse_formula.md)
  : Parse GLLAMM formula
- [`parse_level_weights()`](https://drjoshmcgrane.github.io/gllammr/reference/parse_level_weights.md)
  : Parse observation-level and group-level survey weights
- [`parse_random_formula()`](https://drjoshmcgrane.github.io/gllammr/reference/parse_random_formula.md)
  : Parse Random Effects Formula
- [`parse_random_term()`](https://drjoshmcgrane.github.io/gllammr/reference/parse_random_term.md)
  : Parse a single random effects term
- [`parse_single_random_term()`](https://drjoshmcgrane.github.io/gllammr/reference/parse_single_random_term.md)
  : Parse a single random effects term
- [`plot_ability_distribution_irt()`](https://drjoshmcgrane.github.io/gllammr/reference/plot_ability_distribution_irt.md)
  : Plot Person Ability Distribution
- [`plot_category_probs_ordinal()`](https://drjoshmcgrane.github.io/gllammr/reference/plot_category_probs_ordinal.md)
  : Plot Category Probabilities
- [`plot_class_profiles_lca()`](https://drjoshmcgrane.github.io/gllammr/reference/plot_class_profiles_lca.md)
  : Plot Class Profiles
- [`plot_classification_lca()`](https://drjoshmcgrane.github.io/gllammr/reference/plot_classification_lca.md)
  : Plot Classification Summary
- [`plot_covariate_effects_ordinal()`](https://drjoshmcgrane.github.io/gllammr/reference/plot_covariate_effects_ordinal.md)
  : Plot Covariate Effects
- [`plot_cumulative_probs_ordinal()`](https://drjoshmcgrane.github.io/gllammr/reference/plot_cumulative_probs_ordinal.md)
  : Plot Cumulative Probabilities
- [`plot_icc_irt()`](https://drjoshmcgrane.github.io/gllammr/reference/plot_icc_irt.md)
  : Plot Item Characteristic Curves
- [`plot_iif_irt()`](https://drjoshmcgrane.github.io/gllammr/reference/plot_iif_irt.md)
  : Plot Item Information Functions
- [`plot_irt`](https://drjoshmcgrane.github.io/gllammr/reference/plot_irt.md)
  : Plotting Functions for IRT Models
- [`plot_item_probabilities_lca()`](https://drjoshmcgrane.github.io/gllammr/reference/plot_item_probabilities_lca.md)
  : Plot Item Probability Heatmap
- [`plot_lca`](https://drjoshmcgrane.github.io/gllammr/reference/plot_lca.md)
  : Plotting Functions for Latent Class Analysis
- [`plot_ordinal`](https://drjoshmcgrane.github.io/gllammr/reference/plot_ordinal.md)
  : Plotting Functions for Ordinal Regression Models
- [`plot_thresholds_ordinal()`](https://drjoshmcgrane.github.io/gllammr/reference/plot_thresholds_ordinal.md)
  : Plot Threshold Parameters
- [`plot_tif_irt()`](https://drjoshmcgrane.github.io/gllammr/reference/plot_tif_irt.md)
  : Plot Test Information Function
- [`plot(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/plot.gllamm.md)
  : Plot diagnostics for GLLAMM models
- [`plot(`*`<gllamm_irt>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/plot.gllamm_irt.md)
  : Plot IRT Model Diagnostics
- [`plot(`*`<gllamm_lca>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/plot.gllamm_lca.md)
  : Plot Latent Class Analysis Results
- [`plot(`*`<gllamm_ordinal>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/plot.gllamm_ordinal.md)
  : Plot Ordinal Regression Model Diagnostics
- [`predict_marginal_gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/predict_marginal_gllamm.md)
  : Internal function for marginal predictions
- [`predict_marginal_irt()`](https://drjoshmcgrane.github.io/gllammr/reference/predict_marginal_irt.md)
  : Internal function for marginal IRT predictions
- [`predict_marginal_ordinal()`](https://drjoshmcgrane.github.io/gllammr/reference/predict_marginal_ordinal.md)
  : Internal function for marginal ordinal predictions
- [`predict(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/predict.gllamm.md)
  : Predict method for GLLAMM models
- [`predict(`*`<gllamm_eirt>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/predict.gllamm_eirt.md)
  : Predict method for EIRT models
- [`predict(`*`<gllamm_irt>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/predict.gllamm_irt.md)
  : Predict method for IRT models
- [`predict(`*`<gllamm_irt_poly>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/predict.gllamm_irt_poly.md)
  : Predicted category probabilities for polytomous IRT fits
- [`predict(`*`<gllamm_mixed>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/predict.gllamm_mixed.md)
  : Predict from a fitted mixed-response model
- [`predict(`*`<gllamm_multinomial>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/predict_multinomial.md)
  : Predict method for multinomial models
- [`predict(`*`<gllamm_npml>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/predict.gllamm_npml.md)
  : Predict from a fitted NPML model
- [`predict(`*`<gllamm_ordinal>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/predict.gllamm_ordinal.md)
  : Predict method for ordinal models
- [`predict(`*`<gllamm_survival>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/predict_survival.md)
  : Predict method for survival models
- [`print(`*`<binomial_family>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/print.binomial_family.md)
  : Print method for binomial family
- [`print(`*`<eirt_comparison>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/print.eirt_comparison.md)
  : Print EIRT comparison
- [`print(`*`<eirt_test>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/print.eirt_test.md)
  : Print EIRT test results
- [`print(`*`<fit_statistics>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/print.fit_statistics.md)
  : Print Method for Fit Statistics
- [`print(`*`<gllamm_eirt>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/print.gllamm_eirt.md)
  : Print EIRT model results
- [`print(`*`<ordinal_family>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/print.ordinal_family.md)
  : Print method for ordinal family
- [`print(`*`<po_test>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/print.po_test.md)
  : Print method for proportional odds test
- [`ranef(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/ranef.gllamm.md)
  : Extract random effects
- [`ranef(`*`<gllamm_irt_multilevel>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/ranef.gllamm_irt_multilevel.md)
  : Extract Random Effects from Multi-Level IRT Models
- [`rmvnorm_chol()`](https://drjoshmcgrane.github.io/gllammr/reference/rmvnorm_chol.md)
  : Draw samples from a zero-mean multivariate normal
- [`sandwich_vcov_gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/sandwich_vcov_gllamm.md)
  : Cluster-robust (sandwich) covariance for GLLAMM fixed effects
- [`simulate(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/simulate.gllamm.md)
  : Simulate from a GLLAMM model
- [`simulate(`*`<gllamm_cdm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/simulate.gllamm_cdm.md)
  : Simulate response matrices from a fitted cognitive diagnosis model
- [`simulate(`*`<gllamm_eirt>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/simulate.gllamm_eirt.md)
  : Simulate response matrices from a fitted explanatory IRT model
- [`simulate(`*`<gllamm_irt>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/simulate.gllamm_irt.md)
  : Simulate response matrices from a fitted IRT model
- [`simulate(`*`<gllamm_lca>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/simulate.gllamm_lca.md)
  : Simulate response matrices from a fitted latent class model
- [`simulate(`*`<gllamm_mixed>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/simulate.gllamm_mixed.md)
  : Simulate from a fitted mixed-response model
- [`simulate(`*`<gllamm_multinomial>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/simulate.gllamm_multinomial.md)
  : Simulate from a fitted multinomial model
- [`simulate(`*`<gllamm_npml>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/simulate.gllamm_npml.md)
  : Simulate from a fitted NPML model
- [`simulate(`*`<gllamm_ordinal>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/simulate.gllamm_ordinal.md)
  : Simulate from a fitted ordinal model
- [`simulate(`*`<gllamm_sem>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/simulate.gllamm_sem.md)
  : Simulate indicator data from a fitted SEM
- [`simulate(`*`<gllamm_survival>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/simulate.gllamm_survival.md)
  : Simulate from a fitted parametric frailty survival model
- [`summary(`*`<gllamm_eirt>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/summary.gllamm_eirt.md)
  : Summary of EIRT model
- [`validate_formula()`](https://drjoshmcgrane.github.io/gllammr/reference/validate_formula.md)
  : Validate formula
- [`validate_poly_responses()`](https://drjoshmcgrane.github.io/gllammr/reference/validate_poly_responses.md)
  : Validate and auto-recode polytomous response matrices (shared by the
  Laplace and EM estimation paths)
- [`VarCorr(`*`<gllamm>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/VarCorr.gllamm.md)
  : Extract variance components
- [`VarCorr(`*`<gllamm_irt_multilevel>`*`)`](https://drjoshmcgrane.github.io/gllammr/reference/VarCorr.gllamm_irt_multilevel.md)
  : Extract Variance Components from Multi-Level IRT Models
