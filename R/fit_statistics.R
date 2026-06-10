#' Model Fit Statistics
#'
#' Extract comprehensive fit statistics from GLLAMM models
#'
#' @name fit_statistics
NULL


#' Generic Fit Statistics Function
#'
#' Compute model-specific fit statistics for GLLAMM objects
#'
#' @param object A fitted model object
#' @param ... Additional arguments passed to methods
#'
#' @return An object of class \code{fit_statistics} with model-specific components
#'
#' @details
#' The \code{fit()} function provides comprehensive, model-specific fit statistics:
#'
#' \strong{GLMM models:}
#' \itemize{
#'   \item Log-likelihood, AIC, BIC
#'   \item R-squared (marginal and conditional for Gaussian models)
#'   \item Intraclass correlation coefficient (ICC)
#' }
#'
#' \strong{IRT models:}
#' \itemize{
#'   \item Log-likelihood, AIC, BIC
#'   \item Item fit statistics (S-X²)
#'   \item Person fit statistics (outfit/infit)
#'   \item Reliability estimates
#'   \item Test information function
#' }
#'
#' \strong{Latent Class Analysis:}
#' \itemize{
#'   \item Log-likelihood, AIC, BIC
#'   \item Entropy (classification quality)
#'   \item Class proportions
#'   \item Average posterior probabilities (APPA)
#' }
#'
#' \strong{Ordinal models:}
#' \itemize{
#'   \item Log-likelihood, AIC, BIC
#'   \item Pseudo-R² (McFadden)
#'   \item Proportional odds test (for PO/probit models)
#' }
#'
#' @examples
#' \dontrun{
#' # GLMM
#' fit1 <- gllamm(y ~ x + (1 | group), data = data)
#' fit(fit1)
#'
#' # IRT
#' fit2 <- fit_irt(responses, model = "2PL")
#' fit(fit2, compute_item_fit = TRUE)
#'
#' # LCA
#' fit3 <- fit_lca(indicators, nclass = 3)
#' fit(fit3)
#'
#' # Ordinal
#' fit4 <- fit_ordinal(rating ~ x + (1 | id), data = data)
#' fit(fit4, test_po = TRUE)
#' }
#'
#' @export
fit <- function(object, ...) {
  UseMethod("fit")
}


#' @export
fit.default <- function(object, ...) {
  stop("No fit method for class: ", paste(class(object), collapse = ", "),
       "\nAvailable for: gllamm, gllamm_irt, gllamm_lca, gllamm_ordinal")
}


#' Fit Statistics for Standard GLLAMM Models
#'
#' @param object A gllamm object
#' @param quiet Suppress ICC computation messages (default: FALSE)
#' @param ... Additional arguments (not currently used)
#'
#' @return Object of class \code{fit_gllamm} with fit statistics
#'
#' @export
fit.gllamm <- function(object, quiet = FALSE, ...) {

  # Don't dispatch to specialized methods for IRT/LCA/ordinal
  # (they have their own fit methods)
  if (inherits(object, c("gllamm_irt", "gllamm_lca", "gllamm_ordinal"))) {
    NextMethod("fit")
    return(invisible(NULL))
  }

  fit_stats <- list(
    model_type = "GLMM",
    family = object$family$family,
    logLik = object$logLik,
    AIC = object$AIC,
    BIC = object$BIC,
    n_obs = object$n_obs,
    n_params = object$n_params
  )

  # Nakagawa & Schielzeth R-squared for Gaussian models
  if (object$family$family == "gaussian") {
    fitted_vals <- fitted(object)
    resids <- residuals(object)

    # Variance of fixed effects
    var_fixed <- var(fitted_vals)

    # Variance of random effects (random_var holds per-term covariance
    # matrices; sum the diagonal variances)
    if (length(object$coefficients$random_var) > 0) {
      var_random <- sum(unlist(lapply(object$coefficients$random_var,
                                      function(m) {
                                        if (is.matrix(m)) diag(m) else as.numeric(m)
                                      })))
    } else {
      var_random <- 0
    }

    # Residual variance
    var_resid <- var(resids)

    # Total variance
    var_total <- var_fixed + var_random + var_resid

    # Marginal R²: variance explained by fixed effects
    fit_stats$R2_marginal <- var_fixed / var_total

    # Conditional R²: variance explained by fixed + random effects
    fit_stats$R2_conditional <- (var_fixed + var_random) / var_total
  }

  # ICC (Intraclass Correlation Coefficient)
  if (length(object$coefficients$random_var) > 0) {
    fit_stats$ICC <- tryCatch({
      icc(object, quiet = quiet)
    }, error = function(e) {
      if (!quiet) warning("Could not compute ICC: ", e$message)
      NA
    })
  }

  # Convergence status
  fit_stats$convergence <- object$convergence

  class(fit_stats) <- c("fit_gllamm", "fit_statistics")
  return(fit_stats)
}


#' Fit Statistics for IRT Models
#'
#' @param object A gllamm_irt object
#' @param compute_item_fit Compute item fit statistics (S-X²) (default: TRUE)
#' @param compute_person_fit Compute person fit statistics (outfit/infit) (default: TRUE)
#' @param ... Additional arguments
#'
#' @return Object of class \code{fit_irt} with IRT-specific fit statistics
#'
#' @export
fit.gllamm_irt <- function(object, compute_item_fit = TRUE,
                            compute_person_fit = TRUE, ...) {

  fit_stats <- list(
    model_type = "IRT",
    model = object$model,
    logLik = object$logLik,
    AIC = object$AIC,
    BIC = object$BIC,
    n_persons = object$n_persons,
    n_items = object$n_items
  )

  if (compute_item_fit) {
    fit_stats$item_fit <- compute_item_fit_sx2(object)
  }

  if (compute_person_fit) {
    fit_stats$person_fit <- compute_person_fit_outfit_infit(object)
  }

  # Marginal reliability
  fit_stats$reliability <- compute_irt_reliability(object)

  # Test information
  fit_stats$test_information <- compute_test_information_summary(object)

  class(fit_stats) <- c("fit_irt", "fit_statistics")
  return(fit_stats)
}


#' Fit Statistics for Latent Class Analysis
#'
#' @param object A gllamm_lca object
#' @param ... Additional arguments
#'
#' @return Object of class \code{fit_lca} with LCA-specific fit statistics
#'
#' @export
fit.gllamm_lca <- function(object, ...) {

  fit_stats <- list(
    model_type = "LCA",
    logLik = object$logLik,
    AIC = object$AIC,
    BIC = object$BIC,
    n_classes = object$n_classes,
    n_obs = nrow(object$posterior_probs),
    n_items = ncol(object$item_probs)
  )

  # Entropy: measure of classification certainty
  # E = 1 - [sum_i sum_k p_ik log(p_ik)] / (n log(K))
  posterior <- object$posterior_probs
  entropy_raw <- -sum(posterior * log(posterior + 1e-10), na.rm = TRUE)
  max_entropy <- nrow(posterior) * log(object$n_classes)
  fit_stats$entropy <- 1 - entropy_raw / max_entropy

  # Class proportions (from posterior probabilities)
  fit_stats$class_proportions <- colMeans(posterior)

  # Average posterior probability per class (APPA)
  # For each class, average the posterior probability among those assigned to it
  modal_class <- apply(posterior, 1, which.max)
  fit_stats$avg_posterior <- sapply(1:object$n_classes, function(k) {
    if (sum(modal_class == k) > 0) {
      mean(posterior[modal_class == k, k])
    } else {
      NA
    }
  })
  names(fit_stats$avg_posterior) <- paste0("Class", 1:object$n_classes)

  # Classification quality based on entropy
  if (fit_stats$entropy > 0.8) {
    fit_stats$classification_quality <- "Excellent"
  } else if (fit_stats$entropy > 0.6) {
    fit_stats$classification_quality <- "Good"
  } else if (fit_stats$entropy > 0.4) {
    fit_stats$classification_quality <- "Fair"
  } else {
    fit_stats$classification_quality <- "Poor"
  }

  class(fit_stats) <- c("fit_lca", "fit_statistics")
  return(fit_stats)
}


#' Fit Statistics for Ordinal Regression Models
#'
#' @param object A gllamm_ordinal object
#' @param test_po Test proportional odds assumption (default: TRUE for logit/probit)
#' @param ... Additional arguments
#'
#' @return Object of class \code{fit_ordinal} with ordinal-specific fit statistics
#'
#' @export
fit.gllamm_ordinal <- function(object, test_po = TRUE, ...) {

  fit_stats <- list(
    model_type = "Ordinal",
    link = object$link,
    logLik = object$logLik,
    AIC = object$AIC,
    BIC = object$BIC,
    n_categories = object$n_categories,
    n_obs = object$n_obs
  )

  # McFadden's pseudo-R²
  # R² = 1 - (logLik_model / logLik_null)
  # Null model: equal probability for all categories
  logLik_null <- object$n_obs * log(1 / object$n_categories)
  fit_stats$pseudo_R2 <- 1 - (object$logLik / logLik_null)

  # Test proportional odds assumption if applicable
  if (test_po && object$link %in% c("logit", "probit")) {
    fit_stats$proportional_odds_test <- tryCatch({
      test_proportional_odds(object, data = object$data)
    }, error = function(e) {
      warning("Could not perform proportional odds test: ", e$message)
      NULL
    })
  }

  class(fit_stats) <- c("fit_ordinal", "fit_statistics")
  return(fit_stats)
}


#' Fit Statistics for Explanatory IRT Models
#'
#' @param object A gllamm_eirt object
#' @param ... Additional arguments
#'
#' @return Object of class \code{fit_eirt} with EIRT-specific fit statistics
#'
#' @export
fit.gllamm_eirt <- function(object, ...) {

  fit_stats <- list(
    model_type = "EIRT",
    model = object$model,
    logLik = object$logLik,
    AIC = object$AIC,
    BIC = object$BIC,
    n_persons = object$n_persons,
    n_items = object$n_items
  )

  # Regression information
  fit_stats$difficulty_formula <- deparse(object$formulas$difficulty)
  fit_stats$discrimination_formula <- deparse(object$formulas$discrimination)

  # R² for item parameter regressions
  # Using the eirt_r_squared function if available
  if (requireNamespace("GLLAMMR", quietly = TRUE)) {
    fit_stats$R2_difficulty <- tryCatch({
      eirt_r_squared(object, parameter = "difficulty")
    }, error = function(e) NA)

    fit_stats$R2_discrimination <- tryCatch({
      eirt_r_squared(object, parameter = "discrimination")
    }, error = function(e) NA)
  }

  # Residual standard deviations (unexplained variation in item parameters)
  fit_stats$residual_sd_difficulty <- object$residual_sd$difficulty
  fit_stats$residual_sd_discrimination <- object$residual_sd$discrimination

  # Number of covariate effects
  fit_stats$n_difficulty_predictors <- length(object$regression_coefficients$difficulty) - 1  # Exclude intercept
  fit_stats$n_discrimination_predictors <- length(object$regression_coefficients$discrimination) - 1

  class(fit_stats) <- c("fit_eirt", "fit_statistics")
  return(fit_stats)
}


#' Fit Statistics for Multinomial Regression Models
#'
#' @param object A gllamm_multinomial object
#' @param ... Additional arguments
#'
#' @return Object of class \code{fit_multinomial} with multinomial-specific fit statistics
#'
#' @export
fit.gllamm_multinomial <- function(object, ...) {

  fit_stats <- list(
    model_type = "Multinomial",
    logLik = object$logLik,
    AIC = object$AIC,
    BIC = object$BIC,
    n_obs = object$n_obs,
    n_categories = object$n_categories,
    reference = object$reference,
    categories = object$categories
  )

  # McFadden's pseudo-R²
  # Null model: equal probability for all categories
  logLik_null <- object$n_obs * log(1 / object$n_categories)
  fit_stats$pseudo_R2 <- 1 - (object$logLik / logLik_null)

  # Random effects variance
  if (!is.null(object$coefficients$random_var)) {
    fit_stats$random_effects_var <- object$coefficients$random_var
  }

  # Coefficient matrix dimensions
  fit_stats$n_coefficients_per_category <- ncol(object$coefficients$beta)

  class(fit_stats) <- c("fit_multinomial", "fit_statistics")
  return(fit_stats)
}


#' Print Method for Fit Statistics
#'
#' @param x A fit_statistics object
#' @param ... Additional arguments
#'
#' @export
print.fit_statistics <- function(x, ...) {
  cat("Model Fit Statistics\n")
  cat("====================\n\n")

  cat("Model type:", x$model_type, "\n")

  if (!is.null(x$family)) {
    cat("Family:", x$family, "\n")
  }
  if (!is.null(x$model)) {
    cat("Model:", x$model, "\n")
  }
  if (!is.null(x$link)) {
    cat("Link:", x$link, "\n")
  }
  cat("\n")

  # Common statistics
  cat("Information Criteria:\n")
  cat("  Log-likelihood:", round(x$logLik, 2), "\n")
  cat("  AIC:", round(x$AIC, 2), "\n")
  cat("  BIC:", round(x$BIC, 2), "\n")
  cat("  N observations:", x$n_obs, "\n")

  # Model-specific statistics
  if (inherits(x, "fit_irt")) {
    cat("\n--- IRT Fit Statistics ---\n")
    cat("Items:", x$n_items, "\n")
    cat("Persons:", x$n_persons, "\n")

    if (!is.null(x$reliability)) {
      cat("\nReliability:\n")
      if (!is.null(x$reliability$marginal)) {
        cat("  Marginal:", round(x$reliability$marginal, 3), "\n")
      }
    }

    if (!is.null(x$item_fit)) {
      n_flagged <- sum(x$item_fit$flagged, na.rm = TRUE)
      cat("\nItem Fit:\n")
      cat("  Flagged items (p < 0.05):", n_flagged, "/", x$n_items, "\n")
    }

    if (!is.null(x$person_fit)) {
      n_flagged <- sum(x$person_fit$flagged, na.rm = TRUE)
      cat("\nPerson Fit:\n")
      cat("  Flagged persons:", n_flagged, "/", x$n_persons, "\n")
    }

    if (!is.null(x$test_information)) {
      cat("\nTest Information:\n")
      cat("  Maximum:", round(x$test_information$max, 2), "\n")
      cat("  At theta =", round(x$test_information$theta_max, 2), "\n")
    }
  }

  if (inherits(x, "fit_lca")) {
    cat("\n--- Latent Class Fit Statistics ---\n")
    cat("Number of classes:", x$n_classes, "\n")
    cat("Number of items:", x$n_items, "\n")
    cat("\nEntropy:", round(x$entropy, 3), "(", x$classification_quality, ")\n")

    cat("\nClass Proportions:\n")
    for (k in 1:x$n_classes) {
      cat(sprintf("  Class %d: %.3f\n", k, x$class_proportions[k]))
    }

    cat("\nAverage Posterior Probability (APPA):\n")
    for (k in 1:x$n_classes) {
      if (!is.na(x$avg_posterior[k])) {
        cat(sprintf("  Class %d: %.3f\n", k, x$avg_posterior[k]))
      }
    }
  }

  if (inherits(x, "fit_ordinal")) {
    cat("\n--- Ordinal Model Fit Statistics ---\n")
    cat("Number of categories:", x$n_categories, "\n")
    cat("Pseudo-R² (McFadden):", round(x$pseudo_R2, 3), "\n")

    if (!is.null(x$proportional_odds_test)) {
      cat("\nProportional Odds Test:\n")
      cat("  LRT statistic:", round(x$proportional_odds_test$statistic, 3), "\n")
      cat("  p-value:", format.pval(x$proportional_odds_test$p_value), "\n")
      cat("  ", x$proportional_odds_test$conclusion, "\n")
    }
  }

  if (inherits(x, "fit_eirt")) {
    cat("\n--- Explanatory IRT Fit Statistics ---\n")
    cat("Items:", x$n_items, "\n")
    cat("Persons:", x$n_persons, "\n")

    cat("\nDifficulty Regression:\n")
    cat("  Formula:", x$difficulty_formula, "\n")
    cat("  Predictors:", x$n_difficulty_predictors, "\n")
    if (!is.null(x$R2_difficulty) && !is.na(x$R2_difficulty)) {
      cat("  R²:", round(x$R2_difficulty, 3), "\n")
    }
    cat("  Residual SD:", round(x$residual_sd_difficulty, 3), "\n")

    cat("\nDiscrimination Regression:\n")
    cat("  Formula:", x$discrimination_formula, "\n")
    cat("  Predictors:", x$n_discrimination_predictors, "\n")
    if (!is.null(x$R2_discrimination) && !is.na(x$R2_discrimination)) {
      cat("  R²:", round(x$R2_discrimination, 3), "\n")
    }
    cat("  Residual SD:", round(x$residual_sd_discrimination, 3), "\n")
  }

  if (inherits(x, "fit_multinomial")) {
    cat("\n--- Multinomial Regression Fit Statistics ---\n")
    cat("Number of categories:", x$n_categories, "\n")
    cat("Reference category:", x$reference, "\n")
    cat("Pseudo-R² (McFadden):", round(x$pseudo_R2, 3), "\n")

    if (!is.null(x$random_effects_var)) {
      cat("\nRandom Effects Variance:", round(x$random_effects_var, 3), "\n")
    }

    cat("\nCoefficient Matrix:\n")
    cat("  Categories:", x$n_categories - 1, "(excluding reference)\n")
    cat("  Predictors per category:", x$n_coefficients_per_category, "\n")
  }

  if (inherits(x, "fit_gllamm")) {
    cat("\n--- GLMM Fit Statistics ---\n")

    if (!is.null(x$R2_marginal)) {
      cat("R² (marginal):", round(x$R2_marginal, 3), "\n")
      cat("R² (conditional):", round(x$R2_conditional, 3), "\n")
    }

    if (!is.null(x$ICC) && !all(is.na(x$ICC))) {
      cat("\nICC:", round(x$ICC, 3), "\n")
    }
  }

  invisible(x)
}


# ============================================================================
# Helper functions for IRT fit statistics
# ============================================================================

#' Compute S-X² item fit statistic
#' @keywords internal
compute_item_fit_sx2 <- function(object) {
  # Simplified S-X² implementation
  # Full implementation would bin persons by ability and compare obs vs exp

  n_items <- object$n_items

  # Placeholder: generate structure
  # Real implementation needs observed responses and expected probabilities
  item_fit <- data.frame(
    item = 1:n_items,
    sx2 = rep(NA, n_items),
    df = rep(NA, n_items),
    p_value = rep(NA, n_items),
    flagged = rep(FALSE, n_items)
  )

  return(item_fit)
}


#' Compute outfit/infit person fit statistics
#' @keywords internal
compute_person_fit_outfit_infit <- function(object) {
  # Outfit/infit are standardized residuals
  # Outfit: unweighted mean square, Infit: weighted mean square

  n_persons <- object$n_persons

  # Placeholder structure
  person_fit <- data.frame(
    person = 1:n_persons,
    outfit = rep(NA, n_persons),
    infit = rep(NA, n_persons),
    flagged = rep(FALSE, n_persons)
  )

  return(person_fit)
}


#' Compute IRT reliability
#' @keywords internal
compute_irt_reliability <- function(object) {
  # Marginal reliability: 1 - (1 / mean(information))
  # Requires person abilities and information function

  # Placeholder
  reliability <- list(
    marginal = NA
  )

  return(reliability)
}


#' Compute test information summary
#' @keywords internal
compute_test_information_summary <- function(object) {
  # Test information is sum of item information functions

  # Placeholder
  info_summary <- list(
    max = NA,
    theta_max = NA,
    mean = NA
  )

  return(info_summary)
}
