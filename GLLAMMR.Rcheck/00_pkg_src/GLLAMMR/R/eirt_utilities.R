#' Compare Explanatory IRT Models
#'
#' Compare two or more EIRT models using likelihood ratio tests and information criteria
#'
#' @param ... Two or more gllamm_eirt objects to compare
#' @param test Type of test: "LRT" for likelihood ratio test, "none" for just IC comparison
#'
#' @return A data frame with model comparison statistics
#'
#' @details
#' This function compares EIRT models to test the importance of item covariates.
#' Models are compared using:
#' \itemize{
#'   \item Log-likelihood
#'   \item AIC and BIC
#'   \item Likelihood ratio test (if nested models)
#' }
#'
#' For nested models (e.g., with and without a predictor), the LRT tests
#' whether adding the predictor significantly improves fit.
#'
#' @examples
#' \dontrun{
#' # Fit model with predictor
#' fit1 <- fit_eirt(responses, item_data,
#'                  difficulty_formula = ~ word_freq,
#'                  model = "Rasch")
#'
#' # Fit model without predictor
#' fit0 <- fit_eirt(responses, item_data,
#'                  difficulty_formula = ~ 1,
#'                  model = "Rasch")
#'
#' # Compare models
#' compare_eirt(fit0, fit1)
#' }
#'
#' @export
compare_eirt <- function(..., test = c("LRT", "none")) {
  test <- match.arg(test)

  models <- list(...)

  # Validate inputs
  if (length(models) < 2) {
    stop("Need at least 2 models to compare")
  }

  for (i in seq_along(models)) {
    if (!inherits(models[[i]], "gllamm_eirt")) {
      stop("All models must be of class 'gllamm_eirt'")
    }
  }

  # Extract model information
  n_models <- length(models)
  model_names <- as.character(substitute(list(...))[-1])

  comparison <- data.frame(
    Model = model_names,
    npar = sapply(models, function(m) length(m$tmb_obj$par)),
    logLik = sapply(models, function(m) m$logLik),
    AIC = sapply(models, function(m) m$AIC),
    BIC = sapply(models, function(m) m$BIC),
    stringsAsFactors = FALSE
  )

  # Compute differences from first model
  comparison$delta_AIC <- comparison$AIC - comparison$AIC[1]
  comparison$delta_BIC <- comparison$BIC - comparison$BIC[1]

  # Add LRT if requested
  if (test == "LRT" && n_models == 2) {
    # Check if models are nested
    npar1 <- comparison$npar[1]
    npar2 <- comparison$npar[2]

    if (npar1 != npar2) {
      # Assume nested
      if (npar1 < npar2) {
        smaller <- 1
        larger <- 2
      } else {
        smaller <- 2
        larger <- 1
      }

      lrt_stat <- 2 * (comparison$logLik[larger] - comparison$logLik[smaller])
      df_diff <- abs(npar2 - npar1)
      p_value <- pchisq(lrt_stat, df_diff, lower.tail = FALSE)

      comparison$LRT_stat <- c(NA, lrt_stat)[c(smaller, larger)]
      comparison$LRT_df <- c(NA, df_diff)[c(smaller, larger)]
      comparison$LRT_p <- c(NA, p_value)[c(smaller, larger)]
    }
  }

  class(comparison) <- c("eirt_comparison", "data.frame")
  return(comparison)
}


#' Print EIRT comparison
#' @keywords internal
#' @export
print.eirt_comparison <- function(x, ...) {
  cat("EIRT Model Comparison\n")
  cat("=====================\n\n")

  # Print model info
  print(as.data.frame(x), row.names = FALSE)

  cat("\n")

  # Interpret results
  best_aic <- which.min(x$AIC)
  best_bic <- which.min(x$BIC)

  cat("Best model by AIC:", x$Model[best_aic], "\n")
  cat("Best model by BIC:", x$Model[best_bic], "\n")

  if ("LRT_p" %in% names(x) && !all(is.na(x$LRT_p))) {
    cat("\nLikelihood Ratio Test:\n")
    p_val <- x$LRT_p[!is.na(x$LRT_p)]
    if (p_val < 0.001) {
      cat("  p < 0.001 (highly significant improvement)\n")
    } else if (p_val < 0.05) {
      cat("  p =", round(p_val, 4), "(significant improvement)\n")
    } else {
      cat("  p =", round(p_val, 4), "(no significant improvement)\n")
    }
  }

  invisible(x)
}


#' Test Item Covariate Effects
#'
#' Test whether adding item covariates significantly improves model fit
#'
#' @param response_matrix Matrix of item responses
#' @param item_data Data frame of item-level covariates
#' @param difficulty_formula Formula for difficulty regression
#' @param discrimination_formula Formula for discrimination regression
#' @param model IRT model type
#' @param ... Additional arguments passed to fit_eirt
#'
#' @return A list with:
#'   \item{full_model}{EIRT model with covariates}
#'   \item{null_model}{EIRT model without covariates (intercept only)}
#'   \item{comparison}{Model comparison results}
#'
#' @examples
#' \dontrun{
#' item_data <- data.frame(
#'   word_freq = rnorm(20),
#'   length = rpois(20, 5)
#' )
#'
#' # Test if word frequency matters
#' result <- test_item_covariates(
#'   responses,
#'   item_data,
#'   difficulty_formula = ~ word_freq + length,
#'   model = "Rasch"
#' )
#'
#' print(result$comparison)
#' }
#'
#' @export
test_item_covariates <- function(response_matrix,
                                 item_data,
                                 difficulty_formula = ~ 1,
                                 discrimination_formula = ~ 1,
                                 model = c("Rasch", "2PL", "GRM"),
                                 ...) {

  model <- match.arg(model)

  cat("Fitting null model (intercept only)...\n")
  null_model <- fit_eirt(
    response_matrix = response_matrix,
    item_data = item_data,
    difficulty_formula = ~ 1,
    discrimination_formula = ~ 1,
    model = model,
    ...
  )

  cat("Fitting full model with covariates...\n")
  full_model <- fit_eirt(
    response_matrix = response_matrix,
    item_data = item_data,
    difficulty_formula = difficulty_formula,
    discrimination_formula = discrimination_formula,
    model = model,
    ...
  )

  cat("\nComparing models...\n")
  comparison <- compare_eirt(null_model, full_model, test = "LRT")

  result <- list(
    null_model = null_model,
    full_model = full_model,
    comparison = comparison
  )

  class(result) <- "eirt_test"
  return(result)
}


#' Print EIRT test results
#' @keywords internal
#' @export
print.eirt_test <- function(x, ...) {
  cat("Item Covariate Effects Test\n")
  cat("============================\n\n")

  cat("Null model (intercept only):\n")
  cat("  logLik:", round(x$null_model$logLik, 2), "\n")
  cat("  AIC:", round(x$null_model$AIC, 2), "\n\n")

  cat("Full model:\n")
  cat("  Difficulty formula:", deparse(x$full_model$formulas$difficulty), "\n")
  cat("  Discrimination formula:", deparse(x$full_model$formulas$discrimination), "\n")
  cat("  logLik:", round(x$full_model$logLik, 2), "\n")
  cat("  AIC:", round(x$full_model$AIC, 2), "\n\n")

  print(x$comparison)

  invisible(x)
}


#' Extract Item Parameters from EIRT Model
#'
#' Extract item parameters with standard errors (if available)
#'
#' @param object A gllamm_eirt object
#' @param se Logical; include standard errors?
#' @param ... Additional arguments (ignored)
#'
#' @return A data frame with item parameters
#'
#' @export
coef.gllamm_eirt <- function(object, se = FALSE, ...) {
  result <- data.frame(
    item = paste0("Item", 1:object$n_items),
    difficulty = object$item_parameters$difficulty,
    discrimination = object$item_parameters$discrimination
  )

  if (se && !inherits(object$tmb_sdr, "try-error")) {
    sdr_summary <- summary(object$tmb_sdr, "report")

    diff_se <- sdr_summary[rownames(sdr_summary) == "difficulty", "Std. Error"]
    disc_se <- sdr_summary[rownames(sdr_summary) == "discrimination", "Std. Error"]

    result$difficulty_se <- diff_se
    result$discrimination_se <- disc_se
  }

  return(result)
}


#' Predict Item Difficulties from Covariates
#'
#' Compute predicted difficulties based on item covariate model
#'
#' @param object A gllamm_eirt object
#' @param newdata Optional new item data for predictions
#'
#' @return Vector of predicted difficulties (fixed effects only, no residuals)
#'
#' @examples
#' \dontrun{
#' fit <- fit_eirt(responses, item_data,
#'                 difficulty_formula = ~ word_freq)
#'
#' # Predicted difficulties for fitted data
#' pred_diff <- predict_difficulty(fit)
#'
#' # Predicted difficulties for new items
#' new_items <- data.frame(word_freq = c(-1, 0, 1))
#' pred_new <- predict_difficulty(fit, newdata = new_items)
#' }
#'
#' @export
predict_difficulty <- function(object, newdata = NULL) {
  if (!inherits(object, "gllamm_eirt")) {
    stop("object must be of class 'gllamm_eirt'")
  }

  gamma <- object$regression_coefficients$difficulty

  if (is.null(newdata)) {
    # Use original item data
    W <- model.matrix(object$formulas$difficulty, data = object$item_data)
  } else {
    # Use new data
    W <- model.matrix(object$formulas$difficulty, data = newdata)
  }

  predicted_difficulty <- as.vector(W %*% gamma)

  return(predicted_difficulty)
}


#' Plot Item Covariate Effects
#'
#' Visualize the relationship between item covariates and item parameters
#'
#' @param object A gllamm_eirt object
#' @param covariate Name of covariate to plot
#' @param parameter Which parameter to plot: "difficulty" or "discrimination"
#' @param ... Additional arguments passed to plot
#'
#' @examples
#' \dontrun{
#' fit <- fit_eirt(responses, item_data,
#'                 difficulty_formula = ~ word_freq)
#'
#' plot_item_covariates(fit, covariate = "word_freq")
#' }
#'
#' @export
plot_item_covariates <- function(object,
                                covariate,
                                parameter = c("difficulty", "discrimination"),
                                ...) {

  parameter <- match.arg(parameter)

  if (!inherits(object, "gllamm_eirt")) {
    stop("object must be of class 'gllamm_eirt'")
  }

  if (!covariate %in% names(object$item_data)) {
    stop("Covariate '", covariate, "' not found in item_data")
  }

  x_vals <- object$item_data[[covariate]]

  if (parameter == "difficulty") {
    y_vals <- object$item_parameters$difficulty
    y_lab <- "Item Difficulty"
    formula_used <- object$formulas$difficulty
    coef_used <- object$regression_coefficients$difficulty
  } else {
    y_vals <- object$item_parameters$discrimination
    y_lab <- "Item Discrimination"
    formula_used <- object$formulas$discrimination
    coef_used <- object$regression_coefficients$discrimination
  }

  # Create plot
  plot(x_vals, y_vals,
       xlab = covariate,
       ylab = y_lab,
       main = paste(y_lab, "vs", covariate),
       pch = 19,
       col = "darkblue",
       las = 1,
       ...)

  # Add regression line if covariate is in the formula
  if (covariate %in% names(coef_used)) {
    # Get coefficient
    if ("(Intercept)" %in% names(coef_used)) {
      intercept <- coef_used[["(Intercept)"]]
      slope <- coef_used[[covariate]]

      # Predicted line
      x_seq <- seq(min(x_vals), max(x_vals), length.out = 100)
      y_pred <- intercept + slope * x_seq

      if (parameter == "discrimination") {
        y_pred <- exp(y_pred)  # Discrimination is on log scale
      }

      lines(x_seq, y_pred, col = "red", lwd = 2)

      # Add legend
      legend("topleft",
             legend = c("Observed", "Predicted"),
             pch = c(19, NA),
             lty = c(NA, 1),
             col = c("darkblue", "red"),
             lwd = c(NA, 2),
             bty = "n")
    }
  }

  # Add grid
  grid(col = "gray90")
}


#' Compute R-squared for Item Parameter Regression
#'
#' Calculate proportion of variance in item parameters explained by covariates
#'
#' @param object A gllamm_eirt object
#' @param parameter Which parameter: "difficulty" or "discrimination"
#'
#' @return R-squared value
#'
#' @export
eirt_r_squared <- function(object, parameter = c("difficulty", "discrimination")) {
  parameter <- match.arg(parameter)

  if (parameter == "difficulty") {
    observed <- object$item_parameters$difficulty
    predicted <- predict_difficulty(object)
  } else {
    # For discrimination (on log scale in model)
    observed <- log(object$item_parameters$discrimination)
    # Predict on log scale
    delta <- object$regression_coefficients$discrimination
    W <- model.matrix(object$formulas$discrimination, data = object$item_data)
    predicted <- as.vector(W %*% delta)
  }

  ss_total <- sum((observed - mean(observed))^2)
  ss_residual <- sum((observed - predicted)^2)

  r_squared <- 1 - ss_residual / ss_total

  return(r_squared)
}
