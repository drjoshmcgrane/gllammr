#' Fit Ordinal Regression Models with Random Effects
#'
#' Fit proportional odds or cumulative probit models for ordinal responses.
#' This function can be called directly or through \code{gllamm()} with
#' \code{family = ordinal()}. The \code{gllamm()} interface is recommended
#' for consistency with other model types.
#'
#' @param formula Formula with syntax: y ~ x + (terms | group)
#' @param data Data frame
#' @param link Link function: "logit" (proportional odds) or "probit"
#' @param weights Optional vector of case weights (one per observation)
#' @param start Optional starting values
#' @param control Control parameters
#'
#' @return An object of class \code{gllamm_ordinal}
#'
#' @details
#' For ordinal response Y with K categories (1, 2, ..., K), the model is:
#'
#' Proportional odds (logit link):
#' \deqn{P(Y \le k | x) = \frac{1}{1 + \exp(-(\tau_k - x'\beta))}}
#'
#' Cumulative probit:
#' \deqn{P(Y \le k | x) = \Phi(\tau_k - x'\beta)}
#'
#' where \eqn{\tau_1 < \tau_2 < ... < \tau_{K-1}} are threshold parameters.
#'
#' @note
#' The recommended interface is \code{gllamm(formula, data, family = ordinal(link))}.
#' This function is also available for direct use with the \code{link} argument.
#'
#' @examples
#' \dontrun{
#' # Simulate ordinal data
#' data$satisfaction <- factor(sample(1:5, 100, replace = TRUE),
#'                             ordered = TRUE,
#'                             levels = 1:5)
#'
#' # Recommended: Use gllamm() with ordinal() family
#' fit1 <- gllamm(satisfaction ~ age + (1 | clinic),
#'                data = data,
#'                family = ordinal(link = "logit"))
#' summary(fit1)
#'
#' # Alternative: Call fit_ordinal() directly
#' fit2 <- fit_ordinal(satisfaction ~ age + (1 | clinic),
#'                     data = data,
#'                     link = "logit")
#' summary(fit2)
#' }
#'
#' @export
fit_ordinal <- function(formula, data, link = c("logit", "probit", "acl",
                                                "crl_forward", "crl_backward", "ppo"),
                        weights = NULL,
                        start = NULL, control = list()) {

  link <- match.arg(link)

  # Map link function to numeric code for TMB
  link_code <- switch(link,
    logit = 1L,
    probit = 2L,
    acl = 3L,
    crl_forward = 4L,
    crl_backward = 5L,
    ppo = 6L
  )

  # Validate weights if provided
  if (!is.null(weights)) {
    if (length(weights) != nrow(data)) {
      stop("Length of weights (", length(weights), ") must match number of observations (", nrow(data), ")")
    }
    if (any(weights < 0, na.rm = TRUE)) {
      stop("All weights must be non-negative")
    }
    if (any(is.na(weights))) {
      stop("weights cannot contain missing values")
    }
  }

  # Parse formula
  parsed <- parse_formula(formula, data)
  model_data <- make_model_matrices(parsed, data)

  # Get response as factor
  y_response <- model_data$y

  # Convert to numeric if factor/ordered
  if (is.factor(y_response)) {
    y_numeric <- as.integer(y_response)
    n_categories <- nlevels(y_response)
    category_labels <- levels(y_response)
  } else {
    y_numeric <- as.integer(y_response)
    n_categories <- length(unique(y_numeric))
    category_labels <- sort(unique(y_numeric))
  }

  # Check that responses are in 1:K
  if (min(y_numeric) != 1) {
    stop("Ordinal responses must start at 1")
  }

  # Only handle single random effect term for now
  if (model_data$n_random_terms != 1) {
    stop("Currently only single random effects term supported")
  }

  n_random <- model_data$n_random_coefs[1]
  correlated <- !parsed$random_terms[[1]]$uncorrelated

  # Convert Z to sparse matrix
  Z_sparse <- Matrix::Matrix(model_data$Z[[1]], sparse = TRUE)

  # Prepare weights vector (default to 1.0 if NULL)
  if (is.null(weights)) {
    weights_vec <- rep(1.0, model_data$n_obs)
  } else {
    weights_vec <- as.numeric(weights)
  }

  # Cumulative-link models: the thresholds carry the location, so the
  # intercept must be dropped from X (a free intercept and free thresholds
  # are jointly unidentified)
  X_fixed <- drop_intercept_column(as.matrix(model_data$X))
  n_fixed <- ncol(X_fixed)

  # Prepare TMB data
  tmb_data <- list(
    y = as.integer(y_numeric),
    X = X_fixed,
    Z = Z_sparse,
    groups = as.integer(model_data$groups[[1]]),
    n_groups = as.integer(model_data$n_groups[1]),
    n_obs = as.integer(model_data$n_obs),
    n_fixed = as.integer(n_fixed),
    n_random = as.integer(n_random),
    n_categories = as.integer(n_categories),
    link = as.integer(link_code),
    correlated = as.integer(correlated),
    weights = weights_vec
  )

  # Initialize parameters
  if (is.null(start)) {
    # Initialize thresholds evenly spaced
    threshold_init <- seq(-1, 1, length.out = n_categories - 1)
    # Transform to ensure ordering (using differences)
    threshold_init[1] <- threshold_init[1]
    for (k in 2:length(threshold_init)) {
      threshold_init[k] <- log(threshold_init[k] - threshold_init[k-1])
    }

    # Initialize other parameters
    u_init <- rep(0, tmb_data$n_groups * n_random)
    log_sigma_u_init <- rep(log(0.5), n_random)

    n_theta <- n_random * (n_random - 1) / 2
    theta_init <- rep(0, max(n_theta, 1))

    # PPO requires beta_ppo matrix, others use beta vector
    if (link == "ppo") {
      beta_ppo_init <- matrix(0, nrow = n_categories - 1, ncol = n_fixed)
      beta_init <- rep(0, n_fixed)  # Still need beta for compatibility

      tmb_params <- list(
        beta = beta_init,
        u = u_init,
        threshold = threshold_init,
        log_sigma_u = log_sigma_u_init,
        theta = theta_init,
        beta_ppo = beta_ppo_init
      )
    } else {
      beta_init <- rep(0, n_fixed)
      # For non-PPO models, beta_ppo is still needed but not used
      beta_ppo_init <- matrix(0, nrow = n_categories - 1, ncol = n_fixed)

      tmb_params <- list(
        beta = beta_init,
        u = u_init,
        threshold = threshold_init,
        log_sigma_u = log_sigma_u_init,
        theta = theta_init,
        beta_ppo = beta_ppo_init
      )
    }
  } else {
    tmb_params <- start
  }

  # Fix parameters the chosen link never reads (dead parameters create flat
  # likelihood directions): beta_ppo for links 1-5, beta for PPO, theta
  # without a correlation structure
  tmb_map <- list()
  if (link == "ppo") {
    if (length(tmb_params$beta) > 0) {
      tmb_map$beta <- factor(rep(NA, length(tmb_params$beta)))
    }
  } else {
    tmb_map$beta_ppo <- factor(rep(NA, length(tmb_params$beta_ppo)))
  }
  if (!(correlated && n_random > 1)) {
    tmb_map$theta <- factor(rep(NA, length(tmb_params$theta)))
  }

  # Create TMB object
  tmb_data$model_name <- "ordinal"
  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = "u",
    map = tmb_map,
    DLL = "GLLAMMR",
    silent = TRUE
  )

  # Optimize. Box bounds keep the model finite under separation (a
  # deterministic predictor drives beta and threshold spacings to infinity;
  # bounded estimates with huge SEs match clmm-style behavior).
  control_defaults <- list(eval.max = 2000, iter.max = 1000, trace = 0)
  control <- modifyList(control_defaults, control)

  par_names_opt <- names(obj$par)
  lower <- rep(-30, length(par_names_opt))
  upper <- rep(30, length(par_names_opt))
  lower[par_names_opt == "threshold"] <- -15
  upper[par_names_opt == "threshold"] <- 15   # log-spacing beyond this overflows
  lower[par_names_opt == "log_sigma_u"] <- -10
  upper[par_names_opt == "log_sigma_u"] <- 10

  opt <- nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    lower = lower,
    upper = upper,
    control = control
  )

  if (any(abs(opt$par[par_names_opt == "beta"]) >= 29.5)) {
    warning("Some coefficients reached the optimization boundary; the data ",
            "may exhibit complete separation.")
  }

  # Get standard errors
  sdr <- try(TMB::sdreport(obj), silent = TRUE)

  # Extract parameters
  par_full <- obj$env$last.par.best

  # Extract beta (for non-PPO) or beta_ppo (for PPO)
  if (link == "ppo") {
    beta_ppo_hat <- matrix(par_full[names(par_full) == "beta_ppo"],
                           nrow = n_categories - 1,
                           ncol = n_fixed)
    rownames(beta_ppo_hat) <- paste0("Threshold", 1:(n_categories - 1))
    colnames(beta_ppo_hat) <- colnames(X_fixed)
    beta_hat <- NULL  # Not used for PPO
  } else {
    beta_hat <- par_full[names(par_full) == "beta"]
    names(beta_hat) <- colnames(X_fixed)
    beta_ppo_hat <- NULL
  }

  threshold_hat <- par_full[names(par_full) == "threshold"]

  # Reconstruct ordered thresholds
  ordered_threshold <- numeric(n_categories - 1)
  ordered_threshold[1] <- threshold_hat[1]
  for (k in 2:(n_categories - 1)) {
    ordered_threshold[k] <- ordered_threshold[k-1] + exp(threshold_hat[k])
  }
  names(ordered_threshold) <- paste0("tau", 1:(n_categories - 1))

  log_sigma_u_hat <- par_full[names(par_full) == "log_sigma_u"]
  sigma_u_hat <- exp(log_sigma_u_hat)

  # Random effects
  u_hat <- par_full[names(par_full) == "u"]
  random_effects <- list()
  for (g in 1:tmb_data$n_groups) {
    idx_start <- (g - 1) * n_random + 1
    idx_end <- g * n_random
    random_effects[[g]] <- u_hat[idx_start:idx_end]
  }

  # Construct result
  coef_list <- list(
    thresholds = ordered_threshold,
    random_var = sigma_u_hat^2
  )

  # Add fixed effects (beta for non-PPO, beta_ppo for PPO)
  if (link == "ppo") {
    coef_list$beta_ppo <- beta_ppo_hat
  } else {
    coef_list$fixed <- beta_hat
  }

  result <- list(
    coefficients = coef_list,
    random_effects = random_effects,
    link = link,
    n_categories = n_categories,
    category_labels = category_labels,
    logLik = -opt$objective,
    AIC = 2 * opt$objective + 2 * length(obj$par),
    BIC = 2 * opt$objective + log(model_data$n_obs) * length(obj$par),
    convergence = list(
      converged = (opt$convergence == 0),
      message = opt$message
    ),
    n_obs = model_data$n_obs,
    formula = formula,
    data = data,
    X = model_data$X,
    tmb_obj = obj,
    tmb_opt = opt,
    tmb_sdr = sdr
  )

  class(result) <- c("gllamm_ordinal", "gllamm")

  return(result)
}


#' @export
print.gllamm_ordinal <- function(x, ...) {
  cat("Ordinal Regression Model\n")
  cat("Link:", x$link, "\n\n")

  cat("Number of observations:", x$n_obs, "\n")
  cat("Number of categories:", x$n_categories, "\n\n")

  if (x$link == "ppo") {
    cat("Partial Proportional Odds Coefficients (by threshold):\n")
    print(round(x$coefficients$beta_ppo, 3))
  } else {
    cat("Fixed effects:\n")
    print(round(x$coefficients$fixed, 3))
  }

  cat("\nThreshold parameters:\n")
  print(round(x$coefficients$thresholds, 3))

  cat("\nRandom effects variance:", round(x$coefficients$random_var, 3), "\n")

  cat("\nLog-likelihood:", round(x$logLik, 2), "\n")
  cat("AIC:", round(x$AIC, 2), "\n")

  invisible(x)
}


#' @export
summary.gllamm_ordinal <- function(object, ...) {
  print(object)

  cat("\nCategory distribution:\n")
  # Would need original data to show this
  cat("(Use table(data$response) to see category frequencies)\n")

  invisible(object)
}


#' Fit Multinomial Regression Models with Random Effects
#'
#' Fit baseline category logit models for nominal (unordered) responses
#'
#' @param formula Formula with syntax: y ~ x + (terms | group)
#' @param data Data frame
#' @param reference Reference category (default: first level)
#' @param start Optional starting values
#' @param control Control parameters
#'
#' @return An object of class \code{gllamm_multinomial}
#'
#' @details
#' For nominal response Y with K categories, using baseline category logit:
#'
#' \deqn{P(Y = k | x) = \frac{\exp(x'\beta_k)}{1 + \sum_{j=1}^{K-1} \exp(x'\beta_j)}}
#'
#' where category 0 is the reference with \eqn{\beta_0 = 0}.
#'
#' @examples
#' \dontrun{
#' # Simulate multinomial data
#' data$choice <- factor(sample(c("A", "B", "C"), 100, replace = TRUE))
#'
#' # Fit multinomial model
#' fit <- fit_multinomial(choice ~ price + quality + (1 | person),
#'                        data = data)
#' summary(fit)
#' }
#'
#' @export
fit_multinomial <- function(formula, data, reference = NULL,
                            start = NULL, control = list()) {

  # Parse formula
  parsed <- parse_formula(formula, data)
  model_data <- make_model_matrices(parsed, data)

  # Get response as factor
  y_response <- model_data$y

  if (!is.factor(y_response)) {
    y_response <- factor(y_response)
  }

  category_labels <- levels(y_response)
  n_categories <- length(category_labels)

  # Set reference category
  if (is.null(reference)) {
    reference <- category_labels[1]
  }

  # Recode so reference is category 0
  y_numeric <- as.integer(y_response) - 1

  # Only handle single random effect term for now
  if (model_data$n_random_terms != 1) {
    stop("Currently only single random effects term supported")
  }

  n_random <- model_data$n_random_coefs[1]
  correlated <- !parsed$random_terms[[1]]$uncorrelated

  Z_sparse <- Matrix::Matrix(model_data$Z[[1]], sparse = TRUE)

  # Prepare TMB data
  tmb_data <- list(
    y = as.integer(y_numeric),
    X = as.matrix(model_data$X),
    Z = Z_sparse,
    groups = as.integer(model_data$groups[[1]]),
    n_groups = as.integer(model_data$n_groups[1]),
    n_obs = as.integer(model_data$n_obs),
    n_fixed = as.integer(model_data$n_fixed),
    n_random = as.integer(n_random),
    n_categories = as.integer(n_categories),
    correlated = as.integer(correlated),
    weights = rep(1.0, model_data$n_obs)
  )

  # Initialize parameters
  if (is.null(start)) {
    # Beta is a matrix: (n_categories - 1) x n_fixed
    beta_init <- matrix(0, n_categories - 1, model_data$n_fixed)

    u_init <- rep(0, tmb_data$n_groups * n_random)
    log_sigma_u_init <- rep(log(0.5), n_random)

    n_theta <- n_random * (n_random - 1) / 2
    theta_init <- rep(0, max(n_theta, 1))

    tmb_params <- list(
      beta = beta_init,
      u = u_init,
      log_sigma_u = log_sigma_u_init,
      theta = theta_init
    )
  } else {
    tmb_params <- start
  }

  # Fix theta when the template never reads it (no correlation structure)
  tmb_map <- list()
  if (!(correlated && n_random > 1)) {
    tmb_map$theta <- factor(rep(NA, length(tmb_params$theta)))
  }

  # Create TMB object
  tmb_data$model_name <- "multinomial"
  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = "u",
    map = tmb_map,
    DLL = "GLLAMMR",
    silent = TRUE
  )

  # Optimize
  control_defaults <- list(eval.max = 2000, iter.max = 1000, trace = 0)
  control <- modifyList(control_defaults, control)

  opt <- nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = control
  )

  # Get standard errors
  sdr <- try(TMB::sdreport(obj), silent = TRUE)

  # Extract parameters
  par_full <- obj$env$last.par.best

  beta_hat <- matrix(par_full[names(par_full) == "beta"],
                     nrow = n_categories - 1,
                     ncol = model_data$n_fixed)
  rownames(beta_hat) <- category_labels[-1]  # Exclude reference
  colnames(beta_hat) <- colnames(model_data$X)

  # Random effects variance
  log_sigma_u_hat <- par_full[names(par_full) == "log_sigma_u"]
  sigma_u_hat <- exp(log_sigma_u_hat)

  # Construct result
  result <- list(
    coefficients = list(
      beta = beta_hat,
      random_var = sigma_u_hat^2
    ),
    reference = reference,
    categories = category_labels,
    n_categories = n_categories,
    logLik = -opt$objective,
    AIC = 2 * opt$objective + 2 * length(obj$par),
    BIC = 2 * opt$objective + log(model_data$n_obs) * length(obj$par),
    convergence = list(
      converged = (opt$convergence == 0),
      message = opt$message
    ),
    n_obs = model_data$n_obs,
    formula = formula,
    data = data,
    X = as.matrix(model_data$X),
    tmb_obj = obj,
    tmb_opt = opt,
    tmb_sdr = sdr
  )

  class(result) <- c("gllamm_multinomial", "gllamm")

  return(result)
}


#' @export
print.gllamm_multinomial <- function(x, ...) {
  cat("Multinomial Regression Model\n")
  cat("Reference category:", x$reference, "\n\n")

  cat("Number of observations:", x$n_obs, "\n")
  cat("Number of categories:", x$n_categories, "\n\n")

  cat("Coefficients (reference =", x$reference, "):\n")
  print(round(x$coefficients$beta, 3))

  cat("\nRandom effects variance:", round(x$coefficients$random_var, 3), "\n")

  cat("\nLog-likelihood:", round(x$logLik, 2), "\n")
  cat("AIC:", round(x$AIC, 2), "\n")

  invisible(x)
}


#' @export
summary.gllamm_multinomial <- function(object, ...) {
  print(object)
  invisible(object)
}


#' Test Proportional Odds Assumption
#'
#' Perform a likelihood ratio test of the proportional odds assumption
#' by comparing a proportional odds model to a partial proportional odds model
#'
#' @param object A fitted ordinal regression model (gllamm_ordinal)
#'
#' @return An object of class \code{po_test} with components:
#'   \item{statistic}{Likelihood ratio test statistic}
#'   \item{df}{Degrees of freedom for the test}
#'   \item{p_value}{P-value from chi-squared distribution}
#'   \item{conclusion}{Text interpretation of the test result}
#'   \item{models}{List containing the base and PPO models}
#'
#' @details
#' The proportional odds assumption states that the effect of covariates
#' is the same across all thresholds. This function fits a partial proportional
#' odds (PPO) model where each threshold can have different covariate effects
#' and tests whether this provides a significantly better fit.
#'
#' The test statistic is:
#' \deqn{LRT = 2(logLik_{PPO} - logLik_{PO})}
#'
#' which follows a chi-squared distribution with degrees of freedom equal to
#' the difference in number of parameters.
#'
#' If p < 0.05, the proportional odds assumption is rejected, suggesting
#' that covariate effects vary across thresholds.
#'
#' @note This function currently only works for models with logit or probit links.
#' It will not work with ACL, CRL, or already-PPO models.
#'
#' @examples
#' \dontrun{
#' # Fit proportional odds model
#' fit_po <- fit_ordinal(rating ~ temp + (1 | judge),
#'                       data = wine, link = "logit")
#'
#' # Test proportional odds assumption
#' po_test <- test_proportional_odds(fit_po)
#' print(po_test)
#' }
#'
#' @param data The data frame used to fit the model (required when the
#'   fitted object does not store its data)
#' @export
test_proportional_odds <- function(object, data = NULL) {

  # Check that input is an ordinal model
  if (!inherits(object, "gllamm_ordinal")) {
    stop("Object must be of class 'gllamm_ordinal'")
  }

  # Check that link is logit or probit (can be relaxed to PPO)
  if (!object$link %in% c("logit", "probit")) {
    stop("Proportional odds test only applies to logit or probit link models")
  }

  # Need data to refit model
  if (is.null(data)) {
    stop("Data must be provided to test proportional odds assumption.\n",
         "Usage: test_proportional_odds(model, data = your_data)")
  }

  # Extract original formula
  formula <- object$formula

  # Refit as PPO model
  cat("Fitting partial proportional odds model...\n")
  ppo_fit <- tryCatch({
    fit_ordinal(formula, data, link = "ppo",
                control = list(trace = 0))
  }, error = function(e) {
    stop("Failed to fit PPO model: ", e$message,
         "\nThis may indicate convergence issues or data problems.")
  })

  # Extract information for test
  logLik_base <- object$logLik
  logLik_ppo <- ppo_fit$logLik

  n_params_base <- length(object$tmb_obj$par)
  n_params_ppo <- length(ppo_fit$tmb_obj$par)

  # Likelihood ratio test
  lrt_stat <- 2 * (logLik_ppo - logLik_base)
  df <- n_params_ppo - n_params_base
  p_value <- pchisq(lrt_stat, df, lower.tail = FALSE)

  # Conclusion
  if (p_value < 0.01) {
    conclusion <- "Strong evidence against proportional odds (p < 0.01). Use PPO model."
  } else if (p_value < 0.05) {
    conclusion <- "Reject proportional odds assumption (p < 0.05). Consider PPO model."
  } else if (p_value < 0.10) {
    conclusion <- "Weak evidence against proportional odds (p < 0.10). PO may be adequate."
  } else {
    conclusion <- "Proportional odds assumption is reasonable (p >= 0.10)."
  }

  result <- structure(
    list(
      statistic = lrt_stat,
      df = df,
      p_value = p_value,
      conclusion = conclusion,
      base_logLik = logLik_base,
      ppo_logLik = logLik_ppo,
      base_model = object,
      ppo_model = ppo_fit
    ),
    class = "po_test"
  )

  return(result)
}


#' Print method for proportional odds test
#' @keywords internal
#' @export
print.po_test <- function(x, ...) {
  cat("Proportional Odds Assumption Test\n")
  cat("==================================\n\n")

  cat("Likelihood Ratio Test:\n")
  cat("  LRT statistic:", round(x$statistic, 3), "\n")
  cat("  Degrees of freedom:", x$df, "\n")
  cat("  P-value:", format.pval(x$p_value), "\n\n")

  cat("Conclusion:\n")
  cat(" ", x$conclusion, "\n\n")

  cat("Model comparison:\n")
  cat("  PO logLik:", round(x$base_logLik, 2), "\n")
  cat("  PPO logLik:", round(x$ppo_logLik, 2), "\n")

  invisible(x)
}
