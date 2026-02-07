#' Fit Ordinal Regression Models with Random Effects
#'
#' Fit proportional odds or cumulative probit models for ordinal responses
#'
#' @param formula Formula with syntax: y ~ x + (terms | group)
#' @param data Data frame
#' @param link Link function: "logit" (proportional odds) or "probit"
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
#' @examples
#' \dontrun{
#' # Simulate ordinal data
#' data$satisfaction <- factor(sample(1:5, 100, replace = TRUE),
#'                             ordered = TRUE,
#'                             levels = 1:5)
#'
#' # Fit proportional odds model
#' fit <- fit_ordinal(satisfaction ~ age + (1 | clinic),
#'                    data = data,
#'                    link = "logit")
#' summary(fit)
#' }
#'
#' @export
fit_ordinal <- function(formula, data, link = c("logit", "probit"),
                        start = NULL, control = list()) {

  link <- match.arg(link)
  link_code <- ifelse(link == "logit", 1L, 2L)

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
    link = as.integer(link_code),
    correlated = as.integer(correlated)
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
    beta_init <- rep(0, model_data$n_fixed)
    u_init <- rep(0, tmb_data$n_groups * n_random)
    log_sigma_u_init <- rep(log(0.5), n_random)

    n_theta <- n_random * (n_random - 1) / 2
    theta_init <- rep(0, max(n_theta, 1))

    tmb_params <- list(
      beta = beta_init,
      u = u_init,
      threshold = threshold_init,
      log_sigma_u = log_sigma_u_init,
      theta = theta_init
    )
  } else {
    tmb_params <- start
  }

  # Create TMB object
  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = "u",
    DLL = "gllamm_ordinal",
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

  beta_hat <- par_full[names(par_full) == "beta"]
  names(beta_hat) <- colnames(model_data$X)

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
  result <- list(
    coefficients = list(
      fixed = beta_hat,
      thresholds = ordered_threshold,
      random_var = sigma_u_hat^2
    ),
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

  cat("Fixed effects:\n")
  print(round(x$coefficients$fixed, 3))

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
    correlated = as.integer(correlated)
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

  # Create TMB object
  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = "u",
    DLL = "gllamm_multinomial",
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
