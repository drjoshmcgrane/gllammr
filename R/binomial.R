#' Fit Binomial Regression Models with Random Effects
#'
#' Fit logistic, probit, or complementary log-log models for binary or binomial
#' responses. This function can be called directly or through \code{gllamm()} with
#' \code{family = binomial()}. The \code{gllamm()} interface is recommended
#' for consistency with other model types.
#'
#' @param formula Formula with syntax: y ~ x + (terms | group)
#' @param data Data frame
#' @param link Link function: "logit" (default), "probit", or "cloglog"
#' @param start Optional starting values
#' @param control Control parameters for optimization
#'
#' @return An object of class \code{gllamm_binomial} with components:
#'   \item{coefficients}{List of fixed effects and random effect variances}
#'   \item{logLik}{Log-likelihood at convergence}
#'   \item{AIC}{Akaike Information Criterion}
#'   \item{BIC}{Bayesian Information Criterion}
#'   \item{convergence}{Convergence information}
#'   \item{link}{Link function used}
#'   \item{n_obs}{Number of observations}
#'   \item{fitted_values}{Fitted probabilities}
#'   \item{tmb_obj}{TMB object for further inference}
#'   \item{tmb_opt}{Optimization result}
#'   \item{tmb_sdr}{Standard errors via sdreport}
#'
#' @details
#' For binary response Y in {0, 1} or binomial response Y/n, the model is:
#'
#' \strong{Logit link (default - logistic regression):}
#' \deqn{P(Y=1|x) = \frac{1}{1 + \exp(-x'\beta - Z'u)}}
#'
#' \strong{Probit link:}
#' \deqn{P(Y=1|x) = \Phi(x'\beta + Z'u)}
#' where \eqn{\Phi} is the standard normal CDF.
#'
#' \strong{Complementary log-log (cloglog) link:}
#' \deqn{P(Y=1|x) = 1 - \exp(-\exp(x'\beta + Z'u))}
#'
#' The cloglog link is asymmetric and particularly useful for:
#' \itemize{
#'   \item Rare events (when baseline probability is low)
#'   \item Survival analysis with discrete time
#'   \item When modeling hazards that are proportional
#'   \item Grouped survival data from an underlying Poisson process
#' }
#'
#' Random effects u are assumed multivariate normal with mean 0.
#'
#' @note
#' The recommended interface is \code{gllamm(formula, data, family = binomial(link))}.
#' This function is also available for direct use with the \code{link} argument.
#'
#' @examples
#' \dontrun{
#' # Simulate binary data
#' set.seed(123)
#' n_groups <- 20
#' n_per_group <- 10
#' data <- data.frame(
#'   group = rep(1:n_groups, each = n_per_group),
#'   x = rnorm(n_groups * n_per_group),
#'   y = rbinom(n_groups * n_per_group, 1, 0.5)
#' )
#'
#' # Recommended: Use gllamm() with binomial() family
#' fit1 <- gllamm(y ~ x + (1 | group),
#'                data = data,
#'                family = binomial(link = "logit"))
#' summary(fit1)
#'
#' # Probit link
#' fit2 <- gllamm(y ~ x + (1 | group),
#'                data = data,
#'                family = binomial(link = "probit"))
#'
#' # Complementary log-log for rare events
#' data$rare_event <- rbinom(nrow(data), 1, 0.05)
#' fit3 <- gllamm(rare_event ~ x + (1 | group),
#'                data = data,
#'                family = binomial(link = "cloglog"))
#' summary(fit3)
#'
#' # Alternative: Call fit_binomial() directly
#' fit4 <- fit_binomial(y ~ x + (1 | group),
#'                      data = data,
#'                      link = "logit")
#' }
#'
#' @seealso \code{\link{gllamm}}, \code{\link{binomial}}, \code{\link{ordinal}}
#'
#' @export
fit_binomial <- function(formula, data, link = c("logit", "probit", "cloglog"),
                         weights = NULL,
                         start = NULL, control = list()) {

  link <- match.arg(link)

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

  # Map link to numeric code for TMB
  link_code <- switch(link,
    logit = 1L,
    probit = 2L,
    cloglog = 3L
  )

  # Parse formula
  parsed <- parse_formula(formula, data)
  model_data <- make_model_matrices(parsed, data)

  # Get response
  y_response <- model_data$y

  # Check binary (0/1)
  if (!all(y_response %in% c(0, 1, NA))) {
    stop("Response must be binary (0/1). For binomial counts, use appropriate syntax.")
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

  # Prepare TMB data
  tmb_data <- list(
    y = as.numeric(y_response),
    X = as.matrix(model_data$X),
    Z = Z_sparse,
    groups = as.integer(model_data$groups[[1]]),
    n_groups = as.integer(model_data$n_groups[1]),
    n_obs = as.integer(model_data$n_obs),
    n_fixed = as.integer(model_data$n_fixed),
    n_random = as.integer(n_random),
    link = as.integer(link_code),
    correlated = as.integer(correlated),
    weights = weights_vec,
    group_weights = rep(1.0, model_data$n_groups[1]),
    model_name = "binomial"
  )

  # Initialize parameters
  if (is.null(start)) {
    beta_init <- rep(0, model_data$n_fixed)
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

  # Fix theta when the template never reads it (no correlation structure);
  # leaving it free creates a flat likelihood direction
  tmb_map <- list()
  if (!(correlated && n_random > 1)) {
    tmb_map$theta <- factor(rep(NA, length(tmb_params$theta)))
  }

  # Create TMB object
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

  # Fixed effects
  beta_hat <- par_full[names(par_full) == "beta"]
  names(beta_hat) <- colnames(model_data$X)

  # Random effects standard deviations
  log_sigma_u_hat <- par_full[names(par_full) == "log_sigma_u"]
  sigma_u_hat <- exp(log_sigma_u_hat)
  names(sigma_u_hat) <- paste0("sigma_", seq_along(sigma_u_hat))

  # Full random-effects covariance (reconstruct normalized correlation from
  # the Cholesky parameters when present, matching the template)
  if (correlated && n_random > 1) {
    theta_hat <- par_full[names(par_full) == "theta"]
    L <- diag(n_random)
    idx <- 1
    for (i in 2:n_random) {       # row-major fill, matching the template
      for (j in 1:(i - 1)) {
        L[i, j] <- theta_hat[idx]
        idx <- idx + 1
      }
    }
    R <- L %*% t(L)
    R <- R / sqrt(diag(R) %o% diag(R))
    Sigma_u_hat <- outer(sigma_u_hat, sigma_u_hat) * R
  } else {
    Sigma_u_hat <- diag(sigma_u_hat^2, nrow = n_random)
  }

  # Get fitted values
  fitted_vals <- obj$report()$fitted

  # Construct result
  result <- list(
    coefficients = list(
      fixed = beta_hat,
      random_sd = sigma_u_hat,
      random_var = list(Sigma_u_hat)
    ),
    logLik = -opt$objective,
    AIC = 2 * opt$objective + 2 * length(obj$par),
    BIC = 2 * opt$objective + log(model_data$n_obs) * length(obj$par),
    convergence = list(
      converged = (opt$convergence == 0),
      message = opt$message,
      iterations = opt$iterations
    ),
    link = link,
    family = stats::binomial(link = link),
    n_obs = model_data$n_obs,
    n_groups = model_data$n_groups[1],
    fitted_values = fitted_vals,
    residuals = y_response - fitted_vals,
    formula = formula,
    data = data,
    tmb_obj = obj,
    tmb_opt = opt,
    tmb_sdr = sdr
  )

  class(result) <- c("gllamm_binomial", "gllamm")

  return(result)
}


#' @export
print.gllamm_binomial <- function(x, ...) {
  cat("Binomial Regression Model with Random Effects\n")
  cat("Link function:", x$link, "\n\n")

  cat("Number of observations:", x$n_obs, "\n")
  cat("Number of groups:", x$n_groups, "\n\n")

  cat("Fixed effects:\n")
  print(round(x$coefficients$fixed, 4))

  cat("\nRandom effects standard deviations:\n")
  print(round(x$coefficients$random_sd, 4))

  cat("\nLog-likelihood:", round(x$logLik, 2), "\n")
  cat("AIC:", round(x$AIC, 2), "\n")
  cat("BIC:", round(x$BIC, 2), "\n")

  if (x$convergence$converged) {
    cat("\nConvergence: successful\n")
  } else {
    cat("\nConvergence: FAILED -", x$convergence$message, "\n")
  }

  invisible(x)
}


#' @export
summary.gllamm_binomial <- function(object, ...) {
  cat("Binomial Regression Model with Random Effects\n")
  cat("==============================================\n\n")

  cat("Link function:", object$link, "\n")
  cat("Number of observations:", object$n_obs, "\n")
  cat("Number of groups:", object$n_groups, "\n\n")

  cat("Fixed Effects:\n")
  if (!inherits(object$tmb_sdr, "try-error")) {
    sdr_summary <- summary(object$tmb_sdr, "fixed")
    beta_summary <- sdr_summary[rownames(sdr_summary) == "beta", ]
    beta_df <- data.frame(
      Estimate = beta_summary[, "Estimate"],
      Std.Error = beta_summary[, "Std. Error"],
      z.value = beta_summary[, "Estimate"] / beta_summary[, "Std. Error"],
      p.value = 2 * pnorm(-abs(beta_summary[, "Estimate"] / beta_summary[, "Std. Error"]))
    )
    rownames(beta_df) <- names(object$coefficients$fixed)
    print(beta_df)
  } else {
    print(object$coefficients$fixed)
    cat("\n(Standard errors not available)\n")
  }

  cat("\nRandom Effects:\n")
  cat("  Standard deviations:\n")
  for (i in seq_along(object$coefficients$random_sd)) {
    cat(sprintf("    %s: %.4f\n",
                names(object$coefficients$random_sd)[i],
                object$coefficients$random_sd[i]))
  }

  cat("\nModel Fit:\n")
  cat("  Log-likelihood:", round(object$logLik, 2), "\n")
  cat("  AIC:", round(object$AIC, 2), "\n")
  cat("  BIC:", round(object$BIC, 2), "\n")

  cat("\nConvergence:", ifelse(object$convergence$converged,
                               "successful",
                               paste("FAILED -", object$convergence$message)), "\n")

  invisible(object)
}
