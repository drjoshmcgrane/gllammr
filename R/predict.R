#' Predict method for GLLAMM models
#'
#' Obtain predictions from a fitted GLLAMM model
#'
#' @param object A fitted \code{gllamm} object
#' @param newdata Optional new data frame for predictions. If omitted, fitted
#'   values from the original data are returned.
#' @param type Type of prediction:
#'   \describe{
#'     \item{response}{Predictions on the response scale (default)}
#'     \item{link}{Predictions on the link scale (linear predictor)}
#'     \item{random}{Random effects only}
#'   }
#' @param re.form Formula for random effects to include. Use \code{NA} or
#'   \code{~0} to exclude all random effects (population-level predictions).
#' @param ... Additional arguments (currently unused)
#'
#' @return A vector of predictions
#'
#' @examples
#' \dontrun{
#' fit <- gllamm(y ~ x + (1 | group), data = mydata)
#'
#' # Fitted values (default)
#' pred1 <- predict(fit)
#'
#' # Population-level predictions (no random effects)
#' pred2 <- predict(fit, re.form = NA)
#'
#' # Predictions for new data
#' pred3 <- predict(fit, newdata = newdata)
#' }
#'
#' @export
predict.gllamm <- function(object,
                           newdata = NULL,
                           type = c("response", "link", "random"),
                           re.form = NULL,
                           ...) {

  type <- match.arg(type)

  # If no new data, return fitted values or random effects
  if (is.null(newdata)) {
    if (type == "random") {
      return(object$random_effects)
    } else if (type == "response" || type == "link") {
      if (is.null(re.form) || (!is.na(re.form) && !identical(re.form, ~0))) {
        # Include random effects
        return(object$fitted.values)
      } else {
        # Exclude random effects - fixed effects only
        fixed_pred <- as.numeric(object$X %*% object$coefficients$fixed)
        if (type == "response" && object$family$family != "gaussian") {
          return(object$family$linkinv(fixed_pred))
        } else {
          return(fixed_pred)
        }
      }
    }
  }

  # For new data
  stop("Prediction on new data not yet implemented")
}


#' Simulate from a GLLAMM model
#'
#' Simulate response data from a fitted GLLAMM model
#'
#' @param object A fitted \code{gllamm} object
#' @param nsim Number of simulations (default: 1)
#' @param seed Optional random seed
#' @param newdata Optional new data frame. If omitted, uses original data structure.
#' @param ... Additional arguments (currently unused)
#'
#' @return If \code{nsim = 1}, a vector of simulated responses. If \code{nsim > 1},
#'   a matrix with \code{nsim} columns.
#'
#' @examples
#' \dontrun{
#' fit <- gllamm(y ~ x + (1 | group), data = mydata)
#'
#' # Single simulation
#' sim1 <- simulate(fit)
#'
#' # Multiple simulations
#' sim10 <- simulate(fit, nsim = 10)
#' }
#'
#' @export
simulate.gllamm <- function(object,
                            nsim = 1,
                            seed = NULL,
                            newdata = NULL,
                            ...) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.null(newdata)) {
    stop("Simulation with new data not yet implemented")
  }

  # Get model components
  n_obs <- object$n_obs
  X <- object$X
  beta <- object$coefficients$fixed

  # For now, only support single random effect term
  n_groups <- object$n_groups[1]
  sigma_u <- sqrt(object$coefficients$random_var[[1]])

  # Get grouping structure from original data
  parsed <- parse_formula(object$formula, object$data)
  model_data <- make_model_matrices(parsed, object$data)
  groups <- model_data$groups[[1]]
  Z <- model_data$Z[[1]]

  # Get residual SD
  # For Gaussian: sigma is stored in the model
  if (object$family$family == "gaussian") {
    # Extract from TMB fit
    sigma <- sqrt(var(object$residuals))  # Approximate for now
  } else {
    stop("Simulation for non-Gaussian families not yet implemented")
  }

  # Simulate
  sims <- matrix(NA, nrow = n_obs, ncol = nsim)

  for (s in 1:nsim) {
    # Simulate random effects
    u_sim <- rnorm(n_groups * ncol(Z), mean = 0, sd = sigma_u)

    # Compute linear predictor
    eta <- as.numeric(X %*% beta)
    for (i in 1:n_obs) {
      g <- groups[i]
      for (k in 1:ncol(Z)) {
        u_idx <- g * ncol(Z) + k
        eta[i] <- eta[i] + Z[i, k] * u_sim[u_idx]
      }
    }

    # Simulate response
    if (object$family$family == "gaussian") {
      y_sim <- rnorm(n_obs, mean = eta, sd = sigma)
    } else {
      stop("Simulation for non-Gaussian families not yet implemented")
    }

    sims[, s] <- y_sim
  }

  if (nsim == 1) {
    return(as.numeric(sims[, 1]))
  } else {
    colnames(sims) <- paste0("sim_", 1:nsim)
    return(sims)
  }
}
