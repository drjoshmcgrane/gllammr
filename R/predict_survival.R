#' Predict method for survival models
#'
#' Obtain predictions from a fitted survival model with random effects
#'
#' @param object A fitted survival model
#' @param newdata Optional new data frame for predictions
#' @param type Type of prediction:
#'   \describe{
#'     \item{lp}{Linear predictor (η)}
#'     \item{risk}{Relative risk exp(η)}
#'     \item{survival}{Survival probability at specified times}
#'     \item{hazard}{Hazard at specified times}
#'     \item{marginal_survival}{Marginal survival probability (population-averaged)}
#'     \item{marginal_hazard}{Marginal hazard (population-averaged)}
#'   }
#' @param times Times at which to evaluate survival/hazard (required for survival/hazard types)
#' @param n_sim Number of Monte Carlo samples for marginal predictions (default: 1000)
#' @param ... Additional arguments (currently unused)
#'
#' @return Depends on \code{type}
#'
#' @keywords internal
predict_survival <- function(object,
                            newdata = NULL,
                            type = c("lp", "risk", "survival", "hazard",
                                    "marginal_survival", "marginal_hazard"),
                            times = NULL,
                            n_sim = 1000,
                            ...) {
  type <- match.arg(type)

  # Extract model components
  beta <- object$coefficients$fixed
  distribution <- object$distribution  # "exponential" or "Weibull"

  # Get model matrices
  if (is.null(newdata)) {
    X <- object$X
    parsed <- parse_formula(object$formula, object$data)
    model_data <- make_model_matrices(parsed, object$data)
    Z <- model_data$Z[[1]]
  } else {
    parsed <- parse_formula(object$formula, newdata)
    new_mats <- make_model_matrices(parsed, newdata)
    X <- new_mats$X
    Z <- new_mats$Z[[1]]
  }

  n_obs <- nrow(X)

  # Linear predictor (conditional at u=0)
  if (type == "lp") {
    eta <- as.vector(X %*% beta)
    return(eta)
  }

  if (type == "risk") {
    eta <- as.vector(X %*% beta)
    return(exp(eta))
  }

  # For survival/hazard predictions, need times
  if (type %in% c("survival", "hazard", "marginal_survival", "marginal_hazard")) {
    if (is.null(times)) {
      stop("'times' must be specified for survival/hazard predictions")
    }
  }

  # Marginal predictions
  if (type %in% c("marginal_survival", "marginal_hazard")) {
    Sigma_u <- extract_random_vcov(object)
    n_random <- ncol(Sigma_u)

    # Draw samples
    if (n_random == 1) {
      sigma_u <- sqrt(Sigma_u[1, 1])
      u_samples <- matrix(rnorm(n_sim, 0, sigma_u), n_sim, 1)
    } else {
      u_samples <- rmvnorm_chol(n_sim, Sigma_u)
    }

    # Extract shape parameter for Weibull
    if (distribution == "Weibull") {
      shape <- object$shape_parameter
    }

    # Storage for marginal predictions
    n_times <- length(times)
    marginal_pred <- matrix(0, n_obs, n_times)

    # Fixed effects part
    eta_fixed <- as.vector(X %*% beta)

    # For each MC sample
    for (s in 1:n_sim) {
      u <- u_samples[s, ]

      # Random effects contribution
      if (is.matrix(Z)) {
        eta_random <- as.vector(Z %*% u)
      } else {
        eta_random <- Z * u
      }

      eta <- eta_fixed + eta_random
      lambda <- exp(eta)

      # Compute survival or hazard for each time
      for (t_idx in seq_along(times)) {
        t <- times[t_idx]

        if (type == "marginal_survival") {
          if (distribution == "exponential") {
            # S(t) = exp(-λt)
            surv <- exp(-lambda * t)
          } else {
            # S(t) = exp(-(λt)^shape)
            surv <- exp(-(lambda * t)^shape)
          }
          marginal_pred[, t_idx] <- marginal_pred[, t_idx] + surv
        } else {
          # marginal_hazard
          if (distribution == "exponential") {
            # h(t) = λ (constant hazard)
            haz <- lambda
          } else {
            # h(t) = shape * λ * (λt)^(shape-1)
            haz <- shape * lambda * (lambda * t)^(shape - 1)
          }
          marginal_pred[, t_idx] <- marginal_pred[, t_idx] + haz
        }
      }
    }

    # Average
    marginal_pred <- marginal_pred / n_sim

    if (n_times == 1) {
      return(as.vector(marginal_pred))
    } else {
      colnames(marginal_pred) <- paste0("t=", times)
      return(marginal_pred)
    }
  }

  # Conditional survival/hazard (at u=0)
  eta <- as.vector(X %*% beta)
  lambda <- exp(eta)

  if (distribution == "Weibull") {
    shape <- object$shape_parameter
  }

  n_times <- length(times)
  pred <- matrix(NA, n_obs, n_times)

  for (t_idx in seq_along(times)) {
    t <- times[t_idx]

    if (type == "survival") {
      if (distribution == "exponential") {
        pred[, t_idx] <- exp(-lambda * t)
      } else {
        pred[, t_idx] <- exp(-(lambda * t)^shape)
      }
    } else if (type == "hazard") {
      if (distribution == "exponential") {
        pred[, t_idx] <- lambda
      } else {
        pred[, t_idx] <- shape * lambda * (lambda * t)^(shape - 1)
      }
    }
  }

  if (n_times == 1) {
    return(as.vector(pred))
  } else {
    colnames(pred) <- paste0("t=", times)
    return(pred)
  }
}
