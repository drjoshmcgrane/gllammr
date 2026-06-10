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
#'     \item{marginal}{Population-averaged predictions (integrating over random effects)}
#'   }
#' @param re.form Formula for random effects to include. Use \code{NA} or
#'   \code{~0} to exclude all random effects (population-level predictions).
#'   Ignored when \code{type = "marginal"}.
#' @param n_sim Number of Monte Carlo samples for marginal predictions (default: 1000).
#'   Only used when \code{type = "marginal"}.
#' @param se.fit Logical; return standard errors for marginal predictions? (default: FALSE).
#'   Only used when \code{type = "marginal"}.
#' @param ... Additional arguments (currently unused)
#'
#' @return A vector of predictions
#'
#' @examples
#' \dontrun{
#' fit <- gllamm(y ~ x + (1 | group), data = mydata)
#'
#' # Fitted values (default - conditional on random effects)
#' pred1 <- predict(fit)
#'
#' # Population-level predictions (fixed effects only, u=0)
#' pred2 <- predict(fit, re.form = NA)
#'
#' # Marginal predictions (population-averaged, integrating over u)
#' pred3 <- predict(fit, type = "marginal")
#'
#' # Marginal predictions with standard errors
#' pred4 <- predict(fit, type = "marginal", se.fit = TRUE)
#'
#' # Marginal predictions for new data
#' pred5 <- predict(fit, newdata = newdata, type = "marginal")
#' }
#'
#' @export
predict.gllamm <- function(object,
                           newdata = NULL,
                           type = c("response", "link", "random", "marginal"),
                           re.form = NULL,
                           n_sim = 1000,
                           se.fit = FALSE,
                           ...) {

  type <- match.arg(type)

  # Handle marginal predictions separately
  if (type == "marginal") {
    return(predict_marginal_gllamm(object, newdata, n_sim, se.fit))
  }

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

  # For new data - implement prediction
  # Parse formula to get structure
  parsed <- parse_formula(object$formula, newdata)
  new_mats <- make_model_matrices(parsed, newdata)

  # Fixed effects predictions
  fixed_pred <- as.numeric(new_mats$X %*% object$coefficients$fixed)

  if (is.null(re.form) || (!is.na(re.form) && !identical(re.form, ~0))) {
    # Include random effects - use posterior means (BLUPs) from training
    # For new groups not in training, use 0 (population average)

    # Get group identifiers from new data
    if (length(object$random_terms) > 0) {
      rt <- object$random_terms[[1]]
      new_groups <- if (rt$nested) {
        interaction(newdata[, rt$grouping], drop = TRUE)
      } else {
        factor(newdata[[rt$grouping]])
      }

      # Match to training groups
      training_groups <- names(object$random_effects) %||%
                        paste0("Group", seq_along(object$random_effects))

      random_contrib <- numeric(nrow(newdata))
      for (i in seq_along(new_groups)) {
        g_name <- as.character(new_groups[i])
        if (g_name %in% training_groups) {
          g_idx <- which(training_groups == g_name)
          random_contrib[i] <- sum(new_mats$Z[[1]][i, ] * object$random_effects[[g_idx]])
        }
        # else: use 0 (population average) for new groups
      }

      pred <- fixed_pred + random_contrib
    } else {
      pred <- fixed_pred
    }
  } else {
    pred <- fixed_pred
  }

  # Apply inverse link if needed
  if (type == "response" && object$family$family != "gaussian") {
    pred <- object$family$linkinv(pred)
  }

  return(pred)
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


#' Internal function for marginal predictions
#'
#' Computes population-averaged predictions by integrating over random effects
#'
#' @param object Fitted gllamm object
#' @param newdata New data for predictions (NULL = use original data)
#' @param n_sim Number of Monte Carlo samples
#' @param se.fit Return standard errors?
#'
#' @return Vector of marginal predictions, or list with fit and se.fit
#' @keywords internal
predict_marginal_gllamm <- function(object, newdata = NULL, n_sim = 1000, se.fit = FALSE) {

  # Special case: Gaussian with identity link
  # Marginal = conditional for linear models
  if (object$family$family == "gaussian" && object$family$link == "identity") {
    # Just return fixed effects predictions
    if (is.null(newdata)) {
      pred <- as.numeric(object$X %*% object$coefficients$fixed)
    } else {
      parsed <- parse_formula(object$formula, newdata)
      new_mats <- make_model_matrices(parsed, newdata)
      pred <- as.numeric(new_mats$X %*% object$coefficients$fixed)
    }

    if (se.fit) {
      # SE is 0 for marginal = conditional case
      return(list(fit = pred, se.fit = rep(0, length(pred))))
    } else {
      return(pred)
    }
  }

  # General case: Nonlinear link requires Monte Carlo integration

  # Get model matrices
  if (is.null(newdata)) {
    X <- object$X
    # Construct Z for original data
    parsed <- parse_formula(object$formula, object$data)
    model_data <- make_model_matrices(parsed, object$data)
    Z <- model_data$Z[[1]]
  } else {
    parsed <- parse_formula(object$formula, newdata)
    new_mats <- make_model_matrices(parsed, newdata)
    X <- new_mats$X
    Z <- new_mats$Z[[1]]
  }

  # Get fixed effects
  beta <- object$coefficients$fixed

  # Extract random effects variance-covariance matrix
  Sigma_u <- extract_random_vcov(object)

  # Get inverse link function
  inv_link <- get_inverse_link(object$family)

  # Monte Carlo integration
  result <- mc_integrate_marginal(
    X = as.matrix(X),
    Z = as.matrix(Z),
    beta = beta,
    Sigma_u = Sigma_u,
    inv_link_fn = inv_link,
    n_sim = n_sim
  )

  if (se.fit) {
    return(list(fit = result$fit, se.fit = result$se))
  } else {
    return(result$fit)
  }
}
