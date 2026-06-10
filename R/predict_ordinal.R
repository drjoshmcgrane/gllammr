#' Predict method for ordinal models
#'
#' Obtain predictions from a fitted ordinal regression model
#'
#' @param object A fitted ordinal model (class gllamm with ordinal_family)
#' @param newdata Optional new data frame for predictions
#' @param type Type of prediction:
#'   \describe{
#'     \item{class}{Predicted class (modal category)}
#'     \item{probs}{Conditional probabilities for each category}
#'     \item{cumprobs}{Conditional cumulative probabilities}
#'     \item{marginal}{Marginal probabilities (population-averaged)}
#'   }
#' @param n_sim Number of Monte Carlo samples for marginal predictions (default: 1000)
#' @param ... Additional arguments (currently unused)
#'
#' @return Depends on \code{type}:
#'   \itemize{
#'     \item{class: Vector of predicted classes}
#'     \item{probs: Matrix of probabilities (n_obs × n_categories)}
#'     \item{cumprobs: Matrix of cumulative probabilities}
#'     \item{marginal: Matrix of marginal probabilities (n_obs × n_categories)}
#'   }
#'
#' @examples
#' \dontrun{
#' # Fit ordinal model
#' fit <- gllamm(rating ~ temp + (1 | judge),
#'               data = wine,
#'               family = ordinal(link = "logit"))
#'
#' # Predicted classes
#' pred_class <- predict(fit, type = "class")
#'
#' # Conditional probabilities
#' pred_probs <- predict(fit, type = "probs")
#'
#' # Marginal probabilities (population-averaged)
#' pred_marg <- predict(fit, type = "marginal")
#' }
#'
#' @export
predict.gllamm_ordinal <- function(object,
                                   newdata = NULL,
                                   type = c("class", "probs", "cumprobs", "marginal"),
                                   n_sim = 1000,
                                   ...) {
  type <- match.arg(type)

  # Extract model components
  n_categories <- object$n_categories
  beta <- object$coefficients$fixed
  thresholds <- object$coefficients$thresholds

  # Get model matrices
  if (is.null(newdata)) {
    X <- object$X
    # Need to reconstruct Z for original data
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

  # Get link function info
  link <- object$family$link

  # Marginal predictions require Monte Carlo
  if (type == "marginal") {
    return(predict_marginal_ordinal(object, X, Z, n_sim))
  }

  # Conditional predictions (at u=0)
  # Compute linear predictor (fixed effects only)
  eta <- as.vector(X %*% beta)

  # Compute cumulative probabilities for each threshold
  cumprobs <- matrix(NA, n_obs, n_categories - 1)

  for (k in 1:(n_categories - 1)) {
    if (link == "logit") {
      cumprobs[, k] <- plogis(thresholds[k] - eta)
    } else if (link == "probit") {
      cumprobs[, k] <- pnorm(thresholds[k] - eta)
    } else {
      stop("Unknown ordinal link:", link)
    }
  }

  if (type == "cumprobs") {
    return(cumprobs)
  }

  # Convert cumulative to category probabilities
  probs <- matrix(NA, n_obs, n_categories)
  probs[, 1] <- cumprobs[, 1]
  for (k in 2:(n_categories - 1)) {
    probs[, k] <- cumprobs[, k] - cumprobs[, k - 1]
  }
  probs[, n_categories] <- 1 - cumprobs[, n_categories - 1]

  # Ensure probabilities sum to 1 (numerical precision)
  probs <- probs / rowSums(probs)

  if (type == "probs") {
    colnames(probs) <- paste0("P(Y=", 1:n_categories, ")")
    return(probs)
  }

  if (type == "class") {
    # Return modal category
    predicted_class <- apply(probs, 1, which.max)
    return(predicted_class)
  }
}


#' Internal function for marginal ordinal predictions
#'
#' @param object Fitted ordinal model
#' @param X Fixed effects design matrix
#' @param Z Random effects design matrix
#' @param n_sim Number of MC samples
#'
#' @return Matrix of marginal probabilities (n_obs × n_categories)
#' @keywords internal
predict_marginal_ordinal <- function(object, X, Z, n_sim = 1000) {
  n_obs <- nrow(X)
  n_categories <- object$n_categories
  beta <- object$coefficients$fixed
  thresholds <- object$coefficients$thresholds
  link <- object$family$link

  # Extract random effects variance
  Sigma_u <- extract_random_vcov(object)
  n_random <- ncol(Sigma_u)

  # Draw samples from random effects distribution
  if (n_random == 1) {
    sigma_u <- sqrt(Sigma_u[1, 1])
    u_samples <- matrix(rnorm(n_sim, 0, sigma_u), n_sim, 1)
  } else {
    u_samples <- rmvnorm_chol(n_sim, Sigma_u)
  }

  # Storage for marginal probabilities
  marginal_probs <- matrix(0, n_obs, n_categories)

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

    # Total linear predictor
    eta <- eta_fixed + eta_random

    # Compute cumulative probabilities
    cumprobs <- matrix(NA, n_obs, n_categories - 1)

    for (k in 1:(n_categories - 1)) {
      if (link == "logit") {
        cumprobs[, k] <- plogis(thresholds[k] - eta)
      } else if (link == "probit") {
        cumprobs[, k] <- pnorm(thresholds[k] - eta)
      } else {
        stop("Unknown ordinal link:", link)
      }
    }

    # Convert to category probabilities
    probs <- matrix(NA, n_obs, n_categories)
    probs[, 1] <- cumprobs[, 1]
    for (k in 2:(n_categories - 1)) {
      probs[, k] <- cumprobs[, k] - cumprobs[, k - 1]
    }
    probs[, n_categories] <- 1 - cumprobs[, n_categories - 1]

    # Accumulate
    marginal_probs <- marginal_probs + probs
  }

  # Average across samples
  marginal_probs <- marginal_probs / n_sim

  # Ensure probabilities sum to 1 (numerical precision)
  marginal_probs <- marginal_probs / rowSums(marginal_probs)

  colnames(marginal_probs) <- paste0("P(Y=", 1:n_categories, ")")

  return(marginal_probs)
}
