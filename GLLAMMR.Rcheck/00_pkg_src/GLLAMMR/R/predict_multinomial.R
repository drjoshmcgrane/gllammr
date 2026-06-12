#' Predict method for multinomial models
#'
#' Obtain predictions from a fitted multinomial logit model with random effects
#'
#' @param object A fitted multinomial model
#' @param newdata Optional new data frame for predictions
#' @param type Type of prediction:
#'   \describe{
#'     \item{class}{Predicted class (modal category)}
#'     \item{probs}{Conditional probabilities for each category}
#'     \item{marginal}{Marginal probabilities (population-averaged)}
#'   }
#' @param n_sim Number of Monte Carlo samples for marginal predictions (default: 1000)
#' @param ... Additional arguments (currently unused)
#'
#' @return Depends on \code{type}:
#'   \itemize{
#'     \item{class: Vector of predicted classes}
#'     \item{probs: Matrix of probabilities (n_obs × n_categories)}
#'     \item{marginal: Matrix of marginal probabilities (n_obs × n_categories)}
#'   }
#'
#' @rdname predict_multinomial
#' @export
predict.gllamm_multinomial <- function(object,
                                       newdata = NULL,
                                       type = c("class", "probs", "marginal"),
                                       n_sim = 1000,
                                       ...) {
  predict_multinomial(object, newdata = newdata, type = type,
                      n_sim = n_sim, ...)
}


#' @keywords internal
predict_multinomial <- function(object,
                               newdata = NULL,
                               type = c("class", "probs", "marginal"),
                               n_sim = 1000,
                               ...) {
  type <- match.arg(type)

  # Extract model components
  n_categories <- object$n_categories
  beta <- object$coefficients$beta  # (n_categories-1) × n_fixed matrix
  n_fixed <- ncol(beta)

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

  # Marginal predictions
  if (type == "marginal") {
    # Extract random effects variance
    Sigma_u <- extract_random_vcov(object)
    n_random <- ncol(Sigma_u)

    # Draw samples
    if (n_random == 1) {
      sigma_u <- sqrt(Sigma_u[1, 1])
      u_samples <- matrix(rnorm(n_sim, 0, sigma_u), n_sim, 1)
    } else {
      u_samples <- rmvnorm_chol(n_sim, Sigma_u)
    }

    # Storage
    marginal_probs <- matrix(0, n_obs, n_categories)

    # For each MC sample
    for (s in 1:n_sim) {
      u <- u_samples[s, ]

      # Random effects contribution
      if (is.matrix(Z)) {
        eta_random <- as.vector(Z %*% u)
      } else {
        eta_random <- Z * u
      }

      # Compute probabilities for this sample
      probs_sample <- compute_multinomial_probs(X, beta, eta_random, n_categories)

      # Accumulate
      marginal_probs <- marginal_probs + probs_sample
    }

    # Average
    marginal_probs <- marginal_probs / n_sim
    marginal_probs <- marginal_probs / rowSums(marginal_probs)  # Normalize

    colnames(marginal_probs) <- paste0("P(Y=", 0:(n_categories-1), ")")
    return(marginal_probs)
  }

  # Conditional predictions (at u=0)
  probs <- compute_multinomial_probs(X, beta, rep(0, n_obs), n_categories)

  if (type == "probs") {
    colnames(probs) <- paste0("P(Y=", 0:(n_categories-1), ")")
    return(probs)
  }

  if (type == "class") {
    predicted_class <- apply(probs, 1, which.max) - 1  # 0-indexed
    return(predicted_class)
  }
}


#' Compute multinomial probabilities
#'
#' @param X Fixed effects design matrix
#' @param beta Beta matrix (K-1) × p
#' @param eta_random Random effects contribution (vector of length n_obs)
#' @param n_categories Number of categories
#'
#' @return Matrix of probabilities (n_obs × n_categories)
#' @keywords internal
compute_multinomial_probs <- function(X, beta, eta_random, n_categories) {
  n_obs <- nrow(X)
  n_fixed <- ncol(X)

  # Compute linear predictors for each category
  eta <- matrix(0, n_obs, n_categories)

  # Category 0 is reference (eta = 0)
  # Categories 1 to K-1 have their own parameters
  for (cat in 2:n_categories) {
    beta_cat <- beta[cat - 1, ]
    eta[, cat] <- as.vector(X %*% beta_cat) + eta_random
  }

  # Softmax
  exp_eta <- exp(eta)
  probs <- exp_eta / rowSums(exp_eta)

  return(probs)
}
