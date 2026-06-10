#' Predict method for IRT models
#'
#' Obtain predictions from a fitted IRT model
#'
#' @param object A fitted IRT model (class gllamm_irt)
#' @param newdata Optional specification of items/persons for predictions.
#'   Can be:
#'   \itemize{
#'     \item{NULL: Predict for all original items}
#'     \item{Numeric vector: Item indices to predict}
#'     \item{Data frame with 'item' column: Specific items}
#'   }
#' @param type Type of prediction:
#'   \describe{
#'     \item{probability}{Item response probabilities P(Y=1|θ)}
#'     \item{ability}{Person ability estimates (θ)}
#'     \item{marginal}{Marginal item response probabilities E[P(Y=1|θ)]}
#'   }
#' @param ability Optional vector of ability values for which to compute probabilities.
#'   If NULL and type="probability", uses estimated abilities from fitted model.
#' @param n_sim Number of Monte Carlo samples for marginal predictions (default: 1000)
#' @param ... Additional arguments (currently unused)
#'
#' @return Depends on \code{type}:
#'   \itemize{
#'     \item{probability: Matrix of probabilities (n_persons × n_items) or vector if ability specified}
#'     \item{ability: Vector of person abilities}
#'     \item{marginal: Vector of marginal probabilities (one per item)}
#'   }
#'
#' @examples
#' \dontrun{
#' # Fit 2PL model
#' responses <- matrix(rbinom(1000, 1, 0.6), 100, 10)
#' fit <- fit_irt(responses, model = "2PL")
#'
#' # Person abilities
#' abilities <- predict(fit, type = "ability")
#'
#' # Item response probabilities for each person
#' probs <- predict(fit, type = "probability")
#'
#' # Marginal item response probabilities (population-level)
#' marg_probs <- predict(fit, type = "marginal")
#' }
#'
#' @export
predict.gllamm_irt <- function(object,
                               newdata = NULL,
                               type = c("probability", "ability", "marginal"),
                               ability = NULL,
                               n_sim = 1000,
                               ...) {
  type <- match.arg(type)

  # Type: ability
  if (type == "ability") {
    return(object$person_abilities)
  }

  # Extract item parameters
  n_items <- object$n_items
  difficulty <- object$item_parameters$difficulty
  discrimination <- object$item_parameters$discrimination

  # Determine which items to predict
  if (is.null(newdata)) {
    items <- 1:n_items
  } else if (is.numeric(newdata)) {
    items <- newdata
  } else if (is.data.frame(newdata) && "item" %in% names(newdata)) {
    items <- newdata$item
  } else {
    stop("newdata must be NULL, numeric vector of item indices, or data frame with 'item' column")
  }

  # Type: marginal
  if (type == "marginal") {
    return(predict_marginal_irt(object, items, n_sim))
  }

  # Type: probability (conditional)
  if (is.null(ability)) {
    # Use estimated person abilities
    theta <- object$person_abilities
    n_persons <- length(theta)
  } else {
    # Use provided ability values
    theta <- ability
    n_persons <- length(theta)
  }

  # Compute probabilities
  probs <- matrix(NA, n_persons, length(items))

  for (j in seq_along(items)) {
    item_idx <- items[j]
    b <- difficulty[item_idx]
    a <- discrimination[item_idx]

    if (object$model == "Rasch") {
      # P(Y=1) = logit^{-1}(θ - b)
      probs[, j] <- plogis(theta - b)
    } else if (object$model == "2PL") {
      # P(Y=1) = logit^{-1}(a*(θ - b))
      probs[, j] <- plogis(a * (theta - b))
    } else if (object$model == "3PL") {
      # P(Y=1) = c + (1-c)*logit^{-1}(a*(θ - b))
      c_param <- object$item_parameters$guessing[item_idx]
      probs[, j] <- c_param + (1 - c_param) * plogis(a * (theta - b))
    } else if (object$model %in% c("GRM", "PCM", "GPCM", "NRM")) {
      # Polytomous models - return probability of highest category
      # (For now - full polytomous prediction requires more work)
      stop("Polytomous IRT prediction not yet fully implemented. ",
           "Use predict on the response matrix directly.")
    } else {
      stop("Unknown IRT model:", object$model)
    }
  }

  if (length(items) == 1) {
    return(as.vector(probs))
  } else {
    colnames(probs) <- paste0("Item", items)
    return(probs)
  }
}


#' Internal function for marginal IRT predictions
#'
#' Computes E[P(Y=1|θ)] where θ ~ N(0, σ²_θ)
#'
#' @param object Fitted IRT model
#' @param items Item indices to predict
#' @param n_sim Number of MC samples
#'
#' @return Vector of marginal probabilities (one per item)
#' @keywords internal
predict_marginal_irt <- function(object, items, n_sim = 1000) {
  # Extract parameters
  difficulty <- object$item_parameters$difficulty[items]
  discrimination <- object$item_parameters$discrimination[items]
  sigma_theta <- object$ability_sd

  # Draw ability samples from population distribution
  theta_samples <- rnorm(n_sim, mean = 0, sd = sigma_theta)

  # Storage for marginal probabilities
  marginal_probs <- numeric(length(items))

  # For each item
  for (j in seq_along(items)) {
    b <- difficulty[j]
    a <- discrimination[j]

    # Compute P(Y=1|θ) for each sampled θ
    if (object$model == "Rasch") {
      probs_samples <- plogis(theta_samples - b)
    } else if (object$model == "2PL") {
      probs_samples <- plogis(a * (theta_samples - b))
    } else if (object$model == "3PL") {
      item_idx <- items[j]
      c_param <- object$item_parameters$guessing[item_idx]
      probs_samples <- c_param + (1 - c_param) * plogis(a * (theta_samples - b))
    } else {
      stop("Marginal predictions for model", object$model, "not yet implemented")
    }

    # Average across samples
    marginal_probs[j] <- mean(probs_samples)
  }

  names(marginal_probs) <- paste0("Item", items)
  return(marginal_probs)
}
