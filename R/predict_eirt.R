#' Predict method for EIRT models
#'
#' Obtain predictions from a fitted Explanatory IRT model
#'
#' @param object A fitted EIRT model (class gllamm_eirt)
#' @param newdata Optional data frame with item covariates for new items.
#'   Must include all variables used in difficulty_formula and discrimination_formula.
#' @param type Type of prediction:
#'   \describe{
#'     \item{probability}{Item response probabilities P(Y=1|θ, item_covariates)}
#'     \item{ability}{Person ability estimates (θ)}
#'     \item{difficulty}{Predicted item difficulties}
#'     \item{discrimination}{Predicted item discriminations}
#'     \item{marginal}{Marginal item response probabilities E[P(Y=1|θ)]}
#'   }
#' @param ability Optional vector of ability values. If NULL, uses estimated abilities.
#' @param n_sim Number of Monte Carlo samples for marginal predictions (default: 1000)
#' @param ... Additional arguments (currently unused)
#'
#' @return Depends on \code{type}
#'
#' @examples
#' \dontrun{
#' # Fit EIRT model
#' responses <- matrix(rbinom(500, 1, 0.6), 50, 10)
#' item_data <- data.frame(
#'   item_id = 1:10,
#'   difficulty_pred = rnorm(10)
#' )
#'
#' fit <- fit_eirt(responses, item_data,
#'                 difficulty_formula = ~ difficulty_pred,
#'                 model = "2PL")
#'
#' # Predicted item difficulties
#' pred_diff <- predict(fit, type = "difficulty")
#'
#' # Marginal item response probabilities
#' marg_probs <- predict(fit, type = "marginal")
#'
#' # For new items
#' new_items <- data.frame(difficulty_pred = c(-1, 0, 1))
#' pred_new <- predict(fit, newdata = new_items, type = "marginal")
#' }
#'
#' @export
predict.gllamm_eirt <- function(object,
                                newdata = NULL,
                                type = c("probability", "ability", "difficulty",
                                        "discrimination", "marginal"),
                                ability = NULL,
                                n_sim = 1000,
                                ...) {
  type <- match.arg(type)

  # Type: ability
  if (type == "ability") {
    return(object$person_abilities)
  }

  # Type: difficulty (for original or new items)
  if (type == "difficulty") {
    if (is.null(newdata)) {
      return(object$item_parameters$difficulty)
    } else {
      # Predict difficulty for new items
      gamma <- object$regression_coefficients$difficulty
      W <- model.matrix(object$difficulty_formula, data = newdata)
      pred_diff <- as.vector(W %*% gamma)
      return(pred_diff)
    }
  }

  # Type: discrimination
  if (type == "discrimination") {
    if (is.null(newdata)) {
      return(object$item_parameters$discrimination)
    } else {
      # Predict discrimination for new items
      delta <- object$regression_coefficients$discrimination
      W <- model.matrix(object$discrimination_formula, data = newdata)
      log_disc_pred <- as.vector(W %*% delta)
      pred_disc <- exp(log_disc_pred)
      return(pred_disc)
    }
  }

  # For probability and marginal predictions, need item parameters
  if (is.null(newdata)) {
    difficulty <- object$item_parameters$difficulty
    discrimination <- object$item_parameters$discrimination
    n_items <- length(difficulty)
  } else {
    # Predict parameters for new items
    gamma <- object$regression_coefficients$difficulty
    delta <- object$regression_coefficients$discrimination

    W_diff <- model.matrix(object$difficulty_formula, data = newdata)
    W_disc <- model.matrix(object$discrimination_formula, data = newdata)

    difficulty <- as.vector(W_diff %*% gamma)
    discrimination <- exp(as.vector(W_disc %*% delta))
    n_items <- length(difficulty)
  }

  # Type: marginal
  if (type == "marginal") {
    sigma_theta <- object$ability_sd
    theta_samples <- rnorm(n_sim, mean = 0, sd = sigma_theta)

    marginal_probs <- numeric(n_items)

    for (j in 1:n_items) {
      b <- difficulty[j]
      a <- discrimination[j]

      if (object$model %in% c("Rasch", "2PL")) {
        if (a == 1) {  # Rasch
          probs_samples <- plogis(theta_samples - b)
        } else {  # 2PL
          probs_samples <- plogis(a * (theta_samples - b))
        }
      } else {
        stop("Marginal predictions for polytomous EIRT not yet implemented")
      }

      marginal_probs[j] <- mean(probs_samples)
    }

    names(marginal_probs) <- paste0("Item", 1:n_items)
    return(marginal_probs)
  }

  # Type: probability (conditional)
  if (is.null(ability)) {
    theta <- object$person_abilities
  } else {
    theta <- ability
  }

  n_persons <- length(theta)
  probs <- matrix(NA, n_persons, n_items)

  for (j in 1:n_items) {
    b <- difficulty[j]
    a <- discrimination[j]

    if (object$model == "Rasch" || (object$model == "2PL" && a == 1)) {
      probs[, j] <- plogis(theta - b)
    } else if (object$model == "2PL") {
      probs[, j] <- plogis(a * (theta - b))
    } else {
      stop("Polytomous EIRT probability predictions not yet fully implemented")
    }
  }

  colnames(probs) <- paste0("Item", 1:n_items)
  return(probs)
}
