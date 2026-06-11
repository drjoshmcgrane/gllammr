#' Reconstruct per-item threshold parameters from a fitted polytomous EIRT model
#'
#' Mirrors the threshold construction in the TMB template (gllamm_eirt.hpp):
#' GRM uses ordered sum-to-zero deviations around the item location,
#' PCM/GPCM use sum-to-zero step deviations around the item location,
#' LPCM uses the threshold regression plus step residuals.
#'
#' @param object Fitted gllamm_eirt model (polytomous)
#' @return List of per-item threshold vectors on the absolute scale
#' @keywords internal
eirt_item_thresholds <- function(object) {
  par_full <- object$tmb_obj$env$last.par.best
  n_items <- object$n_items
  K_vec <- object$n_categories_per_item
  n_step_cols <- max(1L, object$max_categories - 1L)
  b <- object$item_parameters$difficulty
  pmt <- object$poly_model_type

  thresholds <- vector("list", n_items)

  if (pmt %in% c(1L, 2L, 3L)) {
    # parList restores the full matrix shape including map-fixed dead cells
    sp <- object$tmb_obj$env$parList(par = par_full)$step_param

    for (j in seq_len(n_items)) {
      K <- K_vec[j]
      if (pmt == 1L) {
        # GRM: ordered sum-to-zero deviations around the item location:
        # u_1 = 0, u_k = u_{k-1} + exp(sp[j, k-1]); tau = b + u - mean(u)
        u <- numeric(K - 1)
        if (K > 2) {
          for (k in 2:(K - 1)) {
            u[k] <- u[k - 1] + exp(sp[j, k - 1])
          }
        }
        tau <- b[j] + u - mean(u)
      } else {
        # PCM/GPCM: delta_m = b + s_m with sum-to-zero step deviations
        if (K > 2) {
          s <- c(sp[j, seq_len(K - 2)], -sum(sp[j, seq_len(K - 2)]))
        } else {
          s <- 0
        }
        tau <- b[j] + s
      }
      thresholds[[j]] <- tau
    }
  } else if (pmt == 4L) {
    # LPCM: delta_m = b + sum_p W[j,p] * xi[p,m] + e_step[j,m]
    xi <- object$regression_coefficients$threshold
    W_thresh <- model.matrix(object$formulas$threshold, data = object$item_data)
    e_step <- matrix(par_full[names(par_full) == "e_step"],
                     n_items, n_step_cols)

    for (j in seq_len(n_items)) {
      K <- K_vec[j]
      tau <- numeric(K - 1)
      for (m in seq_len(K - 1)) {
        tau[m] <- b[j] + sum(W_thresh[j, ] * xi[, m]) + e_step[j, m]
      }
      thresholds[[j]] <- tau
    }
  } else {
    stop("Threshold reconstruction requires a polytomous EIRT model")
  }

  names(thresholds) <- paste0("Item", seq_len(n_items))
  thresholds
}


#' Map an EIRT poly_model_type code to the shared probability helper's model name
#' @keywords internal
eirt_poly_model_name <- function(poly_model_type) {
  switch(as.character(poly_model_type),
         "1" = "GRM",
         "2" = "PCM",
         "3" = "GPCM",
         "4" = "PCM",   # LPCM likelihood is adjacent-categories without discrimination
         stop("Unknown poly_model_type: ", poly_model_type))
}


#' Predict method for EIRT models
#'
#' Obtain predictions from a fitted Explanatory IRT model.
#'
#' @param object A fitted EIRT model (class gllamm_eirt)
#' @param newdata Optional data frame with item covariates for new items.
#'   Must include all variables used in the difficulty and discrimination
#'   formulas. For polytomous models, new-item predictions use the predicted
#'   item location with step deviations at their population value of zero
#'   (PCM/GPCM) or the threshold regression (LPCM); GRM threshold spacing is
#'   item-specific and cannot be predicted for new items.
#' @param type Type of prediction:
#'   \describe{
#'     \item{probability}{Item response probabilities. Dichotomous: persons x
#'       items matrix of P(Y=1). Polytomous: list of persons x categories
#'       matrices, one per item.}
#'     \item{ability}{Person ability estimates (theta)}
#'     \item{difficulty}{Predicted item difficulties}
#'     \item{discrimination}{Predicted item discriminations}
#'     \item{marginal}{Marginal response probabilities, integrating ability
#'       over its population distribution. Dichotomous: vector of E[P(Y=1)].
#'       Polytomous: items x categories matrix of E[P(Y=k)].}
#'   }
#' @param ability Optional vector of ability values. If NULL, uses estimated abilities.
#' @param n_sim Number of Monte Carlo samples for marginal predictions (default: 1000)
#' @param ... Additional arguments (currently unused)
#'
#' @return Depends on \code{type}
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
    }
    gamma <- object$regression_coefficients$difficulty
    W <- model.matrix(object$formulas$difficulty, data = newdata)
    return(as.vector(W %*% gamma))
  }

  # Type: discrimination
  if (type == "discrimination") {
    if (is.null(newdata)) {
      return(object$item_parameters$discrimination)
    }
    delta <- object$regression_coefficients$discrimination
    W <- model.matrix(object$formulas$discrimination, data = newdata)
    return(exp(as.vector(W %*% delta)))
  }

  # For probability and marginal predictions, assemble item parameters
  is_poly <- isTRUE(object$is_polytomous)

  if (is.null(newdata)) {
    difficulty <- object$item_parameters$difficulty
    discrimination <- object$item_parameters$discrimination
    n_items <- length(difficulty)
    if (is_poly) {
      thresholds <- eirt_item_thresholds(object)
    }
  } else {
    gamma <- object$regression_coefficients$difficulty
    W_diff <- model.matrix(object$formulas$difficulty, data = newdata)
    difficulty <- as.vector(W_diff %*% gamma)
    n_items <- length(difficulty)

    delta <- object$regression_coefficients$discrimination
    if (all(is.na(delta))) {
      discrimination <- rep(1, n_items)
    } else {
      W_disc <- model.matrix(object$formulas$discrimination, data = newdata)
      discrimination <- exp(as.vector(W_disc %*% delta))
    }

    if (is_poly) {
      pmt <- object$poly_model_type
      K <- object$max_categories
      if (pmt == 1L) {
        stop("GRM threshold spacing is item-specific and cannot be predicted ",
             "for new items. Use PCM/GPCM, or a threshold_formula (LPCM).")
      }
      thresholds <- vector("list", n_items)
      if (pmt == 4L) {
        xi <- object$regression_coefficients$threshold
        W_thresh <- model.matrix(object$formulas$threshold, data = newdata)
        for (j in seq_len(n_items)) {
          thresholds[[j]] <- difficulty[j] +
            as.vector(t(xi) %*% W_thresh[j, ])  # e_step at population value 0
        }
      } else {
        # PCM/GPCM: step deviations at their population value of zero
        for (j in seq_len(n_items)) {
          thresholds[[j]] <- rep(difficulty[j], K - 1)
        }
      }
    }
  }

  # Type: marginal — integrate ability over N(0, sigma_theta^2)
  if (type == "marginal") {
    sigma_theta <- object$ability_sd
    theta_samples <- rnorm(n_sim, mean = 0, sd = sigma_theta)

    if (is_poly) {
      model_name <- eirt_poly_model_name(object$poly_model_type)
      max_K <- max(vapply(thresholds, length, integer(1))) + 1L
      out <- matrix(NA_real_, n_items, max_K)
      for (j in seq_len(n_items)) {
        p <- irt_category_probs(model_name, theta_samples,
                                thresholds[[j]], discrimination[j])
        out[j, seq_len(ncol(p))] <- colMeans(p)
      }
      rownames(out) <- paste0("Item", seq_len(n_items))
      colnames(out) <- paste0("Category", seq_len(max_K))
      return(out)
    }

    marginal_probs <- numeric(n_items)
    for (j in seq_len(n_items)) {
      eta <- discrimination[j] * (theta_samples - difficulty[j])
      marginal_probs[j] <- mean(plogis(eta))
    }
    names(marginal_probs) <- paste0("Item", seq_len(n_items))
    return(marginal_probs)
  }

  # Type: probability (conditional on ability)
  theta <- if (is.null(ability)) object$person_abilities else ability
  n_persons <- length(theta)

  if (is_poly) {
    model_name <- eirt_poly_model_name(object$poly_model_type)
    probs <- vector("list", n_items)
    for (j in seq_len(n_items)) {
      p <- irt_category_probs(model_name, theta,
                              thresholds[[j]], discrimination[j])
      colnames(p) <- paste0("Category", seq_len(ncol(p)))
      probs[[j]] <- p
    }
    names(probs) <- paste0("Item", seq_len(n_items))
    return(probs)
  }

  probs <- matrix(NA_real_, n_persons, n_items)
  for (j in seq_len(n_items)) {
    probs[, j] <- plogis(discrimination[j] * (theta - difficulty[j]))
  }
  colnames(probs) <- paste0("Item", seq_len(n_items))
  probs
}
