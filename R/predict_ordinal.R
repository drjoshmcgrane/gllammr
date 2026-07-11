#' Category probabilities for every ordinal link
#'
#' Mirrors the likelihood construction in gllamm_ordinal.hpp. For the
#' PPO link \code{eta} must be an n x (K-1) matrix of per-threshold
#' linear predictors; for all other links a length-n vector.
#'
#' @keywords internal
.ordinal_category_probs <- function(eta, thresholds, link, K) {
  n <- if (is.matrix(eta)) nrow(eta) else length(eta)
  probs <- matrix(NA_real_, n, K)

  if (link %in% c("logit", "probit", "ppo")) {
    pfun <- if (link == "probit") stats::pnorm else stats::plogis
    cum <- vapply(seq_len(K - 1), function(k) {
      e <- if (is.matrix(eta)) eta[, k] else eta
      pfun(thresholds[k] - e)
    }, numeric(n))
    cum <- matrix(cum, n, K - 1)
    probs[, 1] <- cum[, 1]
    if (K > 2) {
      for (k in 2:(K - 1)) probs[, k] <- cum[, k] - cum[, k - 1]
    }
    probs[, K] <- 1 - cum[, K - 1]
    # PPO curves can cross at extreme covariate values
    probs <- pmax(probs, 1e-12)
    probs <- probs / rowSums(probs)
  } else if (link == "acl") {
    lp <- matrix(0, n, K)
    for (k in 2:K) lp[, k] <- lp[, k - 1] + thresholds[k - 1] + eta
    m <- lp[, 1]
    for (k in 2:K) m <- pmax(m, lp[, k])
    el <- exp(lp - m)
    probs <- el / rowSums(el)
  } else if (link == "crl_forward") {
    surv <- rep(1, n)
    for (k in seq_len(K - 1)) {
      h <- stats::plogis(thresholds[k] - eta)
      probs[, k] <- surv * h
      surv <- surv * (1 - h)
    }
    probs[, K] <- surv
  } else if (link == "crl_backward") {
    # P(c) = b_c prod_{j>c}(1 - b_j), b_c = plogis(tau_{c-1} - eta)
    surv <- rep(1, n)
    for (c in K:2) {
      b <- stats::plogis(thresholds[c - 1] - eta)
      probs[, c] <- b * surv
      surv <- surv * (1 - b)
    }
    probs[, 1] <- surv
  } else {
    stop("Unknown ordinal link: ", link)
  }
  probs
}


#' Per-term random-effects pieces of a fitted ordinal model
#'
#' Returns aligned lists of Z matrices, 0-based group indices, group
#' counts, and covariance matrices - one element per random-effects term,
#' for both single- and multi-term fits.
#'
#' @keywords internal
.ordinal_re_parts <- function(object, data) {
  parsed <- parse_formula(object$formula, data)
  model_data <- make_model_matrices(parsed, data)
  rv <- object$coefficients$random_var
  Sigmas <- if (is.list(rv)) {
    lapply(rv, as.matrix)
  } else {
    list(as.matrix(extract_random_vcov(object)))
  }
  list(Z = model_data$Z, groups = model_data$groups,
       n_groups = model_data$n_groups, Sigmas = Sigmas,
       X = drop_intercept_column(model_data$X), n_obs = model_data$n_obs)
}


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
  warn_not_converged(object)
  type <- match.arg(type)

  # Extract model components
  n_categories <- object$n_categories
  beta <- object$coefficients$fixed
  thresholds <- object$coefficients$thresholds

  # Get model matrices
  if (is.null(newdata)) {
    X <- drop_intercept_column(object$X)  # thresholds carry the location
    # Need to reconstruct Z for original data
    parsed <- parse_formula(object$formula, object$data)
    model_data <- make_model_matrices(parsed, object$data)
    Z <- model_data$Z[[1]]
  } else {
    parsed <- parse_formula(object$formula, newdata)
    new_mats <- make_model_matrices(parsed, newdata)
    X <- drop_intercept_column(new_mats$X)  # thresholds carry the location
    Z <- new_mats$Z[[1]]
  }

  n_obs <- nrow(X)
  link <- object$link

  # Linear predictor at u = 0: a vector for links 1-5, an n x (K-1)
  # matrix of per-threshold predictors for PPO
  eta <- if (link == "ppo") {
    X %*% t(object$coefficients$beta_ppo)
  } else {
    as.vector(X %*% beta)
  }

  # Marginal predictions require Monte Carlo over all RE terms
  if (type == "marginal") {
    return(predict_marginal_ordinal(object, n_sim = n_sim,
                                    newdata = newdata))
  }

  if (type == "cumprobs") {
    if (!link %in% c("logit", "probit", "ppo")) {
      stop("Cumulative probabilities are defined for the cumulative ",
           "links (logit, probit, ppo); use type = \"probs\" for ", link)
    }
    pfun <- if (link == "probit") pnorm else plogis
    cumprobs <- vapply(seq_len(n_categories - 1), function(k) {
      e <- if (is.matrix(eta)) eta[, k] else eta
      pfun(thresholds[k] - e)
    }, numeric(n_obs))
    return(matrix(cumprobs, n_obs, n_categories - 1))
  }

  probs <- .ordinal_category_probs(eta, thresholds, link, n_categories)

  if (type == "probs") {
    colnames(probs) <- paste0("P(Y=", 1:n_categories, ")")
    return(probs)
  }
  apply(probs, 1, which.max)
}


#' Internal function for marginal ordinal predictions
#'
#' Population-averaged category probabilities by Monte Carlo over the
#' estimated random-effects distributions of every term (fresh draws per
#' replicate; the same draw applies to all members of a group only in
#' expectation, which is what the marginal quantity requires).
#'
#' @keywords internal
predict_marginal_ordinal <- function(object, n_sim = 1000,
                                     newdata = NULL) {
  data <- if (is.null(newdata)) object$data else newdata
  parts <- .ordinal_re_parts(object, data)
  X <- parts$X
  n_obs <- parts$n_obs
  n_categories <- object$n_categories
  thresholds <- object$coefficients$thresholds
  link <- object$link

  eta_fixed <- if (link == "ppo") {
    X %*% t(object$coefficients$beta_ppo)
  } else {
    as.vector(X %*% object$coefficients$fixed)
  }

  marginal_probs <- matrix(0, n_obs, n_categories)
  for (s in seq_len(n_sim)) {
    # One population draw per RE term, mapped to observations
    eta_random <- rep(0, n_obs)
    for (t in seq_along(parts$Sigmas)) {
      u_t <- rmvnorm_chol(1, parts$Sigmas[[t]])
      eta_random <- eta_random +
        as.vector(parts$Z[[t]] %*% as.numeric(u_t))
    }
    eta <- if (is.matrix(eta_fixed)) eta_fixed + eta_random
           else eta_fixed + eta_random
    marginal_probs <- marginal_probs +
      .ordinal_category_probs(eta, thresholds, link, n_categories)
  }
  marginal_probs <- marginal_probs / n_sim
  marginal_probs <- marginal_probs / rowSums(marginal_probs)
  colnames(marginal_probs) <- paste0("P(Y=", 1:n_categories, ")")
  marginal_probs
}


#' Simulate from a fitted ordinal model
#'
#' Draws fresh random effects from their estimated distribution for each
#' replicate (population-level simulation) and samples categories from the
#' implied response distribution.
#'
#' @param object A fitted \code{gllamm_ordinal} object
#' @param nsim Number of simulations (default: 1)
#' @param seed Optional random seed (stored in the \code{"seed"} attribute)
#' @param newdata Optional new data frame with the model covariates and
#'   grouping variables
#' @param ... Additional arguments (currently unused)
#'
#' @return A data frame with \code{nsim} columns of simulated category
#'   responses (integer codes 1..K), with a \code{"seed"} attribute.
#'
#' @export
simulate.gllamm_ordinal <- function(object,
                                    nsim = 1,
                                    seed = NULL,
                                    newdata = NULL,
                                    ...) {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    runif(1)
  }
  if (is.null(seed)) {
    rng_state <- get(".Random.seed", envir = globalenv())
  } else {
    set.seed(seed)
    rng_state <- structure(seed, kind = as.list(RNGkind()))
  }

  link <- object$link
  n_categories <- object$n_categories
  thresholds <- object$coefficients$thresholds

  sim_data <- if (is.null(newdata)) object$data else newdata
  parts <- .ordinal_re_parts(object, sim_data)
  X <- parts$X
  n_obs <- parts$n_obs

  eta_fixed <- if (link == "ppo") {
    X %*% t(object$coefficients$beta_ppo)
  } else {
    as.vector(X %*% object$coefficients$fixed)
  }

  sims <- matrix(NA_integer_, n_obs, nsim)
  for (s in seq_len(nsim)) {
    # Fresh group-level draws for every RE term (population semantics)
    eta_random <- rep(0, n_obs)
    for (t in seq_along(parts$Sigmas)) {
      u_t <- rmvnorm_chol(parts$n_groups[t], parts$Sigmas[[t]])
      eta_random <- eta_random +
        rowSums(parts$Z[[t]] *
                  u_t[parts$groups[[t]] + 1, , drop = FALSE])
    }
    eta <- eta_fixed + eta_random
    probs <- .ordinal_category_probs(eta, thresholds, link, n_categories)
    cum <- t(apply(probs, 1, cumsum))
    r <- runif(n_obs)
    sims[, s] <- 1L + rowSums(r > cum[, -n_categories, drop = FALSE])
  }

  out <- as.data.frame(sims)
  names(out) <- paste0("sim_", seq_len(nsim))
  attr(out, "seed") <- rng_state
  out
}
