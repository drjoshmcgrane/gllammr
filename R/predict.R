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

  warn_not_converged(object)
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
        # Include random effects (some fitters store fitted_values)
        return(fitted(object))
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

    # Get group identifiers from new data. Use the freshly parsed
    # formula, not object$random_terms - some fitters (fit_binomial)
    # do not store the parsed terms on the result.
    if (length(parsed$random_terms) > 0) {
      parts <- .gllamm_re_parts(object, object$data)
      random_contrib <- numeric(nrow(new_mats$X))
      for (t in seq_len(parts$n_terms)) {
        rt <- parsed$random_terms[[t]]
        gv <- rt$grouping_vars %||% rt$grouping
        train_f <- if (length(gv) > 1) {
          interaction(object$data[, gv], drop = TRUE)
        } else {
          factor(object$data[[gv]])
        }
        new_f <- if (length(gv) > 1) {
          interaction(newdata[, gv], drop = TRUE)
        } else {
          factor(newdata[[gv]])
        }
        if (!is.null(new_mats$complete_idx)) {
          new_f <- new_f[new_mats$complete_idx]
        }
        idx <- match(as.character(new_f), levels(train_f))
        u_t <- parts$u[[t]]
        if (is.null(u_t)) next
        contrib <- rowSums(new_mats$Z[[t]] * u_t[idx, , drop = FALSE])
        contrib[is.na(idx)] <- 0      # new groups: population average
        random_contrib <- random_contrib + contrib
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


#' Per-term random-effects pieces of a fitted GLMM
#'
#' Normalizes the two stored shapes (single-term: per-group list of
#' coefficient vectors; multi-term: per-term list of group x coef
#' matrices) into aligned per-term lists of design matrices, group
#' indices, covariance matrices, and BLUP matrices.
#'
#' @keywords internal
.gllamm_re_parts <- function(object, data) {
  parsed <- parse_formula(object$formula, data)
  md <- make_model_matrices(parsed, data)
  rv <- object$coefficients$random_var
  to_mat <- function(m, nr) {
    if (is.matrix(m)) m else diag(as.numeric(m), nrow = nr)
  }
  Sigmas <- if (is.list(rv)) {
    lapply(seq_along(rv), function(t) {
      to_mat(rv[[t]], ncol(md$Z[[min(t, length(md$Z))]]))
    })
  } else {
    list(to_mat(rv, ncol(md$Z[[1]])))
  }
  n_terms <- length(Sigmas)

  re <- object$random_effects
  u_list <- if (is.list(re) && n_terms > 1 && length(re) == n_terms &&
                all(vapply(re, is.matrix, TRUE))) {
    re                                        # per-term matrices
  } else if (is.list(re) && length(re) > 0 &&
             !any(vapply(re, is.matrix, TRUE))) {
    list(do.call(rbind, re))                  # per-group vectors
  } else if (is.matrix(re)) {
    list(re)
  } else {
    vector("list", n_terms)
  }

  list(md = md, parsed = parsed, Sigmas = Sigmas, u = u_list,
       n_terms = n_terms)
}


#' Simulate from a GLLAMM model
#'
#' Simulate response data from a fitted GLLAMM model. Random effects are
#' drawn fresh from their estimated distribution for every replicate
#' (population-level simulation), both for the original data and for
#' \code{newdata}.
#'
#' @param object A fitted \code{gllamm} object
#' @param nsim Number of simulations (default: 1)
#' @param seed Optional random seed (stored in the \code{"seed"} attribute,
#'   following the \code{\link[stats]{simulate}} contract)
#' @param newdata Optional new data frame containing the covariates and
#'   grouping variables of the model formula
#' @param ... Additional arguments (currently unused)
#'
#' @return A data frame with \code{nsim} columns, one simulated response
#'   vector per column, with a \code{"seed"} attribute.
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

  # stats::simulate contract: record the RNG state / seed used
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    runif(1)
  }
  if (is.null(seed)) {
    rng_state <- get(".Random.seed", envir = globalenv())
  } else {
    set.seed(seed)
    rng_state <- structure(seed, kind = as.list(RNGkind()))
  }

  beta <- object$coefficients$fixed

  # Design matrices: original data or newdata via the same construction
  sim_data <- if (is.null(newdata)) object$data else newdata
  parts <- .gllamm_re_parts(object, sim_data)
  model_data <- parts$md
  X <- model_data$X
  n_obs <- model_data$n_obs

  family_name <- object$family$family

  # Residual SD for gaussian: the actual estimate, not a residual-variance
  # approximation
  sigma <- NULL
  if (family_name == "gaussian") {
    log_sigma <- object$tmb_opt$par["log_sigma"]
    sigma <- if (!is.na(log_sigma)) {
      exp(unname(log_sigma))
    } else {
      sqrt(var(object$residuals))
    }
  }

  fixed_part <- as.numeric(X %*% beta)

  sims <- matrix(NA_real_, nrow = n_obs, ncol = nsim)
  for (s in seq_len(nsim)) {
    # Fresh random effects each replicate, for EVERY term:
    # u_g ~ MVN(0, Sigma_t)
    eta <- fixed_part
    for (t in seq_len(parts$n_terms)) {
      u_sim <- rmvnorm_chol(model_data$n_groups[t], parts$Sigmas[[t]])
      eta <- eta + rowSums(model_data$Z[[t]] *
                             u_sim[model_data$groups[[t]] + 1, ,
                                   drop = FALSE])
    }

    sims[, s] <- switch(family_name,
      gaussian = rnorm(n_obs, mean = eta, sd = sigma),
      binomial = rbinom(n_obs, size = 1,
                        prob = object$family$linkinv(eta)),
      poisson = rpois(n_obs, lambda = exp(eta)),
      stop("Simulation not implemented for family: ", family_name)
    )
  }

  out <- as.data.frame(sims)
  names(out) <- paste0("sim_", seq_len(nsim))
  attr(out, "seed") <- rng_state
  out
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

  if (!is.numeric(n_sim) || length(n_sim) != 1 || !is.finite(n_sim) || n_sim < 1) {
    stop("'n_sim' must be a positive integer")
  }
  n_sim <- as.integer(n_sim)

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

  # Per-term design pieces (multi-term fits integrate over every term)
  mdata <- if (is.null(newdata)) object$data else newdata
  parts <- .gllamm_re_parts(object, mdata)
  X <- if (is.null(newdata) && !is.null(object$X)) object$X
       else parts$md$X

  beta <- object$coefficients$fixed
  inv_link <- get_inverse_link(object$family)

  # Monte Carlo integration, vectorized: draw all n_sim population-level
  # random effects up front (one draw per term per replicate) and reduce the
  # inverse-link probabilities column-wise. The draws are generated in the
  # exact order the former per-replicate loop consumed them - all terms of
  # replicate 1, then replicate 2, ... - so a fixed seed reproduces the old
  # results bit-for-bit while replacing n_sim R-level iterations with a couple
  # of matrix operations.
  eta_fixed <- as.numeric(as.matrix(X) %*% beta)

  q_terms <- vapply(parts$Sigmas, ncol, integer(1))
  total_q <- sum(q_terms)
  # byrow: row s holds replicate s's standard normals, laid out term by term.
  raw <- matrix(stats::rnorm(n_sim * total_q),
                nrow = n_sim, ncol = total_q, byrow = TRUE)

  col <- 0L
  terms <- vector("list", parts$n_terms)
  for (t in seq_len(parts$n_terms)) {
    qt <- q_terms[t]
    block <- raw[, (col + 1L):(col + qt), drop = FALSE]
    col <- col + qt
    terms[[t]] <- list(
      Z = as.matrix(parts$md$Z[[t]]),
      U = .apply_chol(block, parts$Sigmas[[t]])
    )
  }

  mom <- .mc_integrate_columns(eta_fixed, terms, inv_link, n_sim)
  # Standard error of the marginal (Monte Carlo) mean: sample SD / sqrt(n_sim).
  result <- list(fit = mom$mean,
                 se = sqrt(mom$m2 / (n_sim * pmax(n_sim - 1, 1))))

  if (se.fit) {
    return(list(fit = result$fit, se.fit = result$se))
  } else {
    return(result$fit)
  }
}
