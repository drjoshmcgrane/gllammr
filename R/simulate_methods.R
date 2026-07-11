# simulate() methods for the model classes that previously fell through
# to the formula-based generic (which cannot handle them), plus
# polytomous IRT prediction. All follow the stats::simulate contract
# (fresh population draws, "seed" attribute).

.sim_rng_setup <- function(seed) {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    runif(1)
  }
  if (is.null(seed)) {
    get(".Random.seed", envir = globalenv())
  } else {
    set.seed(seed)
    structure(seed, kind = as.list(RNGkind()))
  }
}

.sim_wrap <- function(sims, rng_state) {
  out <- as.data.frame(sims)
  names(out) <- paste0("sim_", seq_len(ncol(out)))
  attr(out, "seed") <- rng_state
  out
}

# Total latent SD for an IRT/EIRT fit: person deviation plus all group
# random-effect components (independent normals, so variances add). This
# is the correct integration SD for marginal (population-averaged)
# predictions from multilevel fits.
.irt_total_latent_sd <- function(object) {
  v <- (object$ability_sd %||% 1)^2
  re <- object$random_effects
  if (!is.null(re$sigma_random)) {
    v <- v + sum(re$sigma_random^2)
  }
  sqrt(v)
}

# Draw one fresh composite-ability vector for an IRT/EIRT fit: person
# deviation plus, on multilevel fits, fresh group effects for EVERY
# random-effect term (using the grouping structure stored on the fit).
.irt_sim_theta <- function(object) {
  n <- object$n_persons
  theta <- rnorm(n, 0, object$ability_sd %||% 1)
  re <- object$random_effects
  if (!is.null(re$sigma_random) && !is.null(re$group_ids)) {
    gid <- as.matrix(re$group_ids)
    des <- as.matrix(re$re_design)
    for (k in seq_along(re$sigma_random)) {
      u <- rnorm(re$n_groups[k], 0, re$sigma_random[k])
      ok <- gid[, k] >= 0  # -1 marks NA (partial nesting)
      theta[ok] <- theta[ok] + des[ok, k] * u[gid[ok, k] + 1L]
    }
  }
  theta
}


#' Simulate from a fitted multinomial model
#'
#' Draws fresh random effects for every term and samples categories from
#' the baseline-category logit probabilities.
#'
#' @param object Fitted \code{gllamm_multinomial} object
#' @param nsim Number of replicates
#' @param seed Optional seed (stored as the \code{"seed"} attribute)
#' @param newdata Optional new data frame
#' @param ... Unused
#' @export
simulate.gllamm_multinomial <- function(object, nsim = 1, seed = NULL,
                                        newdata = NULL, ...) {
  rng_state <- .sim_rng_setup(seed)
  sim_data <- if (is.null(newdata)) object$data else newdata
  parsed <- parse_formula(object$formula, sim_data)
  md <- make_model_matrices(parsed, sim_data)
  K <- object$n_categories
  beta <- object$coefficients$beta              # (K-1) x p
  eta_fixed <- as.matrix(md$X) %*% t(beta)      # n x (K-1)
  n_obs <- nrow(eta_fixed)

  rv <- object$coefficients$random_var
  Sigmas <- if (is.list(rv)) lapply(rv, as.matrix) else
    list(diag(as.numeric(rv), nrow = ncol(md$Z[[1]])))

  sims <- matrix(NA_integer_, n_obs, nsim)
  for (s in seq_len(nsim)) {
    shift <- rep(0, n_obs)
    for (t in seq_along(Sigmas)) {
      u_t <- rmvnorm_chol(md$n_groups[t], Sigmas[[t]])
      shift <- shift + rowSums(md$Z[[t]] *
                                 u_t[md$groups[[t]] + 1, , drop = FALSE])
    }
    expeta <- cbind(1, exp(eta_fixed + shift))
    probs <- expeta / rowSums(expeta)
    cum <- t(apply(probs, 1, cumsum))
    r <- runif(n_obs)
    sims[, s] <- 1L + rowSums(r > cum[, -K, drop = FALSE])
  }
  out <- .sim_wrap(sims, rng_state)
  # report as category labels when available
  if (!is.null(object$categories)) {
    out[] <- lapply(out, function(col) object$categories[col])
    attr(out, "seed") <- rng_state
  }
  out
}


#' Simulate from a fitted parametric frailty survival model
#'
#' Draws fresh frailties per replicate and generates uncensored event
#' times by inverse transform from the fitted exponential or Weibull
#' hazard (censoring schemes are design choices, so simulated times are
#' uncensored).
#'
#' @inheritParams simulate.gllamm_multinomial
#' @param object Fitted \code{gllamm_survival} object
#' @export
simulate.gllamm_survival <- function(object, nsim = 1, seed = NULL,
                                     newdata = NULL, ...) {
  rng_state <- .sim_rng_setup(seed)
  sim_data <- if (is.null(newdata)) object$data else newdata

  rhs <- deparse(object$formula[[3]], width.cutoff = 500)
  design_formula <- stats::as.formula(
    paste(object$time_var, "~", paste(rhs, collapse = " ")))
  parsed <- parse_formula(design_formula, sim_data)
  md <- make_model_matrices(parsed, sim_data)
  X <- as.matrix(md$X)
  n_obs <- nrow(X)
  beta <- object$coefficients$fixed
  shape <- if (identical(object$distribution, "weibull")) {
    object$shape
  } else {
    1
  }
  sigma_u <- object$coefficients$random_sd
  eta_fixed <- as.numeric(X %*% beta)

  sims <- matrix(NA_real_, n_obs, nsim)
  for (s in seq_len(nsim)) {
    u <- rnorm(md$n_groups[1], 0, sigma_u)
    eta <- eta_fixed + u[md$groups[[1]] + 1]
    # Template parameterization: S(t) = exp(-(exp(eta) * t)^shape)
    #   =>  T = E^(1/shape) / exp(eta),  E ~ Exp(1)
    sims[, s] <- rexp(n_obs)^(1 / shape) / exp(eta)
  }
  .sim_wrap(sims, rng_state)
}


#' Simulate response matrices from a fitted IRT model
#'
#' Draws fresh person abilities from the fitted ability distribution and
#' samples item responses - dichotomous or polytomous - from the
#' estimated item parameters (parametric bootstrap / posterior
#' predictive style). Works for both EM and Laplace fits.
#'
#' @inheritParams simulate.gllamm_multinomial
#' @param object Fitted \code{gllamm_irt} object
#' @return A list of \code{nsim} simulated response matrices, with a
#'   \code{"seed"} attribute.
#' @export
simulate.gllamm_irt <- function(object, nsim = 1, seed = NULL, ...) {
  rng_state <- .sim_rng_setup(seed)
  n <- object$n_persons
  J <- object$n_items
  sd_th <- object$ability_sd %||% 1

  poly <- inherits(object, "gllamm_irt_poly") ||
    isTRUE(object$max_categories > 2)
  out <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    theta <- .irt_sim_theta(object)
    Y <- matrix(NA_integer_, n, J)
    if (!poly) {
      b <- object$item_parameters$difficulty
      a <- object$item_parameters$discrimination %||% rep(1, J)
      g <- object$item_parameters$guessing %||% rep(0, J)
      a[is.na(a)] <- 1
      g[is.na(g)] <- 0
      for (j in seq_len(J)) {
        p <- g[j] + (1 - g[j]) * plogis(a[j] * (theta - b[j]))
        Y[, j] <- rbinom(n, 1, p)
      }
    } else {
      model_name <- if (object$model == "PCM") "GPCM" else object$model
      for (j in seq_len(J)) {
        th_j <- object$item_parameters$thresholds[[j]]
        a_j <- object$item_parameters$discrimination[j]
        pr <- irt_category_probs(model_name, theta, th_j, a_j)
        cum <- t(apply(pr, 1, cumsum))
        r <- runif(n)
        Y[, j] <- 1L + rowSums(r > cum[, -ncol(cum), drop = FALSE])
      }
    }
    out[[s]] <- Y
  }
  attr(out, "seed") <- rng_state
  out
}


#' Simulate response matrices from a fitted latent class model
#'
#' Draws class memberships from the estimated prevalences and item
#' responses from the class-conditional distributions (binary,
#' categorical, and gaussian indicators).
#'
#' @inheritParams simulate.gllamm_irt
#' @param object Fitted \code{gllamm_lca} object
#' @export
simulate.gllamm_lca <- function(object, nsim = 1, seed = NULL, ...) {
  rng_state <- .sim_rng_setup(seed)
  n <- object$n_obs
  J <- object$n_items
  K <- object$nclass
  out <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    cls <- sample.int(K, n, replace = TRUE, prob = object$class_probs)
    Y <- matrix(NA_real_, n, J)
    for (j in seq_len(J)) {
      tp <- object$item_type[j]
      if (tp == 0L) {
        Y[, j] <- rbinom(n, 1, object$item_probs[j, cls])
      } else if (tp == 1L) {
        pj <- object$cat_probs[[rownames(object$item_probs)[j]]]
        Y[, j] <- vapply(cls, function(k) {
          sample.int(nrow(pj), 1, prob = pj[, k])
        }, 1L)
      } else {
        gp <- object$gaussian_params
        gi <- match(rownames(object$item_probs)[j], rownames(gp$means))
        Y[, j] <- rnorm(n, gp$means[gi, cls], gp$sds[gi, cls])
      }
    }
    out[[s]] <- Y
  }
  attr(out, "seed") <- rng_state
  out
}


#' Simulate response matrices from a fitted cognitive diagnosis model
#'
#' Draws attribute profiles from the estimated prevalences and item
#' responses from the per-profile item kernels.
#'
#' @inheritParams simulate.gllamm_irt
#' @param object Fitted \code{gllamm_cdm} object
#' @export
simulate.gllamm_cdm <- function(object, nsim = 1, seed = NULL, ...) {
  rng_state <- .sim_rng_setup(seed)
  n <- object$n_obs
  J <- object$n_items
  profiles <- object$profiles
  out <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    k <- sample.int(nrow(profiles), n, replace = TRUE,
                    prob = object$profile_probs)
    Y <- matrix(NA_integer_, n, J)
    for (j in seq_len(J)) {
      e <- object$item_params[[j]]
      meas <- match(e$measured, colnames(object$Q))
      red <- profiles[k, meas, drop = FALSE]
      key <- apply(red, 1, paste, collapse = "")
      p <- e$prob[key]
      Y[, j] <- rbinom(n, 1, p)
    }
    out[[s]] <- Y
  }
  attr(out, "seed") <- rng_state
  out
}


#' Simulate indicator data from a fitted SEM
#'
#' Draws from the fitted multivariate-normal implied distribution
#' (means + implied covariance) of the indicators.
#'
#' @inheritParams simulate.gllamm_irt
#' @param object Fitted \code{gllamm_sem} object (ML method)
#' @export
simulate.gllamm_sem <- function(object, nsim = 1, seed = NULL, ...) {
  rng_state <- .sim_rng_setup(seed)
  if (is.null(object$latent_covariance)) {
    stop("simulate() requires an ML-method SEM fit")
  }
  ind <- names(object$indicator_residual_sd)
  Lambda <- object$loadings
  V_eta <- object$latent_variance_total[colnames(Lambda), colnames(Lambda)]
  Sigma <- Lambda %*% V_eta %*% t(Lambda) +
    diag(object$indicator_residual_sd^2, length(ind))
  mu <- object$intercepts[ind]
  n <- object$n_obs
  out <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    Y <- rmvnorm_chol(n, Sigma)
    Y <- sweep(Y, 2, mu, "+")
    colnames(Y) <- ind
    out[[s]] <- Y
  }
  attr(out, "seed") <- rng_state
  out
}


#' Predicted category probabilities for polytomous IRT fits
#'
#' Returns the model-implied category probabilities at each person's
#' estimated ability: a list (one element per item) of persons x
#' categories matrices, or the persons x items matrix of expected scores
#' with \code{type = "expected"}.
#'
#' @param object Fitted polytomous \code{gllamm_irt} model
#' @param type "probs" (default) or "expected"
#' @param ... Unused
#' @export
predict.gllamm_irt_poly <- function(object,
                                    type = c("probs", "expected",
                                             "marginal", "probability",
                                             "ability"),
                                    ...) {
  warn_not_converged(object)
  type <- match.arg(type)
  # Marginal/ability/probability predictions are handled by the general
  # IRT method (Monte Carlo over the ability distribution etc.)
  if (type %in% c("marginal", "probability", "ability")) {
    return(NextMethod())
  }
  theta <- object$person_abilities
  J <- object$n_items
  model_name <- if (object$model == "PCM") "GPCM" else object$model
  probs <- lapply(seq_len(J), function(j) {
    irt_category_probs(model_name, theta,
                       object$item_parameters$thresholds[[j]],
                       object$item_parameters$discrimination[j])
  })
  names(probs) <- names(object$item_parameters$thresholds)
  if (type == "probs") return(probs)
  vapply(probs, function(p) as.numeric(p %*% seq_len(ncol(p))),
         numeric(length(theta)))
}


#' Simulate response matrices from a fitted explanatory IRT model
#'
#' Draws fresh person abilities (plus group effects when present is not
#' supported - simulation is at the population level) and samples
#' responses from the estimated item parameters.
#'
#' @inheritParams simulate.gllamm_irt
#' @param object Fitted \code{gllamm_eirt} object
#' @export
simulate.gllamm_eirt <- function(object, nsim = 1, seed = NULL, ...) {
  rng_state <- .sim_rng_setup(seed)
  n <- object$n_persons
  J <- object$n_items
  sd_th <- object$ability_sd %||% 1
  out <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    theta <- .irt_sim_theta(object)
    Y <- matrix(NA_integer_, n, J)
    if (!isTRUE(object$is_polytomous)) {
      b <- object$item_parameters$difficulty
      a <- object$item_parameters$discrimination
      if (is.null(a) || all(is.na(a))) a <- rep(1, J)
      for (j in seq_len(J)) {
        Y[, j] <- rbinom(n, 1, plogis(a[j] * (theta - b[j])))
      }
    } else {
      th_list <- eirt_item_thresholds(object)
      model_name <- eirt_poly_model_name(object$poly_model_type)
      a <- object$item_parameters$discrimination
      if (is.null(a) || all(is.na(a))) a <- rep(1, J)
      for (j in seq_len(J)) {
        pr <- irt_category_probs(model_name, theta, th_list[[j]], a[j])
        cum <- t(apply(pr, 1, cumsum))
        r <- runif(n)
        Y[, j] <- 1L + rowSums(r > cum[, -ncol(cum), drop = FALSE])
      }
    }
    out[[s]] <- Y
  }
  attr(out, "seed") <- rng_state
  out
}


#' Simulate from a fitted NPML model
#'
#' Draws group-level intercepts from the estimated discrete (mass-point)
#' distribution and responses from the fitted family.
#'
#' @inheritParams simulate.gllamm_multinomial
#' @param object Fitted \code{gllamm_npml} object
#' @export
simulate.gllamm_npml <- function(object, nsim = 1, seed = NULL,
                                 newdata = NULL, ...) {
  rng_state <- .sim_rng_setup(seed)
  sim_data <- if (is.null(newdata)) object$data else newdata
  parsed <- parse_formula(object$formula, sim_data)
  md <- make_model_matrices(parsed, sim_data)
  X <- as.matrix(md$X)
  beta <- object$coefficients$fixed
  # The mass-point locations play the role of group intercepts; drop the
  # fixed intercept column contribution if locations carry it
  has_int <- "(Intercept)" %in% colnames(X)
  eta_fixed <- as.numeric(X[, setdiff(colnames(X), "(Intercept)"),
                            drop = FALSE] %*%
                            beta[setdiff(names(beta), "(Intercept)")])
  fam <- object$family$family %||% "gaussian"
  n_obs <- nrow(X)

  sims <- matrix(NA_real_, n_obs, nsim)
  for (s in seq_len(nsim)) {
    k <- sample.int(length(object$masses), md$n_groups[1],
                    replace = TRUE, prob = object$masses)
    eta <- eta_fixed + object$locations[k][md$groups[[1]] + 1]
    sims[, s] <- switch(fam,
      gaussian = rnorm(n_obs, eta, object$sigma %||% 1),
      binomial = rbinom(n_obs, 1, plogis(eta)),
      poisson = rpois(n_obs, exp(eta)),
      stop("simulate not implemented for NPML family: ", fam))
  }
  .sim_wrap(sims, rng_state)
}


#' Simulate from a fitted mixed-response model
#'
#' Draws fresh shared random intercepts and simulates every outcome from
#' its fitted family (gaussian/binomial/poisson), returning a list of
#' data frames.
#'
#' @inheritParams simulate.gllamm_multinomial
#' @param object Fitted \code{gllamm_mixed} object
#' @export
simulate.gllamm_mixed <- function(object, nsim = 1, seed = NULL, ...) {
  rng_state <- .sim_rng_setup(seed)
  data <- object$data
  re_term <- attr(stats::terms(object$random), "term.labels")
  rt <- parse_random_term(re_term, data)
  g <- factor(data[[rt$grouping[1]]])
  n_groups <- nlevels(g)
  gi <- as.integer(g)
  sigma_u <- object$random_sd
  loadings <- object$loadings %||%
    setNames(rep(1, length(object$outcomes)), object$outcomes)
  out <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    u <- rnorm(n_groups, 0, sigma_u)
    sim_d <- list()
    for (fam in object$outcomes) {
      f <- object$formulas[[fam]]
      X <- model.matrix(f, data = data)
      cf <- object$coefficients[[fam]]
      lam <- loadings[[fam]] %||% 1
      eta <- as.numeric(X %*% cf) + lam * u[gi]
      yname <- as.character(f[[2]])
      sim_d[[yname]] <- switch(fam,
        gaussian = rnorm(length(eta), eta, object$residual_sd %||% 1),
        binomial = rbinom(length(eta), 1, plogis(eta)),
        poisson = rpois(length(eta), exp(eta)))
    }
    out[[s]] <- as.data.frame(sim_d)
  }
  attr(out, "seed") <- rng_state
  out
}


#' @export
simulate.gllamm_rank <- function(object, nsim = 1, seed = NULL, ...) {
  stop("simulate() is not implemented for rank-ordered logit models; ",
       "rankings can be simulated by drawing Gumbel utilities from the ",
       "fitted linear predictor (see ?fit_rank)")
}
