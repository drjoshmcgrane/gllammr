#' Cluster-robust (sandwich) covariance for GLLAMM fixed effects
#'
#' Computes the cluster-robust covariance \eqn{H^{-1} M H^{-1}} where
#' \eqn{H} is the observed information of the outer parameters and
#' \eqn{M = \sum_j s_j s_j'} accumulates per-cluster score vectors. Because
#' clusters are independent, the Laplace marginal log-likelihood decomposes
#' by cluster; each cluster's score is evaluated by rebuilding its objective
#' on the cluster's rows and differentiating at the full-data estimates
#' (with the cluster's random effects re-profiled internally).
#'
#' @param object A two-level GLMM fitted via \code{gllamm()} (gaussian,
#'   binomial, poisson, or Gamma family, single random-effects term)
#'
#' @return Covariance matrix for the full outer parameter vector, with a
#'   \code{"fixed"} attribute holding the fixed-effects block.
#'
#' @seealso \code{\link{vcov.gllamm}} with \code{type = "sandwich"}
#' @keywords internal
sandwich_vcov_gllamm <- function(object) {
  if (is.null(object$formula) || is.null(object$data) || is.null(object$family)) {
    stop("Sandwich covariance requires a model fitted via gllamm()")
  }

  parsed <- parse_formula(object$formula, object$data)
  if (length(parsed$random_terms) != 1) {
    stop("Sandwich covariance currently supports a single random-effects term")
  }
  model_data <- make_model_matrices(parsed, object$data)
  groups <- model_data$groups[[1]]
  n_groups <- model_data$n_groups[1]

  theta_hat <- object$tmb_opt$par
  H <- tryCatch(
    optimHess(theta_hat, object$tmb_obj$fn, object$tmb_obj$gr),
    error = function(e) stop("Could not compute the observed information: ",
                             conditionMessage(e))
  )

  # Per-cluster scores at the full-data estimates
  scores <- matrix(0, n_groups, length(theta_hat))
  for (j in seq_len(n_groups)) {
    rows <- which(groups == (j - 1))
    d_j <- object$data[rows, , drop = FALSE]

    parsed_j <- parse_formula(object$formula, d_j)
    md_j <- make_model_matrices(parsed_j, d_j)
    fit_j <- try(
      fit_tmb_objective_only(md_j, object$family, parsed_j$random_terms),
      silent = TRUE
    )
    if (inherits(fit_j, "try-error")) {
      stop("Could not build the cluster ", j, " objective: ",
           attr(fit_j, "condition")$message)
    }
    # nll gradient -> score is its negative
    scores[j, ] <- -as.numeric(fit_j$gr(theta_hat))
  }

  meat <- crossprod(scores)
  H_inv <- solve(H)
  V <- H_inv %*% meat %*% H_inv
  dimnames(V) <- list(names(theta_hat), names(theta_hat))

  beta_idx <- which(names(theta_hat) == "beta")
  V_fixed <- V[beta_idx, beta_idx, drop = FALSE]
  dimnames(V_fixed) <- list(names(object$coefficients$fixed),
                            names(object$coefficients$fixed))
  attr(V, "fixed") <- V_fixed
  V
}


#' Build a TMB objective (no optimization) for one cluster
#'
#' Mirrors the single-term v2 engine's data/parameter construction so the
#' cluster objective is evaluated under exactly the same model.
#'
#' @keywords internal
fit_tmb_objective_only <- function(model_data, family, random_terms) {
  n_random <- model_data$n_random_coefs[1]
  correlated <- !random_terms[[1]]$uncorrelated

  use_slopes <- (n_random > 1)
  model_name <- if (use_slopes || family$family == "Gamma") {
    "glmm_slopes"
  } else if (family$family == "binomial") {
    "binomial"
  } else if (family$family == "poisson") {
    "poisson"
  } else {
    "gaussian"
  }

  Z_sparse <- Matrix::Matrix(model_data$Z[[1]], sparse = TRUE)
  tmb_data <- list(
    y = as.numeric(model_data$y),
    X = as.matrix(model_data$X),
    Z = if (model_name == "gaussian") as.matrix(model_data$Z[[1]]) else Z_sparse,
    groups = as.integer(model_data$groups[[1]]),
    n_groups = as.integer(model_data$n_groups[1]),
    n_obs = as.integer(model_data$n_obs),
    n_fixed = as.integer(model_data$n_fixed),
    n_random = as.integer(n_random),
    correlated = as.integer(correlated),
    weights = rep(1.0, model_data$n_obs),
    group_weights = rep(1.0, model_data$n_groups[1]),
    model_name = model_name
  )
  tmb_data$family <- switch(family$family,
                            gaussian = 0L, binomial = 1L, poisson = 2L,
                            Gamma = 3L, 0L)
  tmb_data$link <- if (family$family == "binomial") {
    switch(family$link, logit = 1L, probit = 2L, cloglog = 3L, 1L)
  } else if (family$family == "Gamma") {
    switch(family$link, log = 1L, inverse = 2L, identity = 3L, 1L)
  } else 1L

  n_theta <- max(1L, n_random * (n_random - 1) %/% 2)
  tmb_params <- switch(model_name,
    gaussian = list(beta = rep(0, model_data$n_fixed),
                    u = rep(0, tmb_data$n_groups * n_random),
                    log_sigma = 0, log_sigma_u = 0),
    glmm_slopes = list(beta = rep(0, model_data$n_fixed),
                       u = rep(0, tmb_data$n_groups * n_random),
                       log_sigma = 0,
                       log_sigma_u = rep(0, n_random),
                       theta = rep(0, n_theta)),
    list(beta = rep(0, model_data$n_fixed),
         u = rep(0, tmb_data$n_groups * n_random),
         log_sigma_u = rep(0, n_random),
         theta = rep(0, n_theta))
  )

  tmb_map <- list()
  has_theta <- !is.null(tmb_params$theta)
  theta_used <- use_slopes && correlated && n_random > 1
  if (has_theta && !theta_used) {
    tmb_map$theta <- factor(rep(NA, length(tmb_params$theta)))
  }
  if (model_name == "glmm_slopes" &&
      !(family$family %in% c("gaussian", "Gamma"))) {
    tmb_map$log_sigma <- factor(NA)
  }

  TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = "u",
    map = tmb_map,
    DLL = "GLLAMMR",
    silent = TRUE
  )
}
