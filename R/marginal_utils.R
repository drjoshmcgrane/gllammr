#' Utility Functions for Marginal Predictions
#'
#' Internal functions for computing population-averaged (marginal) predictions
#' via Monte Carlo integration over random effects distribution
#'
#' @keywords internal
#' @name marginal_utils
NULL


#' Monte Carlo integration for marginal predictions
#'
#' Computes E[g^(-1)(X'beta + Z'u)] by averaging over draws from u ~ N(0, Sigma_u)
#'
#' @param X Fixed effects design matrix (n x p)
#' @param Z Random effects design matrix (n x q) for a single observation
#' @param beta Fixed effects coefficients (p x 1)
#' @param u_samples Matrix of random effects samples (n_sim x q)
#' @param inv_link_fn Inverse link function
#'
#' @return Vector of marginal predictions (length n)
#' @keywords internal
mc_integrate_fixed_samples <- function(X, Z, beta, u_samples, inv_link_fn) {
  n_obs <- nrow(X)
  n_sim <- nrow(u_samples)

  # Compute fixed part (same for all samples)
  eta_fixed <- as.vector(X %*% beta)

  # Running mean and variance (Welford) instead of an n_obs x n_sim matrix
  mean_pred <- numeric(n_obs)
  m2_pred <- numeric(n_obs)

  for (s in 1:n_sim) {
    u <- u_samples[s, ]

    # Random part for all observations
    if (is.matrix(Z)) {
      eta_random <- as.vector(Z %*% u)
    } else {
      # Z is a vector (single random effect)
      eta_random <- Z * u
    }

    pred_s <- inv_link_fn(eta_fixed + eta_random)

    delta <- pred_s - mean_pred
    mean_pred <- mean_pred + delta / s
    m2_pred <- m2_pred + delta * (pred_s - mean_pred)
  }

  list(
    fit = mean_pred,
    se = sqrt(m2_pred / (n_sim - 1))
  )
}


#' Draw samples from a zero-mean multivariate normal
#'
#' Cholesky-based sampler so MASS can remain in Suggests.
#'
#' @param n Number of samples
#' @param Sigma Variance-covariance matrix
#' @return n x ncol(Sigma) matrix of draws
#' @keywords internal
rmvnorm_chol <- function(n, Sigma) {
  q <- ncol(Sigma)
  if (q == 1) {
    v <- Sigma[1, 1]
    if (v < 0) {
      warning("Random-effects variance is negative; using zero.",
              call. = FALSE)
      v <- 0
    }
    return(matrix(rnorm(n, 0, sqrt(max(v, 0))), n, 1))
  }
  L <- safe_chol(Sigma)  # upper triangular: Sigma = t(L) %*% L
  matrix(rnorm(n * q), n, q) %*% L
}


#' Cholesky factor with a positive-definite fallback
#'
#' Returns the upper-triangular Cholesky factor of \code{Sigma}. If
#' \code{Sigma} is not positive definite, warns and retries on the nearest
#' positive-definite matrix (\code{Matrix::nearPD}); if that also fails,
#' stops with an informative message.
#'
#' @param Sigma A covariance matrix.
#' @return Upper-triangular Cholesky factor such that \code{t(L) \%*\% L == Sigma}.
#' @keywords internal
#' @noRd
safe_chol <- function(Sigma) {
  L <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (!is.null(L)) return(L)
  warning("Random-effects covariance is not positive definite; using ",
          "nearest positive-definite approximation.", call. = FALSE)
  Sigma_pd <- tryCatch(
    as.matrix(Matrix::nearPD(Sigma)$mat),
    error = function(e) NULL
  )
  if (!is.null(Sigma_pd)) {
    L <- tryCatch(chol(Sigma_pd), error = function(e) NULL)
    if (!is.null(L)) return(L)
  }
  stop("Random-effects covariance is not positive definite and could not ",
       "be repaired; cannot draw random-effects samples.", call. = FALSE)
}


#' Monte Carlo integration for marginal predictions (one sample at a time)
#'
#' More memory-efficient version that processes samples sequentially
#'
#' @param X Fixed effects design matrix
#' @param Z Random effects design matrix
#' @param beta Fixed effects coefficients
#' @param Sigma_u Random effects variance-covariance matrix
#' @param inv_link_fn Inverse link function
#' @param n_sim Number of Monte Carlo samples
#'
#' @return List with fit and se
#' @keywords internal
mc_integrate_marginal <- function(X, Z, beta, Sigma_u,
                                  inv_link_fn, n_sim = 1000) {
  n_obs <- nrow(X)
  n_random <- ncol(Sigma_u)

  # Draw samples from random effects distribution
  u_samples <- rmvnorm_chol(n_sim, Sigma_u)

  # Use fixed samples version
  mc_integrate_fixed_samples(X, Z, beta, u_samples, inv_link_fn)
}


#' Get inverse link function for a family
#'
#' @param family Family object (standard R family or custom)
#' @return Inverse link function
#' @keywords internal
get_inverse_link <- function(family) {
  # Handle custom family classes
  if (inherits(family, "binomial_family")) {
    return(switch(family$link,
      logit = plogis,
      probit = pnorm,
      cloglog = function(x) 1 - exp(-exp(x)),
      stop("Unknown binomial link:", family$link)
    ))
  }

  if (inherits(family, "ordinal_family")) {
    # Ordinal families need special handling (cumulative probabilities)
    # Will be handled in predict.gllamm_ordinal
    stop("Use predict.gllamm with ordinal_family, not direct inverse link")
  }

  # Standard R families
  if (family$family == "gaussian") {
    return(identity)
  }

  if (family$family == "binomial") {
    return(switch(family$link,
      logit = plogis,
      probit = pnorm,
      cloglog = function(x) 1 - exp(-exp(x)),
      family$linkinv
    ))
  }

  if (family$family == "poisson") {
    return(exp)  # Log link
  }

  if (family$family == "Gamma") {
    return(switch(family$link,
      log = exp,
      inverse = function(x) 1/x,
      identity = identity,
      family$linkinv
    ))
  }

  # Fallback to family's linkinv
  if (!is.null(family$linkinv)) {
    return(family$linkinv)
  }

  stop("Cannot determine inverse link for family:", family$family)
}


#' Extract random effects variance-covariance matrix from fitted model
#'
#' @param object Fitted gllamm object
#' @return Variance-covariance matrix of random effects
#' @keywords internal
extract_random_vcov <- function(object) {
  # Try to extract from TMB sdreport
  if (!is.null(object$tmb_sdr)) {
    sdr <- object$tmb_sdr

    # Look for sigma_u in parameter estimates
    par_names <- names(sdr$par.fixed)

    if ("sigma_u" %in% par_names) {
      # Multiple random effects (vector of SDs)
      sigma_u <- sdr$par.fixed[par_names == "sigma_u"]
      n_random <- length(sigma_u)

      # Check for correlation parameters (theta)
      if ("theta" %in% par_names && n_random > 1) {
        theta <- sdr$par.fixed[par_names == "theta"]

        # Reconstruct correlation matrix from Cholesky parameters
        L <- matrix(0, n_random, n_random)
        diag(L) <- 1

        idx <- 1
        for (i in 2:n_random) {
          for (j in 1:(i-1)) {
            L[i, j] <- theta[idx]
            idx <- idx + 1
          }
        }

        R <- L %*% t(L)
        R <- R / sqrt(diag(R) %o% diag(R))  # normalize to unit diagonal

        # Convert to covariance matrix
        Sigma_u <- outer(sigma_u, sigma_u) * R
        return(Sigma_u)
      } else {
        # Uncorrelated random effects (diagonal)
        return(diag(sigma_u^2, nrow = n_random))
      }
    }

    if ("log_sigma_u" %in% par_names) {
      log_sigma_u <- sdr$par.fixed[par_names == "log_sigma_u"]
      sigma_u <- exp(log_sigma_u)
      n_random <- length(sigma_u)

      # Check for theta
      if ("theta" %in% par_names && n_random > 1) {
        theta <- sdr$par.fixed[par_names == "theta"]

        L <- matrix(0, n_random, n_random)
        diag(L) <- 1

        idx <- 1
        for (i in 2:n_random) {
          for (j in 1:(i-1)) {
            L[i, j] <- theta[idx]
            idx <- idx + 1
          }
        }

        R <- L %*% t(L)
        R <- R / sqrt(diag(R) %o% diag(R))  # normalize to unit diagonal
        Sigma_u <- outer(sigma_u, sigma_u) * R
        return(Sigma_u)
      } else {
        return(diag(sigma_u^2, nrow = n_random))
      }
    }
  }

  # Fallback: try coefficients$random_var (for older interface)
  if (!is.null(object$coefficients$random_var)) {
    random_var <- object$coefficients$random_var

    if (is.list(random_var)) {
      # Multiple random effect terms
      # For now, assume we want the first term
      var_mat <- random_var[[1]]

      if (is.matrix(var_mat)) {
        return(var_mat)
      } else {
        # Scalar variance
        return(matrix(var_mat, 1, 1))
      }
    } else {
      # Single scalar variance
      return(matrix(random_var, 1, 1))
    }
  }

  stop("Cannot extract random effects variance from object. ",
       "Make sure the model was fit with random effects.")
}


#' Construct random effects design matrix for newdata
#'
#' @param newdata New data frame
#' @param random_terms Parsed random effects terms from original model
#' @param group_var Name of grouping variable
#'
#' @return Random effects design matrix Z
#' @keywords internal
construct_Z_matrix <- function(newdata, random_terms, group_var = NULL) {
  if (is.null(random_terms) || length(random_terms) == 0) {
    stop("No random effects in model")
  }

  # Same construction the fitting code uses (make_model_matrices), so the
  # columns line up with the estimated random effects: intercept-only gives
  # a column of 1s, slope terms add their covariate columns.
  term <- random_terms[[1]]
  model.matrix(term$formula, data = newdata)
}
