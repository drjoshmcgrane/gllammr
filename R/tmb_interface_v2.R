#' Enhanced interface to TMB for GLLAMMR models
#'
#' Supports random intercepts, random slopes, and multiple GLM families
#'
#' @param model_data List from make_model_matrices()
#' @param family GLM family object
#' @param random_terms Parsed random effects terms
#' @param start_params Optional starting values
#' @param control Control parameters for optimization
#'
#' @return TMB fit object with parameter estimates
#' @keywords internal
fit_tmb_gllamm_v2 <- function(model_data, family, random_terms, start_params = NULL, control = list(), weights = NULL) {

  # For now, only handle single random effect term
  if (model_data$n_random_terms != 1) {
    stop("Currently only single random effects term supported")
  }

  # Determine model type
  n_random <- model_data$n_random_coefs[1]
  correlated <- !random_terms[[1]]$uncorrelated

  # Convert Z to sparse matrix for efficiency
  Z_sparse <- Matrix::Matrix(model_data$Z[[1]], sparse = TRUE)

  # Prepare weights vector (default to 1.0 if NULL)
  if (is.null(weights)) {
    weights_vec <- rep(1.0, model_data$n_obs)
  } else {
    weights_vec <- as.numeric(weights)
  }

  # Select model template: family + random-effects structure
  use_slopes <- (n_random > 1)
  model_name <- if (use_slopes) {
    "glmm_slopes"
  } else if (family$family == "binomial") {
    "binomial"
  } else if (family$family == "poisson") {
    "poisson"
  } else {
    "gaussian"
  }

  # Prepare TMB data. The plain gaussian template expects a dense Z;
  # the others declare DATA_SPARSE_MATRIX(Z).
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
    weights = weights_vec,
    model_name = model_name
  )

  # Family and link codes (read by the binomial and glmm_slopes templates;
  # unused data entries are ignored by the others)
  tmb_data$family <- switch(family$family,
                            gaussian = 0L, binomial = 1L, poisson = 2L, 0L)
  tmb_data$link <- if (family$family == "binomial") {
    switch(family$link, logit = 1L, probit = 2L, cloglog = 3L, 1L)
  } else {
    1L  # canonical link (identity / log)
  }

  # Initialize parameters
  if (is.null(start_params)) {
    # Fit simple model for initial values
    if (family$family == "gaussian") {
      lm_fit <- lm(model_data$y ~ model_data$X - 1)
      beta_init <- coef(lm_fit)
      sigma_init <- summary(lm_fit)$sigma
    } else {
      glm_fit <- glm(model_data$y ~ model_data$X - 1, family = family)
      beta_init <- coef(glm_fit)
      sigma_init <- 1.0  # Not used for non-Gaussian
    }

    # Initialize random effects to zero
    u_init <- rep(0, tmb_data$n_groups * n_random)

    # Initialize variance components
    log_sigma_u_init <- rep(log(0.5), n_random)

    # Initialize correlation parameters (Cholesky)
    n_theta <- n_random * (n_random - 1) / 2
    theta_init <- rep(0, max(n_theta, 1))

    # Each template reads an exact parameter set; passing extras would add
    # dead entries to the optimization vector.
    tmb_params <- switch(model_name,
      gaussian = list(
        beta = beta_init,
        u = u_init,
        log_sigma = log(max(sigma_init, 0.1)),
        log_sigma_u = log_sigma_u_init[1]
      ),
      glmm_slopes = list(
        beta = beta_init,
        u = u_init,
        log_sigma = log(max(sigma_init, 0.1)),
        log_sigma_u = log_sigma_u_init,
        theta = theta_init
      ),
      # binomial and poisson: no residual SD parameter
      list(
        beta = beta_init,
        u = u_init,
        log_sigma_u = log_sigma_u_init,
        theta = theta_init
      )
    )
  } else {
    tmb_params <- start_params
  }

  # Fix parameters the chosen model never reads: theta when there is no
  # correlation structure, log_sigma for non-gaussian slopes models. Leaving
  # them free creates flat likelihood directions and singular Hessians.
  tmb_map <- list()
  has_theta <- !is.null(tmb_params$theta)
  theta_used <- use_slopes && correlated && n_random > 1
  if (has_theta && !theta_used) {
    tmb_map$theta <- factor(rep(NA, length(tmb_params$theta)))
  }
  if (model_name == "glmm_slopes" && family$family != "gaussian") {
    tmb_map$log_sigma <- factor(NA)
  }

  # Create TMB object
  tryCatch({
    obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = tmb_params,
      random = "u",
      map = tmb_map,
      DLL = "GLLAMMR",
      silent = TRUE
    )
  }, error = function(e) {
    stop("Failed to create TMB object: ", e$message)
  })

  # Optimize
  control_defaults <- list(
    eval.max = 1000,
    iter.max = 500,
    trace = 0
  )
  control <- modifyList(control_defaults, control)

  opt <- try(
    nlminb(
      start = obj$par,
      objective = obj$fn,
      gradient = obj$gr,
      control = control
    ),
    silent = FALSE
  )

  if (inherits(opt, "try-error")) {
    stop("Optimization failed")
  }

  # Check convergence
  converged <- (opt$convergence == 0)

  # Get standard errors
  sdr <- try(TMB::sdreport(obj), silent = TRUE)

  if (inherits(sdr, "try-error")) {
    warning("Failed to compute standard errors")
    sdr <- NULL
  }

  # Extract results
  par_full <- obj$env$last.par.best

  # Fixed effects
  beta_hat <- par_full[names(par_full) == "beta"]
  names(beta_hat) <- colnames(model_data$X)

  # Variance components
  log_sigma_u_hat <- par_full[names(par_full) == "log_sigma_u"]
  sigma_u_hat <- exp(log_sigma_u_hat)

  # Build variance-covariance matrix for random effects
  if (n_random > 1 && correlated) {
    theta_hat <- par_full[names(par_full) == "theta"]

    # Reconstruct Cholesky factor
    L <- matrix(0, n_random, n_random)
    diag(L) <- 1
    idx <- 1
    for (i in 2:n_random) {
      for (j in 1:(i-1)) {
        L[i, j] <- theta_hat[idx]
        idx <- idx + 1
      }
    }

    # Normalized correlation matrix (matches the template: L L' rescaled to
    # unit diagonal so sigma_u are genuine standard deviations)
    R <- L %*% t(L)
    R <- R / sqrt(diag(R) %o% diag(R))

    # Covariance matrix
    D <- diag(sigma_u_hat)
    Sigma_u <- D %*% R %*% D
  } else {
    # Diagonal covariance (nrow guards the scalar case: diag(x) with a
    # length-1 numeric would otherwise build a floor(x)-dimensional identity)
    Sigma_u <- diag(sigma_u_hat^2, nrow = n_random)
  }

  # Random effects
  u_hat <- par_full[names(par_full) == "u"]

  # Organize random effects by group
  random_effects <- list()
  for (g in 1:tmb_data$n_groups) {
    idx_start <- (g - 1) * n_random + 1
    idx_end <- g * n_random
    random_effects[[g]] <- u_hat[idx_start:idx_end]
  }

  # Standard errors
  if (!is.null(sdr)) {
    se_all <- summary(sdr, "fixed")
    se_beta <- se_all[rownames(se_all) == "beta", "Std. Error"]
    vcov_fixed <- matrix(0, length(beta_hat), length(beta_hat))
    diag(vcov_fixed) <- se_beta^2
    dimnames(vcov_fixed) <- list(names(beta_hat), names(beta_hat))
  } else {
    vcov_fixed <- matrix(NA, length(beta_hat), length(beta_hat))
  }

  # Fitted values
  fitted_vals <- as.numeric(model_data$X %*% beta_hat)
  for (i in 1:length(fitted_vals)) {
    g <- tmb_data$groups[i] + 1
    fitted_vals[i] <- fitted_vals[i] + sum(model_data$Z[[1]][i, ] * random_effects[[g]])
  }

  # Apply inverse link for GLMs
  if (family$family != "gaussian") {
    fitted_vals <- family$linkinv(fitted_vals)
  }

  # Log-likelihood
  loglik <- -opt$objective

  # Number of parameters
  n_params <- length(beta_hat) + n_random + ifelse(family$family == "gaussian", 1, 0)
  if (n_random > 1 && correlated) {
    n_params <- n_params + n_random * (n_random - 1) / 2
  }

  list(
    coefficients = list(
      fixed = beta_hat,
      random_var = list(Sigma_u)
    ),
    vcov = list(
      fixed = vcov_fixed,
      all = NULL
    ),
    random_effects = random_effects,
    fitted.values = fitted_vals,
    logLik = loglik,
    n_params = n_params,
    convergence = list(
      converged = converged,
      message = opt$message,
      iterations = opt$iterations
    ),
    tmb_obj = obj,
    tmb_opt = opt,
    tmb_sdr = sdr
  )
}
