#' Interface to TMB for GLLAMM models
#'
#' Prepares data and parameters for TMB, compiles and runs optimization
#'
#' @param model_data List from make_model_matrices()
#' @param family GLM family object
#' @param start_params Optional starting values
#' @param control Control parameters for optimization
#'
#' @return TMB fit object with parameter estimates
#' @keywords internal
fit_tmb_gllamm <- function(model_data, family, start_params = NULL, control = list(), weights = NULL) {

  # Currently only support Gaussian
  if (family$family != "gaussian" || family$link != "identity") {
    stop("Only Gaussian family with identity link currently supported")
  }

  # For now, only handle single random effect term
  if (model_data$n_random_terms != 1) {
    stop("Currently only single random effects term supported")
  }

  # Determine which template to use based on number of random coefficients
  n_random <- model_data$n_random_coefs[1]
  use_slopes <- (n_random > 1)

  # Prepare weights vector (default to 1.0 if NULL)
  if (is.null(weights)) {
    weights_vec <- rep(1.0, model_data$n_obs)
  } else {
    weights_vec <- as.numeric(weights)
  }

  # Prepare TMB data
  tmb_data <- list(
    y = as.numeric(model_data$y),
    X = as.matrix(model_data$X),
    Z = as.matrix(model_data$Z[[1]]),
    groups = as.integer(model_data$groups[[1]]),
    n_groups = as.integer(model_data$n_groups[1]),
    n_obs = as.integer(model_data$n_obs),
    n_fixed = as.integer(model_data$n_fixed),
    n_random = as.integer(model_data$n_random_coefs[1]),
    weights = weights_vec,
    model_name = "gaussian"
  )

  # Initialize parameters
  if (is.null(start_params)) {
    # Fit simple lm for initial values
    lm_fit <- lm(model_data$y ~ model_data$X - 1)

    beta_init <- coef(lm_fit)
    sigma_init <- summary(lm_fit)$sigma

    # Initialize random effects to zero
    u_init <- rep(0, tmb_data$n_groups * tmb_data$n_random)

    # Initialize variance components
    # Use 50% of residual variance for random effect
    sigma_u_init <- sigma_init * 0.7

    tmb_params <- list(
      beta = beta_init,
      u = u_init,
      log_sigma = log(sigma_init),
      log_sigma_u = log(sigma_u_init)
    )
  } else {
    tmb_params <- start_params
  }

  # Create TMB object
  tryCatch({
    obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = tmb_params,
      random = "u",  # Integrate out random effects
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
  par_full <- obj$env$last.par.best  # Includes random effects

  # Fixed effects
  beta_hat <- par_full[names(par_full) == "beta"]
  names(beta_hat) <- colnames(model_data$X)

  # Variance components
  log_sigma_hat <- par_full[names(par_full) == "log_sigma"]
  log_sigma_u_hat <- par_full[names(par_full) == "log_sigma_u"]
  sigma_hat <- exp(log_sigma_hat)
  sigma_u_hat <- exp(log_sigma_u_hat)

  # Random effects
  u_hat <- par_full[names(par_full) == "u"]

  # Organize random effects by group
  n_re_per_group <- tmb_data$n_random
  random_effects <- list()
  for (g in 1:tmb_data$n_groups) {
    idx_start <- (g - 1) * n_re_per_group + 1
    idx_end <- g * n_re_per_group
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
    g <- tmb_data$groups[i] + 1  # Convert back to 1-indexed
    fitted_vals[i] <- fitted_vals[i] + sum(model_data$Z[[1]][i, ] * random_effects[[g]])
  }

  # Log-likelihood
  loglik <- -opt$objective

  # Number of parameters
  n_params <- length(beta_hat) + 2  # beta + sigma + sigma_u

  list(
    coefficients = list(
      fixed = beta_hat,
      random_var = list(sigma_u_hat^2)
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
