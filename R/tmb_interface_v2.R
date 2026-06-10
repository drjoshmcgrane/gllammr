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

  # Multiple random-effects terms (crossed or expanded nested structures)
  # route to the combined-Z engine
  if (model_data$n_random_terms > 1) {
    return(fit_tmb_gllamm_multi(model_data, family, random_terms,
                                start_params = start_params,
                                control = control, weights = weights))
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

  # Select model template: family + random-effects structure. Gamma always
  # routes through the general template (no gamma branch in the dedicated
  # single-term templates).
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
                            gaussian = 0L, binomial = 1L, poisson = 2L,
                            Gamma = 3L, 0L)
  tmb_data$link <- if (family$family == "binomial") {
    switch(family$link, logit = 1L, probit = 2L, cloglog = 3L, 1L)
  } else if (family$family == "Gamma") {
    switch(family$link, log = 1L, inverse = 2L, identity = 3L, 1L)
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
      sigma_init <- if (family$family == "Gamma") {
        max(summary(glm_fit)$dispersion, 0.05)
      } else 1.0
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
  if (model_name == "glmm_slopes" &&
      !(family$family %in% c("gaussian", "Gamma"))) {
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

  # Optimize. control$optimizer selects the outer optimizer:
  # "nlminb" (default) or "BFGS" (stats::optim).
  optimizer <- control$optimizer %||% "nlminb"
  control$optimizer <- NULL
  control_defaults <- list(
    eval.max = 1000,
    iter.max = 500,
    trace = 0
  )
  control <- modifyList(control_defaults, control)

  if (identical(optimizer, "BFGS")) {
    opt <- try(
      optim(obj$par, obj$fn, obj$gr, method = "BFGS",
            control = list(maxit = control$iter.max, trace = control$trace)),
      silent = FALSE
    )
    if (!inherits(opt, "try-error")) {
      # Align with the nlminb result shape used downstream
      opt$objective <- opt$value
      opt$iterations <- opt$counts[["function"]]
    }
  } else {
    opt <- try(
      nlminb(
        start = obj$par,
        objective = obj$fn,
        gradient = obj$gr,
        control = control
      ),
      silent = FALSE
    )
  }

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

  # Fitted values (vectorized: rows of Z times each observation's group RE)
  u_mat <- matrix(u_hat, nrow = tmb_data$n_groups, ncol = n_random, byrow = TRUE)
  fitted_vals <- as.numeric(model_data$X %*% beta_hat) +
    rowSums(model_data$Z[[1]] * u_mat[tmb_data$groups + 1, , drop = FALSE])

  # Apply inverse link for GLMs
  if (family$family != "gaussian") {
    fitted_vals <- family$linkinv(fitted_vals)
  }

  # Log-likelihood
  loglik <- -opt$objective

  # Number of parameters
  n_params <- length(beta_hat) + n_random +
    ifelse(family$family %in% c("gaussian", "Gamma"), 1, 0)
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


#' GLMM engine for multiple random-effects terms (crossed/nested)
#'
#' Builds the lme4-style combined sparse Z mapping the full term-major
#' random-effects vector to observations, and fits via the glmm_multi
#' TMB model.
#'
#' @keywords internal
fit_tmb_gllamm_multi <- function(model_data, family, random_terms,
                                 start_params = NULL, control = list(),
                                 weights = NULL) {
  n_terms <- model_data$n_random_terms
  n_obs <- model_data$n_obs

  term_n_random <- as.integer(model_data$n_random_coefs)
  term_n_groups <- as.integer(model_data$n_groups)
  term_correlated <- vapply(seq_len(n_terms), function(t) {
    nr <- term_n_random[t]
    as.integer(nr > 1 && !isTRUE(random_terms[[t]]$uncorrelated))
  }, integer(1))

  # ---- Combined sparse Z: column block per term, group-major within ----
  q_per_term <- term_n_groups * term_n_random
  offsets <- c(0L, cumsum(q_per_term))
  q_total <- offsets[n_terms + 1]

  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  for (t in seq_len(n_terms)) {
    Z_t <- model_data$Z[[t]]                 # n_obs x n_random_t
    g_t <- model_data$groups[[t]]            # 0-indexed group per obs
    nr <- term_n_random[t]
    for (k in seq_len(nr)) {
      ii <- c(ii, seq_len(n_obs))
      jj <- c(jj, offsets[t] + g_t * nr + k)
      xx <- c(xx, Z_t[, k])
    }
  }
  Z_combined <- Matrix::sparseMatrix(i = ii, j = jj, x = xx,
                                     dims = c(n_obs, q_total))

  if (is.null(weights)) {
    weights_vec <- rep(1.0, n_obs)
  } else {
    weights_vec <- as.numeric(weights)
  }

  tmb_data <- list(
    y = as.numeric(model_data$y),
    X = as.matrix(model_data$X),
    Z = Z_combined,
    n_obs = as.integer(n_obs),
    n_terms = as.integer(n_terms),
    term_n_random = term_n_random,
    term_n_groups = term_n_groups,
    term_correlated = term_correlated,
    family = switch(family$family,
                    gaussian = 0L, binomial = 1L, poisson = 2L, Gamma = 3L,
                    stop("Unsupported family for multi-term GLMM: ",
                         family$family)),
    link = if (family$family == "binomial") {
      switch(family$link, logit = 1L, probit = 2L, cloglog = 3L, 1L)
    } else if (family$family == "Gamma") {
      switch(family$link, log = 1L, inverse = 2L, identity = 3L, 1L)
    } else 1L,
    weights = weights_vec,
    model_name = "glmm_multi"
  )

  # ---- Parameters ----
  n_theta_per_term <- ifelse(term_correlated == 1L,
                             term_n_random * (term_n_random - 1) / 2, 0L)
  n_theta <- sum(n_theta_per_term)

  if (is.null(start_params)) {
    if (family$family == "gaussian") {
      lm_fit <- lm(model_data$y ~ model_data$X - 1)
      beta_init <- coef(lm_fit)
      sigma_init <- summary(lm_fit)$sigma
    } else {
      glm_fit <- glm(model_data$y ~ model_data$X - 1, family = family)
      beta_init <- coef(glm_fit)
      sigma_init <- if (family$family == "Gamma") {
        max(summary(glm_fit)$dispersion, 0.05)
      } else 1.0
    }
    tmb_params <- list(
      beta = beta_init,
      u = rep(0, q_total),
      log_sigma = log(max(sigma_init, 0.1)),
      log_sigma_u = rep(log(0.5), sum(term_n_random)),
      theta = rep(0, max(n_theta, 1L))
    )
  } else {
    tmb_params <- start_params
  }

  # Map off parameters the model never reads
  tmb_map <- list()
  if (n_theta == 0) {
    tmb_map$theta <- factor(rep(NA, length(tmb_params$theta)))
  }
  if (!(family$family %in% c("gaussian", "Gamma"))) {
    tmb_map$log_sigma <- factor(NA)
  }

  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = "u",
    map = tmb_map,
    DLL = "GLLAMMR",
    silent = TRUE
  )

  control_defaults <- list(eval.max = 1000, iter.max = 500, trace = 0)
  control$optimizer <- NULL
  control <- modifyList(control_defaults, control)

  opt <- try(
    nlminb(obj$par, obj$fn, obj$gr, control = control),
    silent = FALSE
  )
  if (inherits(opt, "try-error")) {
    stop("Optimization failed")
  }
  converged <- (opt$convergence == 0)

  sdr <- try(TMB::sdreport(obj), silent = TRUE)
  if (inherits(sdr, "try-error")) {
    warning("Failed to compute standard errors")
    sdr <- NULL
  }

  par_full <- obj$env$last.par.best

  beta_hat <- par_full[names(par_full) == "beta"]
  names(beta_hat) <- colnames(model_data$X)

  # ---- Per-term variance components ----
  log_sigma_u_hat <- par_full[names(par_full) == "log_sigma_u"]
  theta_hat <- par_full[names(par_full) == "theta"]
  sd_offsets <- c(0L, cumsum(term_n_random))
  th_offsets <- c(0L, cumsum(n_theta_per_term))

  Sigma_list <- vector("list", n_terms)
  for (t in seq_len(n_terms)) {
    nr <- term_n_random[t]
    sds <- exp(log_sigma_u_hat[(sd_offsets[t] + 1):sd_offsets[t + 1]])
    if (term_correlated[t] == 1L && nr > 1) {
      th <- theta_hat[(th_offsets[t] + 1):th_offsets[t + 1]]
      L <- diag(nr)
      idx <- 1
      for (i in 2:nr) {
        for (j in 1:(i - 1)) {
          L[i, j] <- th[idx]
          idx <- idx + 1
        }
      }
      R <- L %*% t(L)
      R <- R / sqrt(diag(R) %o% diag(R))
      Sigma_list[[t]] <- outer(sds, sds) * R
    } else {
      Sigma_list[[t]] <- diag(sds^2, nrow = nr)
    }
    dimnames(Sigma_list[[t]]) <- list(colnames(model_data$Z[[t]]),
                                      colnames(model_data$Z[[t]]))
  }
  names(Sigma_list) <- vapply(random_terms,
                              function(rt) paste(rt$grouping, collapse = ":"),
                              character(1))

  # ---- Random effects, organized per term as group x coef matrices ----
  u_hat <- par_full[names(par_full) == "u"]
  random_effects <- vector("list", n_terms)
  for (t in seq_len(n_terms)) {
    nr <- term_n_random[t]
    block <- u_hat[(offsets[t] + 1):offsets[t + 1]]
    random_effects[[t]] <- matrix(block, ncol = nr, byrow = TRUE,
                                  dimnames = list(NULL, colnames(model_data$Z[[t]])))
  }
  names(random_effects) <- names(Sigma_list)

  if (!is.null(sdr)) {
    se_all <- summary(sdr, "fixed")
    se_beta <- se_all[rownames(se_all) == "beta", "Std. Error"]
    vcov_fixed <- matrix(0, length(beta_hat), length(beta_hat))
    diag(vcov_fixed) <- se_beta^2
    dimnames(vcov_fixed) <- list(names(beta_hat), names(beta_hat))
  } else {
    vcov_fixed <- matrix(NA, length(beta_hat), length(beta_hat))
  }

  # Fitted values via the same combined-Z product as the template
  fitted_vals <- as.numeric(model_data$X %*% beta_hat) +
    as.numeric(Z_combined %*% u_hat)
  if (family$family != "gaussian") {
    fitted_vals <- family$linkinv(fitted_vals)
  }

  loglik <- -opt$objective
  n_params <- length(beta_hat) + sum(term_n_random) + n_theta +
    ifelse(family$family %in% c("gaussian", "Gamma"), 1, 0)

  list(
    coefficients = list(
      fixed = beta_hat,
      random_var = Sigma_list
    ),
    vcov = list(fixed = vcov_fixed, all = NULL),
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
