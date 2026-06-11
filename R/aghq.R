#' Adaptive quadrature integration specification
#'
#' Use as \code{gllamm(..., integration = aghq(k))} to integrate the random
#' intercept by adaptive Gauss-Hermite quadrature with \code{k} nodes
#' instead of the Laplace approximation. Currently supports two-level
#' random-intercept models with gaussian, binomial, or poisson families.
#' Laplace (the default) is equivalent to \code{aghq(1)}.
#'
#' @param k Number of quadrature nodes (default 15)
#' @return An object of class \code{gllamm_integration}
#' @export
aghq <- function(k = 15) {
  if (!is.numeric(k) || length(k) != 1 || k < 2) {
    stop("k must be a single integer >= 2")
  }
  structure(list(method = "aghq", k = as.integer(k)),
            class = "gllamm_integration")
}


#' Nonparametric maximum likelihood integration
#'
#' Integration specification for \code{\link{gllamm}}: replace the normal
#' latent distribution with k estimated mass points and masses (NPML;
#' Aitkin 1999). Equivalent to \code{\link{fit_npml}}.
#'
#' @param k Number of mass points (default 2)
#' @return An object of class \code{gllamm_integration}
#' @examples
#' \dontrun{
#' gllamm(y ~ x + (1 | g), data = d, family = binomial(),
#'        integration = npml(2))
#' }
#' @export
npml <- function(k = 2) {
  if (!is.numeric(k) || length(k) != 1 || k < 2) {
    stop("k must be a single integer >= 2")
  }
  structure(list(method = "npml", k = as.integer(k)),
            class = "gllamm_integration")
}


#' Gauss-Hermite nodes and weights (Golub-Welsch)
#' @keywords internal
gauss_hermite <- function(n) {
  i <- seq_len(n - 1)
  J <- matrix(0, n, n)
  off <- sqrt(i / 2)
  J[cbind(i, i + 1)] <- off
  J[cbind(i + 1, i)] <- off
  e <- eigen(J, symmetric = TRUE)
  x <- e$values
  w <- sqrt(pi) * e$vectors[1, ]^2
  ord <- order(x)
  list(nodes = x[ord], weights = w[ord])
}


#' Fit a two-level GLMM by adaptive Gauss-Hermite quadrature
#'
#' R driver for the glmm_aghq TMB objective: alternates parameter
#' optimization with updates of the per-group adaptation centers and scales
#' (posterior modes and curvatures), the classic adaptive quadrature scheme.
#'
#' @keywords internal
fit_tmb_gllamm_aghq <- function(model_data, family, random_terms,
                                k = 15, start_params = NULL,
                                control = list(), weights = NULL,
                                max_adapt = 5, adapt_tol = 1e-4) {
  if (model_data$n_random_terms != 1 || model_data$n_random_coefs[1] != 1) {
    stop("Adaptive quadrature currently supports a single random intercept")
  }
  fam_code <- switch(family$family,
                     gaussian = 0L, binomial = 1L, poisson = 2L,
                     stop("aghq supports gaussian, binomial, poisson"))
  link_code <- if (family$family == "binomial") {
    switch(family$link, logit = 1L, probit = 2L, cloglog = 3L, 1L)
  } else 1L

  n_obs <- model_data$n_obs
  n_groups <- model_data$n_groups[1]
  groups <- model_data$groups[[1]]

  w <- parse_level_weights(weights, n_obs, groups, n_groups)
  gh <- gauss_hermite(k)

  tmb_data <- list(
    y = as.numeric(model_data$y),
    X = as.matrix(model_data$X),
    groups = as.integer(groups),
    n_obs = as.integer(n_obs),
    n_groups = as.integer(n_groups),
    family = fam_code,
    link = link_code,
    weights = w$level1,
    group_weights = w$level2,
    gh_x = gh$nodes,
    gh_logw = log(gh$weights) + gh$nodes^2,
    center = rep(0, n_groups),
    scale = rep(1, n_groups),
    model_name = "glmm_aghq"
  )

  if (is.null(start_params)) {
    if (family$family == "gaussian") {
      lm_fit <- lm(model_data$y ~ model_data$X - 1)
      beta_init <- coef(lm_fit)
      sigma_init <- summary(lm_fit)$sigma
    } else {
      glm_fit <- glm(model_data$y ~ model_data$X - 1, family = family)
      beta_init <- coef(glm_fit)
      sigma_init <- 1.0
    }
    tmb_params <- list(beta = beta_init,
                       log_sigma = log(max(sigma_init, 0.1)),
                       log_sigma_u = log(0.5))
  } else {
    tmb_params <- start_params
  }

  tmb_map <- list()
  if (family$family != "gaussian") tmb_map$log_sigma <- factor(NA)

  control_defaults <- list(eval.max = 2000, iter.max = 1000, trace = 0)
  control$optimizer <- NULL
  ctl <- modifyList(control_defaults, control)

  # Per-group conditional log-likelihood for the adaptation step
  xb_loglik <- function(beta, sigma) {
    xb <- as.numeric(tmb_data$X %*% beta)
    function(g, u) {
      idx <- which(groups == g)
      eta <- xb[idx] + u
      ll <- switch(as.character(fam_code),
        "0" = dnorm(tmb_data$y[idx], eta, sigma, log = TRUE),
        "1" = {
          p <- switch(as.character(link_code),
                      "2" = pnorm(eta),
                      "3" = 1 - exp(-exp(eta)),
                      plogis(eta))
          tmb_data$y[idx] * log(p + 1e-12) +
            (1 - tmb_data$y[idx]) * log(1 - p + 1e-12)
        },
        "2" = dpois(tmb_data$y[idx], exp(eta), log = TRUE))
      sum(w$level1[idx] * ll)
    }
  }

  opt <- NULL; obj <- NULL
  for (round in seq_len(max_adapt)) {
    obj <- TMB::MakeADFun(data = tmb_data, parameters = tmb_params,
                          map = tmb_map, DLL = "GLLAMMR", silent = TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = ctl)
    tmb_params$beta <- opt$par[names(opt$par) == "beta"]
    if (fam_code == 0L) {
      tmb_params$log_sigma <- unname(opt$par[names(opt$par) == "log_sigma"])
    }
    tmb_params$log_sigma_u <- unname(opt$par[names(opt$par) == "log_sigma_u"])

    # Update adaptation: posterior mode and curvature per group
    sigma_u <- exp(tmb_params$log_sigma_u)
    sigma_res <- if (fam_code == 0L) exp(tmb_params$log_sigma) else 1
    ll_fun <- xb_loglik(tmb_params$beta, sigma_res)

    new_center <- numeric(n_groups)
    new_scale <- numeric(n_groups)
    for (g in seq_len(n_groups) - 1L) {
      post <- function(u) {
        dnorm(u, 0, sigma_u, log = TRUE) + ll_fun(g, u)
      }
      o <- optimize(post, interval = c(-6, 6) * max(sigma_u, 0.5),
                    maximum = TRUE)
      m_g <- o$maximum
      h <- 1e-3 * max(abs(m_g), 1)
      curv <- -(post(m_g + h) - 2 * post(m_g) + post(m_g - h)) / h^2
      new_center[g + 1] <- m_g
      new_scale[g + 1] <- 1 / sqrt(max(curv, 1e-8)) / sqrt(2)
    }

    delta <- max(abs(new_center - tmb_data$center),
                 abs(new_scale - tmb_data$scale))
    tmb_data$center <- new_center
    tmb_data$scale <- new_scale
    if (delta < adapt_tol) break
  }

  # Final fit at converged adaptation points
  obj <- TMB::MakeADFun(data = tmb_data, parameters = tmb_params,
                        map = tmb_map, DLL = "GLLAMMR", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = ctl)
  sdr <- try(TMB::sdreport(obj), silent = TRUE)

  par_full <- obj$env$last.par.best
  beta_hat <- par_full[names(par_full) == "beta"]
  names(beta_hat) <- colnames(model_data$X)
  sigma_u_hat <- exp(unname(par_full[names(par_full) == "log_sigma_u"]))

  if (!inherits(sdr, "try-error") && !is.null(sdr)) {
    se_all <- summary(sdr, "fixed")
    se_beta <- se_all[rownames(se_all) == "beta", "Std. Error"]
    vcov_fixed <- diag(se_beta^2, length(beta_hat))
    dimnames(vcov_fixed) <- list(names(beta_hat), names(beta_hat))
  } else {
    vcov_fixed <- matrix(NA, length(beta_hat), length(beta_hat))
  }

  # Empirical Bayes modes serve as random-effect estimates
  random_effects <- as.list(tmb_data$center)

  fitted_vals <- as.numeric(model_data$X %*% beta_hat) +
    tmb_data$center[groups + 1]
  if (family$family != "gaussian") {
    fitted_vals <- family$linkinv(fitted_vals)
  }

  n_params <- length(beta_hat) + 1 + (fam_code == 0L)

  list(
    coefficients = list(fixed = beta_hat,
                        random_var = list(matrix(sigma_u_hat^2, 1, 1))),
    vcov = list(fixed = vcov_fixed, all = NULL),
    random_effects = random_effects,
    fitted.values = fitted_vals,
    logLik = -opt$objective,
    n_params = n_params,
    integration = list(method = "aghq", k = k),
    convergence = list(converged = (opt$convergence == 0),
                       message = opt$message,
                       iterations = opt$iterations),
    tmb_obj = obj,
    tmb_opt = opt,
    tmb_sdr = sdr
  )
}
