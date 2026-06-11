#' Fit Two-Level GLMMs by Nonparametric Maximum Likelihood (NPML)
#'
#' Replaces the normal random-intercept distribution with a discrete
#' distribution on \code{k} estimated mass points (locations and masses).
#' The marginal likelihood is an exact finite mixture - no integral
#' approximation - making NPML robust to misspecification of the
#' random-effects distribution.
#'
#' @param formula Model formula \code{y ~ x + (1 | group)}
#' @param data Data frame
#' @param k Number of mass points (default 2)
#' @param family gaussian(), binomial(), or poisson()
#' @param weights Optional observation weights
#' @param n_starts Random restarts (mixtures are prone to local optima)
#' @param start Optional starting values
#' @param control Optimization control list
#'
#' @details
#' The fixed-effects design drops its intercept: the mass-point locations
#' play the role of component intercepts (the usual NPML identification,
#' as in \pkg{npmlreg}). Masses are softmax-parameterized.
#'
#' @return An object of class \code{gllamm_npml}
#'
#' @examples
#' \dontrun{
#' fit <- fit_npml(y ~ x + (1 | group), data = d, k = 3,
#'                 family = stats::binomial())
#' }
#'
#' @export
fit_npml <- function(formula, data, k = 2, family = stats::gaussian(),
                     weights = NULL, n_starts = 3,
                     start = NULL, control = list()) {

  if (is.function(family)) family <- family()
  fam_code <- switch(family$family,
                     gaussian = 0L, binomial = 1L, poisson = 2L,
                     stop("NPML supports gaussian, binomial, and poisson"))

  parsed <- parse_formula(formula, data)
  if (length(parsed$random_terms) != 1) {
    stop("NPML requires exactly one random-intercept term (1 | group)")
  }
  model_data <- make_model_matrices(parsed, data)
  if (model_data$n_random_coefs[1] != 1) {
    stop("NPML currently supports random intercepts only")
  }

  # Locations replace the intercept
  X <- drop_intercept_column(as.matrix(model_data$X))

  weights_vec <- if (is.null(weights)) rep(1.0, model_data$n_obs)
                 else as.numeric(align_weights(weights, model_data))

  tmb_data <- list(
    y = as.numeric(model_data$y),
    X = X,
    groups = as.integer(model_data$groups[[1]]),
    n_obs = as.integer(model_data$n_obs),
    n_groups = as.integer(model_data$n_groups[1]),
    K = as.integer(k),
    family = fam_code,
    weights = weights_vec,
    model_name = "npml"
  )

  # Spread of the response on the linear-predictor scale for initial locations
  y_center <- if (fam_code == 0L) mean(model_data$y)
              else if (fam_code == 1L) qlogis(pmin(pmax(mean(model_data$y), 0.05), 0.95))
              else log(max(mean(model_data$y), 0.1))
  y_spread <- if (fam_code == 0L) stats::sd(model_data$y) else 1

  make_params <- function(jitter = 0) {
    locs <- y_center + y_spread * seq(-1, 1, length.out = k) +
      rnorm(k, 0, jitter * y_spread)
    list(beta = rep(0, ncol(X)),
         locations = sort(locs),
         mass_logits = rnorm(max(k - 1, 1), 0, jitter),
         log_sigma = if (fam_code == 0L) log(max(y_spread / 2, 0.1)) else 0)
  }

  tmb_map <- list()
  if (fam_code != 0L) tmb_map$log_sigma <- factor(NA)
  if (k == 1) tmb_map$mass_logits <- factor(NA)

  best <- NULL
  for (s in seq_len(n_starts)) {
    tmb_params <- if (!is.null(start)) start else make_params(jitter = (s - 1) * 0.4)
    obj <- TMB::MakeADFun(data = tmb_data, parameters = tmb_params,
                          map = tmb_map, DLL = "GLLAMMR", silent = TRUE)
    control_defaults <- list(eval.max = 3000, iter.max = 1500, trace = 0)
    ctl <- modifyList(control_defaults, control)
    opt <- try(nlminb(obj$par, obj$fn, obj$gr, control = ctl), silent = TRUE)
    if (!inherits(opt, "try-error") &&
        (is.null(best) || opt$objective < best$opt$objective)) {
      best <- list(obj = obj, opt = opt)
    }
    if (!is.null(start)) break
  }
  if (is.null(best)) stop("All NPML optimization attempts failed")
  obj <- best$obj; opt <- best$opt

  sdr <- try(TMB::sdreport(obj), silent = TRUE)
  par_full <- obj$env$last.par.best

  beta_hat <- par_full[names(par_full) == "beta"]
  names(beta_hat) <- colnames(X)
  locations <- unname(par_full[names(par_full) == "locations"])
  ml <- unname(par_full[names(par_full) == "mass_logits"])
  masses <- if (k == 1) 1 else {
    e <- exp(c(ml, 0)); e / sum(e)
  }
  ord <- order(locations)
  locations <- locations[ord]; masses <- masses[ord]

  # Moments of the discrete random-effects distribution
  re_mean <- sum(masses * locations)
  re_sd <- sqrt(sum(masses * (locations - re_mean)^2))

  n_params <- length(beta_hat) + k + (k - 1) + (fam_code == 0L)

  result <- list(
    coefficients = list(fixed = beta_hat),
    locations = locations,
    masses = masses,
    re_mean = re_mean,
    re_sd = re_sd,
    residual_sd = if (fam_code == 0L) {
      exp(unname(par_full[names(par_full) == "log_sigma"]))
    } else NA_real_,
    k = k,
    family = family,
    logLik = -opt$objective,
    AIC = 2 * opt$objective + 2 * n_params,
    BIC = 2 * opt$objective + log(model_data$n_obs) * n_params,
    convergence = list(converged = (opt$convergence == 0),
                       message = opt$message),
    n_obs = model_data$n_obs,
    n_groups = model_data$n_groups[1],
    formula = formula,
    data = data,
    tmb_obj = obj,
    tmb_opt = opt,
    tmb_sdr = sdr
  )
  class(result) <- c("gllamm_npml", "gllamm")
  result
}


#' @export
print.gllamm_npml <- function(x, ...) {
  cat("NPML Two-Level GLMM (", x$family$family, ", ", x$k,
      " mass points)\n\n", sep = "")
  cat("Observations:", x$n_obs, " Groups:", x$n_groups, "\n\n")
  cat("Fixed effects:\n")
  print(round(x$coefficients$fixed, 4))
  cat("\nMass points:\n")
  print(round(rbind(location = x$locations, mass = x$masses), 4))
  cat("\nImplied RE mean:", round(x$re_mean, 4),
      " SD:", round(x$re_sd, 4), "\n")
  if (!is.na(x$residual_sd)) {
    cat("Residual SD:", round(x$residual_sd, 4), "\n")
  }
  cat("Log-likelihood:", round(x$logLik, 2), "\n")
  invisible(x)
}
