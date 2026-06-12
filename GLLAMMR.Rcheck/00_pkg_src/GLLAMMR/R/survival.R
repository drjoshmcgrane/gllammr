#' Fit Parametric Survival Models with Random Effects (Frailty)
#'
#' Fits exponential or Weibull proportional-hazards models with log-normal
#' frailty (normally distributed random effects on the log-hazard scale),
#' supporting right censoring.
#'
#' @param formula Formula of the form \code{Surv(time, event) ~ x + (1 | group)}.
#'   The left-hand side names the time and event (1 = event, 0 = censored)
#'   variables; the \pkg{survival} package is not required.
#' @param data Data frame
#' @param distribution "weibull" (default) or "exponential"
#' @param weights Optional case weights
#' @param start Optional starting values
#' @param control Optimization control list
#'
#' @details
#' The cumulative hazard is \eqn{H(t \mid x, u) = (\lambda t)^{shape}} with
#' \eqn{\lambda = \exp(x'\beta + z'u)}; the exponential model fixes
#' shape = 1. Under this parameterization the accelerated-failure-time
#' coefficients of \code{survival::survreg} correspond to \eqn{-\beta} and
#' its scale to \eqn{1/shape}.
#'
#' The exponential frailty model is likelihood-equivalent to a Poisson GLMM
#' on the event indicator with offset \eqn{\log(t)}, which is used in the
#' package validation suite.
#'
#' @return An object of class \code{gllamm_survival}
#'
#' @examples
#' \dontrun{
#' fit <- fit_survival(Surv(time, status) ~ age + (1 | center),
#'                     data = d, distribution = "weibull")
#' }
#'
#' @export
fit_survival <- function(formula, data,
                         distribution = c("weibull", "exponential"),
                         weights = NULL,
                         start = NULL, control = list()) {

  distribution <- match.arg(distribution)
  dist_code <- switch(distribution, exponential = 1L, weibull = 2L)

  # ---- Parse Surv(time, event) on the LHS without requiring survival ----
  lhs <- formula[[2]]
  if (!(is.call(lhs) && identical(as.character(lhs[[1]]), "Surv"))) {
    stop("The left-hand side must be Surv(time, event)")
  }
  time_var <- as.character(lhs[[2]])
  event_var <- as.character(lhs[[3]])
  if (!all(c(time_var, event_var) %in% names(data))) {
    stop("Variables '", time_var, "' and '", event_var,
         "' must be present in data")
  }
  time <- as.numeric(data[[time_var]])
  event <- as.integer(data[[event_var]])
  if (any(time <= 0, na.rm = TRUE)) {
    stop("Survival times must be positive")
  }
  if (!all(event %in% c(0L, 1L))) {
    stop("Event indicator must be 0 (censored) or 1 (event)")
  }

  # Rebuild a formula with the time variable as response so the standard
  # design-matrix machinery applies
  rhs <- deparse(formula[[3]], width.cutoff = 500)
  design_formula <- as.formula(paste(time_var, "~", paste(rhs, collapse = " ")))

  parsed <- parse_formula(design_formula, data)
  if (length(parsed$random_terms) == 0) {
    stop("No random effects specified; use survival::survreg for fixed-effects",
         " parametric survival models")
  }
  model_data <- make_model_matrices(parsed, data)

  if (model_data$n_random_terms != 1) {
    stop("Currently only a single random-effects term is supported for ",
         "survival models")
  }

  n_random <- model_data$n_random_coefs[1]
  correlated <- !parsed$random_terms[[1]]$uncorrelated
  Z_sparse <- Matrix::Matrix(model_data$Z[[1]], sparse = TRUE)

  if (is.null(weights)) {
    weights_vec <- rep(1.0, model_data$n_obs)
  } else {
    if (length(weights) != model_data$n_obs) {
      stop("weights length must match the number of observations")
    }
    weights_vec <- as.numeric(weights)
  }

  tmb_data <- list(
    time = time,
    event = event,
    X = as.matrix(model_data$X),
    Z = Z_sparse,
    groups = as.integer(model_data$groups[[1]]),
    n_groups = as.integer(model_data$n_groups[1]),
    n_obs = as.integer(model_data$n_obs),
    n_fixed = as.integer(model_data$n_fixed),
    n_random = as.integer(n_random),
    distribution = dist_code,
    correlated = as.integer(correlated),
    weights = weights_vec,
    model_name = "survival"
  )

  if (is.null(start)) {
    n_theta <- max(1L, n_random * (n_random - 1) %/% 2)
    tmb_params <- list(
      beta = rep(0, model_data$n_fixed),
      u = rep(0, tmb_data$n_groups * n_random),
      log_shape = 0,
      log_sigma_u = rep(log(0.5), n_random),
      theta = rep(0, n_theta)
    )
  } else {
    tmb_params <- start
  }

  # Map parameters the model never reads
  tmb_map <- list()
  if (dist_code == 1L) {
    tmb_map$log_shape <- factor(NA)   # exponential: shape fixed at 1
  }
  if (!(correlated && n_random > 1)) {
    tmb_map$theta <- factor(rep(NA, length(tmb_params$theta)))
  }

  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = "u",
    map = tmb_map,
    DLL = "GLLAMMR",
    silent = TRUE
  )

  control_defaults <- list(eval.max = 2000, iter.max = 1000, trace = 0)
  control <- modifyList(control_defaults, control)

  opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)

  sdr <- try(TMB::sdreport(obj), silent = TRUE)
  par_full <- obj$env$last.par.best

  beta_hat <- par_full[names(par_full) == "beta"]
  names(beta_hat) <- colnames(model_data$X)

  shape_hat <- if (dist_code == 2L) {
    exp(unname(par_full[names(par_full) == "log_shape"]))
  } else 1

  log_sigma_u_hat <- par_full[names(par_full) == "log_sigma_u"]
  sigma_u_hat <- exp(unname(log_sigma_u_hat))

  u_hat <- par_full[names(par_full) == "u"]
  random_effects <- matrix(u_hat, ncol = n_random, byrow = TRUE,
                           dimnames = list(NULL, colnames(model_data$Z[[1]])))

  n_params <- length(beta_hat) + n_random + (dist_code == 2L) +
    ifelse(correlated && n_random > 1, n_random * (n_random - 1) / 2, 0)

  result <- list(
    coefficients = list(
      fixed = beta_hat,
      random_sd = sigma_u_hat,
      random_var = list(diag(sigma_u_hat^2, nrow = n_random))
    ),
    shape = shape_hat,
    distribution = distribution,
    random_effects = random_effects,
    logLik = -opt$objective,
    AIC = 2 * opt$objective + 2 * n_params,
    BIC = 2 * opt$objective + log(model_data$n_obs) * n_params,
    convergence = list(
      converged = (opt$convergence == 0),
      message = opt$message
    ),
    n_obs = model_data$n_obs,
    n_groups = model_data$n_groups[1],
    n_events = sum(event),
    formula = formula,
    data = data,
    X = as.matrix(model_data$X),
    time_var = time_var,
    event_var = event_var,
    tmb_obj = obj,
    tmb_opt = opt,
    tmb_sdr = sdr
  )

  class(result) <- c("gllamm_survival", "gllamm")
  result
}


#' @export
print.gllamm_survival <- function(x, ...) {
  cat("Parametric Survival Model with Frailty (",
      x$distribution, ")\n\n", sep = "")
  cat("Observations:", x$n_obs, " Events:", x$n_events,
      " Groups:", x$n_groups, "\n\n")

  cat("Fixed effects (log relative cumulative hazard scale):\n")
  print(round(x$coefficients$fixed, 4))

  if (x$distribution == "weibull") {
    cat("\nWeibull shape:", round(x$shape, 4), "\n")
  }
  cat("\nFrailty SD:", round(x$coefficients$random_sd, 4), "\n")
  cat("Log-likelihood:", round(x$logLik, 2), "\n")
  cat("AIC:", round(x$AIC, 2), "  BIC:", round(x$BIC, 2), "\n")
  invisible(x)
}