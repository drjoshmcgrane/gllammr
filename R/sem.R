#' Fit Structural Equation Models with Latent Variables
#'
#' Fits a SEM with continuous indicators: a measurement model defining each
#' latent variable by its indicators (marker-variable identification: the
#' first indicator's loading is fixed at 1) and an optional recursive
#' structural model of regressions among latent variables.
#'
#' @param measurement Named list of one-sided formulas defining each latent
#'   variable, e.g. \code{list(ability = ~ y1 + y2 + y3,
#'   motivation = ~ y4 + y5 + y6)}.
#' @param structural Optional list of formulas among latent variables, e.g.
#'   \code{list(motivation ~ ability)}. Must be recursive (no cycles).
#' @param data Data frame with the indicator variables
#' @param start Optional starting values
#' @param control Optimization control list
#'
#' @return An object of class \code{gllamm_sem}
#'
#' @examples
#' \dontrun{
#' fit <- fit_sem(
#'   measurement = list(f1 = ~ x1 + x2 + x3, f2 = ~ y1 + y2 + y3),
#'   structural = list(f2 ~ f1),
#'   data = d)
#' }
#'
#' @export
fit_sem <- function(measurement, structural = NULL, data,
                    start = NULL, control = list()) {

  if (!is.list(measurement) || is.null(names(measurement)) ||
      any(names(measurement) == "")) {
    stop("measurement must be a named list of one-sided formulas")
  }
  latent_names <- names(measurement)
  n_latent <- length(latent_names)

  # ---- Indicators and loading pattern ----
  indicators_per_factor <- lapply(measurement, function(f) {
    all.vars(f)
  })
  indicator_names <- unique(unlist(indicators_per_factor))
  n_indicators <- length(indicator_names)
  missing_vars <- setdiff(indicator_names, names(data))
  if (length(missing_vars) > 0) {
    stop("Indicator variables not in data: ", paste(missing_vars, collapse = ", "))
  }

  # 0 = zero, 1 = free, 2 = fixed at 1 (first indicator of each factor)
  lambda_pattern <- matrix(0L, n_indicators, n_latent,
                           dimnames = list(indicator_names, latent_names))
  for (j in seq_len(n_latent)) {
    ind <- indicators_per_factor[[j]]
    lambda_pattern[ind[1], j] <- 2L
    if (length(ind) > 1) {
      lambda_pattern[ind[-1], j] <- 1L
    }
  }

  # ---- Structural pattern (recursive) ----
  beta_pattern <- matrix(0L, n_latent, n_latent,
                         dimnames = list(latent_names, latent_names))
  if (!is.null(structural)) {
    if (inherits(structural, "formula")) structural <- list(structural)
    for (f in structural) {
      outcome <- as.character(f[[2]])
      preds <- all.vars(f[[3]])
      if (!outcome %in% latent_names || !all(preds %in% latent_names)) {
        stop("Structural formulas must relate latent variables: ",
             paste(latent_names, collapse = ", "))
      }
      beta_pattern[outcome, preds] <- 1L
    }
    # Recursivity check: the directed graph must be acyclic
    bp <- beta_pattern; reach <- bp
    for (k in seq_len(n_latent)) reach <- pmin(reach + reach %*% bp, 1L)
    if (any(diag(reach) > 0)) {
      stop("The structural model contains a cycle; only recursive models ",
           "are supported")
    }
  }

  Y <- as.matrix(data[, indicator_names, drop = FALSE])
  if (anyNA(Y)) {
    stop("Indicator variables must be complete (no NA values)")
  }
  n_obs <- nrow(Y)

  tmb_data <- list(
    Y = Y,
    n_obs = as.integer(n_obs),
    n_indicators = as.integer(n_indicators),
    n_latent = as.integer(n_latent),
    lambda_pattern = lambda_pattern,
    beta_pattern = beta_pattern,
    model_name = "sem"
  )

  n_lambda_free <- sum(lambda_pattern == 1L)
  n_beta_free <- sum(beta_pattern == 1L)

  if (is.null(start)) {
    tmb_params <- list(
      nu = colMeans(Y),
      lambda_free = rep(1, n_lambda_free),
      beta_free = rep(0, max(n_beta_free, 1L)),
      log_psi = rep(log(stats::sd(Y[, 1])), n_latent),
      log_theta = log(apply(Y, 2, stats::sd) / 2),
      eta = matrix(0, n_obs, n_latent)
    )
  } else {
    tmb_params <- start
  }

  tmb_map <- list()
  if (n_beta_free == 0) {
    tmb_map$beta_free <- factor(rep(NA, length(tmb_params$beta_free)))
  }

  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = "eta",
    map = tmb_map,
    DLL = "GLLAMMR",
    silent = TRUE
  )

  control_defaults <- list(eval.max = 3000, iter.max = 2000, trace = 0)
  control <- modifyList(control_defaults, control)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)

  sdr <- try(TMB::sdreport(obj), silent = TRUE)
  par_full <- obj$env$last.par.best

  # ---- Reassemble Lambda and Beta ----
  Lambda <- matrix(0, n_indicators, n_latent,
                   dimnames = dimnames(lambda_pattern))
  Lambda[lambda_pattern == 2L] <- 1
  Lambda[lambda_pattern == 1L] <- par_full[names(par_full) == "lambda_free"]

  Beta <- matrix(0, n_latent, n_latent, dimnames = dimnames(beta_pattern))
  if (n_beta_free > 0) {
    Beta[beta_pattern == 1L] <- par_full[names(par_full) == "beta_free"]
  }

  psi_sd <- exp(unname(par_full[names(par_full) == "log_psi"]))
  theta_sd <- exp(unname(par_full[names(par_full) == "log_theta"]))
  names(psi_sd) <- latent_names
  names(theta_sd) <- indicator_names

  eta_hat <- matrix(par_full[names(par_full) == "eta"], n_obs, n_latent,
                    dimnames = list(NULL, latent_names))

  n_params <- n_indicators + n_lambda_free + n_beta_free +
    n_latent + n_indicators

  result <- list(
    loadings = Lambda,
    structural = Beta,
    latent_residual_sd = psi_sd,
    indicator_residual_sd = theta_sd,
    intercepts = setNames(unname(par_full[names(par_full) == "nu"]),
                          indicator_names),
    factor_scores = eta_hat,
    logLik = -opt$objective,
    AIC = 2 * opt$objective + 2 * n_params,
    BIC = 2 * opt$objective + log(n_obs) * n_params,
    convergence = list(converged = (opt$convergence == 0),
                       message = opt$message),
    n_obs = n_obs,
    measurement = measurement,
    structural_formulas = structural,
    tmb_obj = obj,
    tmb_opt = opt,
    tmb_sdr = sdr
  )
  class(result) <- c("gllamm_sem", "gllamm")
  result
}


#' @export
print.gllamm_sem <- function(x, ...) {
  cat("Structural Equation Model (continuous indicators)\n\n")
  cat("Observations:", x$n_obs, "\n\n")
  cat("Loadings (marker indicators fixed at 1):\n")
  print(round(x$loadings, 4))
  if (any(x$structural != 0)) {
    cat("\nStructural paths (row ~ column):\n")
    print(round(x$structural, 4))
  }
  cat("\nLatent residual SDs:\n")
  print(round(x$latent_residual_sd, 4))
  cat("\nLog-likelihood:", round(x$logLik, 2), "\n")
  invisible(x)
}