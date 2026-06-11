#' Fit Structural Equation Models with Latent Variables
#'
#' Fits a SEM with continuous indicators: a measurement model defining each
#' latent variable by its indicators (marker-variable identification: the
#' first indicator's loading is fixed at 1), an optional recursive
#' structural model of regressions among latent variables and on observed
#' covariates (MIMIC models), freely correlated exogenous latent
#' variables, full-information maximum likelihood for missing data, and
#' the standard fit-index battery.
#'
#' @param measurement Named list of one-sided formulas defining each latent
#'   variable, e.g. \code{list(ability = ~ y1 + y2 + y3,
#'   motivation = ~ y4 + y5 + y6)}.
#' @param structural Optional list of formulas regressing latent variables
#'   on other latent variables and/or observed covariates, e.g.
#'   \code{list(motivation ~ ability + ses)}. The latent part must be
#'   recursive (no cycles); observed covariates are carried as
#'   perfectly-measured exogenous variables (the joint-normal formulation,
#'   equivalent to \code{lavaan} with \code{fixed.x = FALSE}).
#' @param data Data frame with the indicator (and covariate) variables
#' @param method Estimation method: "ml" (default; Wishart maximum
#'   likelihood on the sample covariance for complete data, casewise FIML
#'   under \code{missing = "fiml"}) or "laplace" (legacy full-data TMB
#'   path; latent-only structural models, uncorrelated exogenous factors,
#'   complete data).
#' @param missing "listwise" (default) or "fiml" (full-information ML
#'   over all observed values, assuming MAR)
#' @param se Compute standard errors (default TRUE)
#' @param start Optional starting values (laplace method only)
#' @param control Optimization control list
#'
#' @return An object of class \code{gllamm_sem}. Key components:
#'   \code{param_table} (estimates, SEs, z, p), \code{fit_measures}
#'   (chisq, df, CFI, TLI, RMSEA with 90\% CI, SRMR),
#'   \code{latent_covariance} (exogenous block free; disturbances
#'   diagonal), \code{standardized} (std.all solution), \code{loadings},
#'   \code{structural}, \code{factor_scores}, \code{logLik}/\code{AIC}/
#'   \code{BIC}.
#'
#' @details
#' Exogenous latent variables (those with no incoming structural paths)
#' covary freely, as in lavaan. Likelihood-equivalent to lavaan with
#' \code{fixed.x = FALSE} when covariates are present (covariate moments
#' are modeled jointly), and to default lavaan otherwise.
#'
#' @examples
#' \dontrun{
#' fit <- fit_sem(
#'   measurement = list(f1 = ~ x1 + x2 + x3, f2 = ~ y1 + y2 + y3),
#'   structural = list(f2 ~ f1 + w),
#'   data = d, missing = "fiml")
#' summary(fit)
#' }
#'
#' @export
fit_sem <- function(measurement, structural = NULL, data,
                    method = c("ml", "laplace"),
                    missing = c("listwise", "fiml"),
                    se = TRUE,
                    start = NULL, control = list()) {

  method <- match.arg(method)
  missing <- match.arg(missing)

  if (!is.list(measurement) || is.null(names(measurement)) ||
      any(names(measurement) == "")) {
    stop("measurement must be a named list of one-sided formulas")
  }
  latent_names <- names(measurement)
  n_latent <- length(latent_names)

  # ---- Indicators and loading pattern ----
  indicators_per_factor <- lapply(measurement, function(f) all.vars(f))
  indicator_names <- unique(unlist(indicators_per_factor))
  n_indicators <- length(indicator_names)
  missing_vars <- setdiff(indicator_names, names(data))
  if (length(missing_vars) > 0) {
    stop("Indicator variables not in data: ",
         paste(missing_vars, collapse = ", "))
  }

  # ---- Structural pattern: latent and observed predictors ----
  covariate_names <- character(0)
  struct_pairs <- list()
  if (!is.null(structural)) {
    if (inherits(structural, "formula")) structural <- list(structural)
    for (f in structural) {
      outcome <- as.character(f[[2]])
      preds <- all.vars(f[[3]])
      if (!outcome %in% latent_names) {
        stop("Structural outcomes must be latent variables: ",
             paste(latent_names, collapse = ", "))
      }
      for (pr in preds) {
        if (!pr %in% latent_names) {
          if (!pr %in% names(data)) {
            stop("Structural predictor '", pr, "' is neither a latent ",
                 "variable nor a column of data")
          }
          covariate_names <- union(covariate_names, pr)
        }
        struct_pairs[[length(struct_pairs) + 1]] <- c(outcome, pr)
      }
    }
  }
  n_cov <- length(covariate_names)

  if (method == "laplace" && (n_cov > 0 || missing == "fiml")) {
    stop("method = \"laplace\" supports latent-only structural models ",
         "and complete data; use method = \"ml\"")
  }

  # Extended variable space: true latents + covariate pseudo-latents
  all_latent <- c(latent_names, covariate_names)
  q_all <- length(all_latent)
  all_obs <- c(indicator_names, covariate_names)
  p_all <- length(all_obs)

  # 0 = zero, 1 = free, 2 = fixed at 1
  lambda_pattern <- matrix(0L, p_all, q_all,
                           dimnames = list(all_obs, all_latent))
  for (j in seq_len(n_latent)) {
    ind <- indicators_per_factor[[j]]
    lambda_pattern[ind[1], j] <- 2L
    if (length(ind) > 1) lambda_pattern[ind[-1], j] <- 1L
  }
  for (w in covariate_names) lambda_pattern[w, w] <- 2L
  theta_zero <- all_obs %in% covariate_names

  beta_pattern <- matrix(0L, q_all, q_all,
                         dimnames = list(all_latent, all_latent))
  for (pp in struct_pairs) beta_pattern[pp[1], pp[2]] <- 1L

  # Recursivity check on the latent block
  bp <- beta_pattern
  reach <- bp
  for (k in seq_len(q_all)) reach <- pmin(reach + reach %*% bp, 1L)
  if (any(diag(reach) > 0)) {
    stop("The structural model contains a cycle; only recursive models ",
         "are supported")
  }

  Y <- as.matrix(data[, all_obs, drop = FALSE])
  if (!is.numeric(Y)) stop("Indicators and covariates must be numeric")
  if (anyNA(Y) && missing == "listwise" && method == "ml") {
    n_drop <- sum(!stats::complete.cases(Y))
    warning("Removing ", n_drop, " incomplete rows (listwise); use ",
            "missing = \"fiml\" for full-information ML")
  }
  if (anyNA(Y) && method == "laplace") {
    stop("method = \"laplace\" requires complete data")
  }

  if (method == "ml") {
    ml <- fit_sem_ml(Y, lambda_pattern, beta_pattern,
                     theta_zero = theta_zero, missing = missing,
                     se = se, control = control)
    n_obs <- ml$n_obs

    Lambda <- ml$Lambda[seq_len(n_indicators), seq_len(n_latent),
                        drop = FALSE]
    dimnames(Lambda) <- list(indicator_names, latent_names)
    Beta <- ml$B
    dimnames(Beta) <- list(all_latent, all_latent)
    Psi <- ml$Psi
    dimnames(Psi) <- list(all_latent, all_latent)
    V_eta <- ml$V_eta
    dimnames(V_eta) <- list(all_latent, all_latent)

    psi_sd <- setNames(sqrt(diag(Psi))[seq_len(n_latent)], latent_names)
    theta_sd <- setNames(sqrt(ml$theta_var[seq_len(n_indicators)]),
                         indicator_names)
    scores <- ml$factor_scores[, seq_len(n_latent), drop = FALSE]
    colnames(scores) <- latent_names

    # ---- Standardized (std.all) solution ----
    sd_eta <- sqrt(diag(V_eta))
    sd_y <- sqrt(diag(ml$Sigma))
    names(sd_y) <- all_obs
    std <- ml$param_table
    for (r in seq_len(nrow(std))) {
      lab <- std$label[r]
      if (grepl("=~", lab, fixed = TRUE)) {
        parts <- strsplit(lab, "=~", fixed = TRUE)[[1]]
        std$est[r] <- std$est[r] * sd_eta[parts[1]] / sd_y[parts[2]]
      } else if (grepl("~~", lab, fixed = TRUE)) {
        parts <- strsplit(lab, "~~", fixed = TRUE)[[1]]
        s1 <- if (parts[1] %in% all_latent) sd_eta[parts[1]]
              else sd_y[parts[1]]
        s2 <- if (parts[2] %in% all_latent) sd_eta[parts[2]]
              else sd_y[parts[2]]
        std$est[r] <- std$est[r] / (s1 * s2)
      } else if (grepl("~", lab, fixed = TRUE)) {
        parts <- strsplit(lab, "~", fixed = TRUE)[[1]]
        std$est[r] <- std$est[r] * sd_eta[parts[2]] / sd_eta[parts[1]]
      }
    }
    standardized <- data.frame(label = std$label, est_std = std$est,
                               stringsAsFactors = FALSE)

    result <- list(
      method = "ML",
      missing = ml$missing,
      loadings = Lambda,
      structural = Beta,
      latent_covariance = Psi,
      latent_variance_total = V_eta,
      latent_residual_sd = psi_sd,
      indicator_residual_sd = theta_sd,
      intercepts = setNames(ml$intercepts[seq_len(n_indicators)],
                            indicator_names),
      covariates = covariate_names,
      factor_scores = scores,
      param_table = ml$param_table,
      standardized = standardized,
      vcov = ml$vcov,
      fit_measures = ml$fit_measures,
      logLik = ml$logLik,
      logLik_saturated = ml$logLik_saturated,
      AIC = -2 * ml$logLik + 2 * ml$n_params,
      BIC = -2 * ml$logLik + log(n_obs) * ml$n_params,
      n_params = ml$n_params,
      convergence = list(converged = ml$converged, message = ml$message),
      n_obs = n_obs,
      measurement = measurement,
      structural_formulas = structural
    )
    class(result) <- c("gllamm_sem", "gllamm")
    return(result)
  }

  # ---- Legacy Laplace path (latent-only, orthogonal exogenous) ----
  exo_count <- sum(rowSums(beta_pattern) == 0)
  if (exo_count > 1) {
    warning("method = \"laplace\" treats exogenous latent variables as ",
            "uncorrelated; use method = \"ml\" (the default) for freely ",
            "correlated factors")
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

  Lambda <- matrix(0, n_indicators, n_latent,
                   dimnames = list(indicator_names, latent_names))
  Lambda[lambda_pattern == 2L] <- 1
  Lambda[lambda_pattern == 1L] <- par_full[names(par_full) == "lambda_free"]

  Beta <- matrix(0, n_latent, n_latent,
                 dimnames = list(latent_names, latent_names))
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
    method = "Laplace",
    missing = "listwise",
    loadings = Lambda,
    structural = Beta,
    latent_residual_sd = psi_sd,
    indicator_residual_sd = theta_sd,
    intercepts = setNames(unname(par_full[names(par_full) == "nu"]),
                          indicator_names),
    covariates = character(0),
    factor_scores = eta_hat,
    logLik = -opt$objective,
    AIC = 2 * opt$objective + 2 * n_params,
    BIC = 2 * opt$objective + log(n_obs) * n_params,
    n_params = n_params,
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
  cat("Observations:", x$n_obs)
  if (identical(x$missing, "fiml")) cat("  (missing: FIML)")
  cat("\n\n")
  cat("Loadings (marker indicators fixed at 1):\n")
  print(round(x$loadings, 4))
  if (any(x$structural != 0)) {
    cat("\nStructural paths (row ~ column):\n")
    print(round(x$structural, 4))
  }
  if (!is.null(x$latent_covariance)) {
    cat("\nLatent (co)variances:\n")
    print(round(x$latent_covariance, 4))
  } else {
    cat("\nLatent residual SDs:\n")
    print(round(x$latent_residual_sd, 4))
  }
  if (!is.null(x$fit_measures)) {
    fm <- x$fit_measures
    cat(sprintf(
      "\nchisq(%d) = %.3f, p = %.3f | CFI = %.3f | TLI = %.3f\n",
      as.integer(fm["df"]), fm["chisq"], fm["pvalue"], fm["cfi"],
      fm["tli"]))
    cat(sprintf("RMSEA = %.3f [%.3f, %.3f] | SRMR = %.3f\n",
                fm["rmsea"], fm["rmsea_ci_lower"], fm["rmsea_ci_upper"],
                fm["srmr"]))
  }
  cat("\nLog-likelihood:", round(x$logLik, 2), "\n")
  invisible(x)
}


#' @export
summary.gllamm_sem <- function(object, standardized = TRUE, ...) {
  print(object)
  if (!is.null(object$param_table)) {
    tab <- object$param_table
    if (standardized && !is.null(object$standardized)) {
      tab$est_std <- object$standardized$est_std
    }
    cat("\nParameter estimates:\n")
    tab_print <- tab
    tab_print[, -1] <- round(tab_print[, -1], 4)
    print(tab_print, row.names = FALSE)
  }
  invisible(object)
}
