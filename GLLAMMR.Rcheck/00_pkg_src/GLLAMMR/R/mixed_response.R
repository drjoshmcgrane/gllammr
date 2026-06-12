#' Fit Joint Models for Mixed Response Types
#'
#' Fits a joint model for up to three outcomes of different types
#' (gaussian, binomial, poisson) measured on the same observations and
#' sharing a common random effect. The shared random effect induces
#' association between the outcomes.
#'
#' @param formulas Named list of up to three fixed-effects formulas, with
#'   names among "gaussian", "binomial", "poisson", e.g.
#'   \code{list(gaussian = y1 ~ x1, binomial = y2 ~ x2)}.
#' @param random One-sided random-intercept formula, e.g. \code{~ (1 | group)}
#' @param data Data frame containing all outcome and covariate variables;
#'   rows must be complete for every supplied outcome
#' @param start Optional starting values
#' @param control Optimization control list
#'
#' @details
#' The shared random effect enters every outcome's linear predictor with
#' loading 1 (a common-intercept joint model). Outcome-specific loadings
#' are not yet supported.
#'
#' @return An object of class \code{gllamm_mixed}
#'
#' @examples
#' \dontrun{
#' fit <- fit_mixed(
#'   formulas = list(gaussian = biomarker ~ age,
#'                   binomial = event ~ age + treatment),
#'   random = ~ (1 | patient),
#'   data = d)
#' }
#'
#' @export
fit_mixed <- function(formulas, random, data, start = NULL, control = list()) {

  allowed <- c("gaussian", "binomial", "poisson")
  if (!is.list(formulas) || is.null(names(formulas)) ||
      !all(names(formulas) %in% allowed) || length(formulas) < 1) {
    stop("formulas must be a named list with names among: ",
         paste(allowed, collapse = ", "))
  }

  # ---- Shared grouping from the random formula ----
  re_term <- attr(terms(random), "term.labels")
  if (length(re_term) != 1 || !grepl("\\|", re_term)) {
    stop("random must contain exactly one term of the form (1 | group)")
  }
  rt <- parse_random_term(re_term, data)
  if (!identical(deparse(rt$formula), "~1")) {
    stop("Only a shared random intercept (1 | group) is currently supported")
  }
  group_factor <- factor(data[[rt$grouping[1]]])
  groups <- as.integer(group_factor) - 1L
  n_groups <- nlevels(group_factor)
  n_obs <- nrow(data)
  n_random <- 1L

  # ---- Per-outcome responses and design matrices ----
  get_part <- function(fam) {
    if (is.null(formulas[[fam]])) {
      return(list(present = 0L, y = numeric(0), X = matrix(0, 0, 0), beta0 = numeric(0)))
    }
    f <- formulas[[fam]]
    mf <- model.frame(f, data = data, na.action = na.fail)
    list(present = 1L,
         y = model.response(mf),
         X = model.matrix(f, data = data),
         name = as.character(f[[2]]))
  }
  p1 <- get_part("gaussian")
  p2 <- get_part("binomial")
  p3 <- get_part("poisson")

  if (p2$present == 1L && !all(p2$y %in% c(0, 1))) {
    stop("The binomial outcome must be coded 0/1")
  }

  tmb_data <- list(
    y1 = as.numeric(p1$y),
    y2 = as.integer(p2$y),
    y3 = as.numeric(p3$y),
    X1 = as.matrix(p1$X),
    X2 = as.matrix(p2$X),
    X3 = as.matrix(p3$X),
    Z = Matrix::Matrix(matrix(1, n_obs, 1), sparse = TRUE),
    groups = groups,
    n_groups = as.integer(n_groups),
    n1 = as.integer(length(p1$y)),
    n2 = as.integer(length(p2$y)),
    n3 = as.integer(length(p3$y)),
    n_fixed1 = as.integer(ncol(p1$X)),
    n_fixed2 = as.integer(ncol(p2$X)),
    n_fixed3 = as.integer(ncol(p3$X)),
    n_random = n_random,
    has_y1 = p1$present,
    has_y2 = p2$present,
    has_y3 = p3$present,
    model_name = "mixed_response"
  )

  if (is.null(start)) {
    tmb_params <- list(
      beta1 = rep(0, max(ncol(p1$X), 0)),
      beta2 = rep(0, max(ncol(p2$X), 0)),
      beta3 = rep(0, max(ncol(p3$X), 0)),
      u = rep(0, n_groups),
      log_sigma1 = 0,
      log_sigma_u = log(0.5),
      theta = 0
    )
  } else {
    tmb_params <- start
  }

  # Dead parameters: theta (single RE), log_sigma1 without a gaussian outcome,
  # and any absent outcome's beta block
  tmb_map <- list(theta = factor(NA))
  if (p1$present == 0L) {
    tmb_map$log_sigma1 <- factor(NA)
    tmb_map$beta1 <- factor(rep(NA, length(tmb_params$beta1)))
  }
  if (p2$present == 0L) tmb_map$beta2 <- factor(rep(NA, length(tmb_params$beta2)))
  if (p3$present == 0L) tmb_map$beta3 <- factor(rep(NA, length(tmb_params$beta3)))

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

  coefs <- list()
  if (p1$present == 1L) {
    b <- par_full[names(par_full) == "beta1"]; names(b) <- colnames(p1$X)
    coefs$gaussian <- b
  }
  if (p2$present == 1L) {
    b <- par_full[names(par_full) == "beta2"]; names(b) <- colnames(p2$X)
    coefs$binomial <- b
  }
  if (p3$present == 1L) {
    b <- par_full[names(par_full) == "beta3"]; names(b) <- colnames(p3$X)
    coefs$poisson <- b
  }

  sigma1 <- if (p1$present == 1L) {
    exp(unname(par_full[names(par_full) == "log_sigma1"]))
  } else NA_real_
  sigma_u <- exp(unname(par_full[names(par_full) == "log_sigma_u"]))

  n_params <- length(unlist(coefs)) + 1 + (p1$present == 1L)

  result <- list(
    coefficients = coefs,
    residual_sd = sigma1,
    random_sd = sigma_u,
    random_effects = unname(par_full[names(par_full) == "u"]),
    outcomes = names(formulas),
    logLik = -opt$objective,
    AIC = 2 * opt$objective + 2 * n_params,
    BIC = 2 * opt$objective + log(n_obs) * n_params,
    convergence = list(converged = (opt$convergence == 0),
                       message = opt$message),
    n_obs = n_obs,
    n_groups = n_groups,
    formulas = formulas,
    random = random,
    data = data,
    tmb_obj = obj,
    tmb_opt = opt,
    tmb_sdr = sdr
  )
  class(result) <- c("gllamm_mixed", "gllamm")
  result
}


#' @export
print.gllamm_mixed <- function(x, ...) {
  cat("Joint Mixed-Response Model (shared random intercept)\n\n")
  cat("Observations:", x$n_obs, " Groups:", x$n_groups, "\n")
  cat("Outcomes:", paste(x$outcomes, collapse = ", "), "\n\n")
  for (nm in names(x$coefficients)) {
    cat(nm, "coefficients:\n")
    print(round(x$coefficients[[nm]], 4))
  }
  if (!is.na(x$residual_sd)) {
    cat("\nGaussian residual SD:", round(x$residual_sd, 4), "\n")
  }
  cat("Shared random-effect SD:", round(x$random_sd, 4), "\n")
  cat("Log-likelihood:", round(x$logLik, 2), "\n")
  invisible(x)
}