#' GLLAMM model object class
#'
#' S3 class for fitted GLLAMM models
#'
#' @param x A gllamm object
#' @param object A gllamm object
#' @param which Which covariance block to return: "fixed" (default) or "all"
#' @param type Covariance type: "model" (default) or "sandwich"
#'   (cluster-robust)
#' @param ... Additional arguments
#'
#' @name gllamm-class
NULL


#' @export
#' @rdname gllamm-class
print.gllamm <- function(x, ...) {
  cat("Generalized Linear Latent and Mixed Model\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Family:", x$family$family, "\n")
  cat("Link:", x$family$link, "\n\n")

  cat("Random effects:\n")
  for (i in seq_along(x$random_terms)) {
    rt <- x$random_terms[[i]]
    cat("  Groups:", paste(rt$grouping, collapse = "/"), "\n")
    cat("    Number of groups:", x$n_groups[i], "\n")
  }
  cat("\n")

  cat("Fixed effects:\n")
  print(round(x$coefficients$fixed, 4))
  cat("\n")

  cat("Log-likelihood:", round(x$logLik, 2), "\n")
  cat("AIC:", round(x$AIC, 2), "  BIC:", round(x$BIC, 2), "\n")

  if (!x$convergence$converged) {
    cat("\nWarning: Model did not converge!\n")
  }

  invisible(x)
}


#' @export
#' @rdname gllamm-class
summary.gllamm <- function(object, ...) {
  cat("Generalized Linear Latent and Mixed Model\n\n")
  cat("Call:\n")
  print(object$call)
  cat("\n")

  if (!is.null(object$family$family)) {
    cat("Family:", object$family$family, "\n")
    cat("Link:", object$family$link, "\n\n")
  }

  # Fixed effects table (SEs only when a matching vcov is available -
  # never recycle a mismatched matrix into wrong standard errors)
  est <- object$coefficients$fixed
  if (!is.null(est) && length(est)) {
    vf <- if (is.list(object$vcov)) object$vcov$fixed else object$vcov
    se <- rep(NA_real_, length(est))
    if (is.matrix(vf) && nrow(vf) == length(est)) se <- sqrt(diag(vf))
    cat("Fixed effects:\n")
    fe <- cbind(Estimate = est, `Std. Error` = se, `z value` = est / se,
                `Pr(>|z|)` = 2 * pnorm(-abs(est / se)))
    printCoefmat(fe, digits = 4, signif.stars = TRUE, na.print = "-")
    cat("\n")
  }

  # Random effects
  if (!is.null(object$random_terms) && length(object$random_terms) &&
      !is.null(object$coefficients$random_var)) {
    cat("Random effects:\n")
    for (i in seq_along(object$random_terms)) {
      rt <- object$random_terms[[i]]
      cat("  Groups:", paste(rt$grouping, collapse = "/"), "\n")
      if (!is.null(object$n_groups)) {
        cat("    Number of groups:", object$n_groups[i], "\n")
      }
      rv <- object$coefficients$random_var[[i]]
      cat("    Variance:", round(rv, 4), "\n")
      cat("    Std.Dev.:", round(sqrt(rv), 4), "\n")
    }
    cat("\n")
  }

  # Model fit
  cat("Number of observations:", object$n_obs, "\n")
  cat("Log-likelihood:", round(object$logLik, 2), "\n")
  cat("AIC:", round(object$AIC, 2), "\n")
  cat("BIC:", round(object$BIC, 2), "\n")

  if (!object$convergence$converged) {
    cat("\nWarning: Model did not converge!\n")
    cat("  Message:", object$convergence$message, "\n")
  }

  invisible(object)
}


#' @export
#' @rdname gllamm-class
coef.gllamm <- function(object, ...) {
  object$coefficients
}


#' @export
#' @rdname gllamm-class
vcov.gllamm <- function(object, which = "fixed",
                        type = c("model", "sandwich"), ...) {
  type <- match.arg(type)
  if (type == "sandwich") {
    V <- sandwich_vcov_gllamm(object)
    return(attr(V, "fixed"))
  }
  # Models that store one parameter covariance matrix directly (e.g. SEM)
  if (is.matrix(object$vcov)) {
    return(object$vcov)
  }
  if (which == "fixed") {
    return(object$vcov$fixed)
  } else if (which == "all") {
    return(object$vcov$all)
  } else {
    stop("'which' must be 'fixed' or 'all'")
  }
}


#' @export
#' @rdname gllamm-class
logLik.gllamm <- function(object, ...) {
  structure(
    object$logLik,
    df = object$n_params,
    nobs = object$n_obs,
    class = "logLik"
  )
}


#' @export
#' @rdname gllamm-class
fitted.gllamm <- function(object, ...) {
  # Some fitters (e.g. fit_binomial) store fitted_values instead
  if (is.null(object$fitted.values)) object$fitted_values else object$fitted.values
}


#' @export
#' @rdname gllamm-class
residuals.gllamm <- function(object, type = c("response", "pearson", "deviance"), ...) {
  type <- match.arg(type)

  y <- object$y
  mu <- fitted(object)

  switch(type,
    response = y - mu,
    pearson = (y - mu) / sqrt(object$family$variance(mu)),
    deviance = {
      # Deviance residuals
      sign(y - mu) * sqrt(pmax(object$family$dev.resids(y, mu, rep(1, length(y))), 0))
    }
  )
}


#' Extract fixed effects
#'
#' @param object A gllamm object
#' @param ... Additional arguments (ignored)
#'
#' @return Named vector of fixed effects coefficients
#' @export
fixef.gllamm <- function(object, ...) {
  object$coefficients$fixed
}


#' Generic fixef
#' @param object A fitted model object
#' @param ... Additional arguments
#' @export
fixef <- function(object, ...) {
  UseMethod("fixef")
}


.no_re_hint <- function(object) {
  if (inherits(object, "gllamm_irt")) {
    return(paste0("For IRT fits: person abilities are in ",
                  "$person_abilities, the ability SD in $ability_sd."))
  }
  if (inherits(object, "gllamm_eirt")) {
    return(paste0("For explanatory IRT fits: person abilities are in ",
                  "$person_abilities; item residual SDs in $residual_sd."))
  }
  if (inherits(object, "gllamm_lca")) {
    return(paste0("For latent class fits: class posteriors are in ",
                  "$posterior, prevalences in $class_probs."))
  }
  if (inherits(object, "gllamm_cdm")) {
    return(paste0("For CDM fits: attribute posteriors are in ",
                  "$attribute_posteriors, profile prevalences in ",
                  "$profile_probs."))
  }
  if (inherits(object, "gllamm_sem")) {
    return(paste0("For SEM fits: factor scores are in $factor_scores, ",
                  "latent (co)variances in $latent_covariance."))
  }
  if (inherits(object, "gllamm_npml")) {
    return(paste0("For NPML fits the latent distribution is discrete: ",
                  "mass-point locations are in $locations, masses in ",
                  "$masses."))
  }
  if (inherits(object, "gllamm_mixed")) {
    return(paste0("For mixed-response fits: the shared random-intercept ",
                  "SD is in $random_sd."))
  }
  "Model does not contain (normal) random effects."
}


#' Extract random effects
#'
#' @param object A gllamm object
#' @param ... Additional arguments (ignored)
#'
#' @return List of random effects by group
#' @export
ranef.gllamm <- function(object, ...) {
  if (is.null(object$random_effects)) {
    stop("ranef is only available for models with (normal) random ",
         "effects. ", .no_re_hint(object))
  }
  object$random_effects
}


#' Generic ranef
#' @param object A fitted model object
#' @param ... Additional arguments
#' @export
ranef <- function(object, ...) {
  UseMethod("ranef")
}


#' Extract variance components
#'
#' @param x A gllamm object
#' @param ... Additional arguments (ignored)
#'
#' @return List of variance-covariance matrices for random effects
#' @export
VarCorr.gllamm <- function(x, ...) {
  vc <- x$coefficients$random_var
  if (is.null(vc) || length(vc) == 0) {
    stop("VarCorr is only available for models with (normal) random ",
         "effects. ", .no_re_hint(x))
  }
  if (!is.null(x$random_terms)) {
    names(vc) <- sapply(x$random_terms,
                        function(rt) paste(rt$grouping, collapse = "/"))
  }
  class(vc) <- "VarCorr.gllamm"
  vc
}


# (print.VarCorr.gllamm lives in multilevel_methods.R and handles both the
# list form returned by VarCorr.gllamm and the data-frame form returned by
# VarCorr.gllamm_irt_multilevel)
