#' GLLAMM model object class
#'
#' S3 class for fitted GLLAMM models
#'
#' @param x A gllamm object
#' @param object A gllamm object
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

  cat("Family:", object$family$family, "\n")
  cat("Link:", object$family$link, "\n\n")

  # Fixed effects table
  cat("Fixed effects:\n")
  fe <- cbind(
    Estimate = object$coefficients$fixed,
    `Std. Error` = sqrt(diag(object$vcov$fixed)),
    `z value` = object$coefficients$fixed / sqrt(diag(object$vcov$fixed))
  )
  fe <- cbind(fe, `Pr(>|z|)` = 2 * pnorm(-abs(fe[, "z value"])))
  printCoefmat(fe, digits = 4, signif.stars = TRUE)
  cat("\n")

  # Random effects
  cat("Random effects:\n")
  for (i in seq_along(object$random_terms)) {
    rt <- object$random_terms[[i]]
    cat("  Groups:", paste(rt$grouping, collapse = "/"), "\n")
    cat("    Number of groups:", object$n_groups[i], "\n")
    cat("    Variance:", round(object$coefficients$random_var[[i]], 4), "\n")
    cat("    Std.Dev.:", round(sqrt(object$coefficients$random_var[[i]]), 4), "\n")
  }
  cat("\n")

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


#' Extract random effects
#'
#' @param object A gllamm object
#' @param ... Additional arguments (ignored)
#'
#' @return List of random effects by group
#' @export
ranef.gllamm <- function(object, ...) {
  if (is.null(object$random_effects)) {
    stop("ranef is only available for multi-level models. ",
         "Model does not contain random effects.")
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
#' @param object A gllamm object
#' @param ... Additional arguments (ignored)
#'
#' @return List of variance-covariance matrices for random effects
#' @export
VarCorr.gllamm <- function(object, ...) {
  vc <- object$coefficients$random_var
  if (is.null(vc) || length(vc) == 0) {
    stop("VarCorr is only available for multi-level models. ",
         "Model does not contain random effects.")
  }
  if (!is.null(object$random_terms)) {
    names(vc) <- sapply(object$random_terms,
                        function(rt) paste(rt$grouping, collapse = "/"))
  }
  class(vc) <- "VarCorr.gllamm"
  vc
}


#' Generic VarCorr
#' @param object A fitted model object
#' @param ... Additional arguments
#' @export
VarCorr <- function(object, ...) {
  UseMethod("VarCorr")
}


# (print.VarCorr.gllamm lives in multilevel_methods.R and handles both the
# list form returned by VarCorr.gllamm and the data-frame form returned by
# VarCorr.gllamm_irt_multilevel)
