#' Shared robustness helpers for model fitting
#'
#' Internal utilities that give consistent, informative diagnostics across
#' every fitter: validating a TMB sdreport, warning on non-convergence,
#' guarding nlminb optimization, and guarding matrix inversion.
#'
#' @name fit_checks
#' @keywords internal
#' @noRd
NULL


#' Validate a TMB sdreport result
#'
#' Input is the value of \code{try(TMB::sdreport(...), silent = TRUE)}. Emits
#' an informative warning when standard errors could not be computed or when
#' the Hessian is not positive definite, and reports whether the standard
#' errors are trustworthy.
#'
#' @param sdr Result of \code{try(TMB::sdreport(obj), silent = TRUE)}, or
#'   \code{NULL} when standard errors were not requested.
#' @param context Short label describing the model, used in warnings.
#' @return A list with elements \code{sdr} (the original object, or \code{NULL}
#'   on failure) and \code{se_ok} (logical).
#' @keywords internal
#' @noRd
check_sdreport <- function(sdr, context = "model") {
  # Standard errors were not requested (sdr deliberately NULL): nothing to say.
  if (is.null(sdr)) {
    return(list(sdr = NULL, se_ok = FALSE))
  }
  if (inherits(sdr, "try-error")) {
    msg <- conditionMessage(attr(sdr, "condition") %||% simpleError(""))
    warning("Standard error computation failed for ", context, ": ", msg,
            call. = FALSE)
    return(list(sdr = NULL, se_ok = FALSE))
  }
  if (!isTRUE(sdr$pdHess)) {
    warning("Hessian not positive definite for ", context,
            "; standard errors are unreliable. The model may be ",
            "over-parameterized or a variance component may be near zero.",
            call. = FALSE)
    return(list(sdr = sdr, se_ok = FALSE))
  }
  list(sdr = sdr, se_ok = TRUE)
}


#' Warn when a fitted model failed to converge
#'
#' Safe to call on any fitted object: looks for a convergence flag at
#' \code{object$converged} or \code{object$convergence$converged} and warns
#' only when it is explicitly \code{FALSE}.
#'
#' @param object A fitted model object.
#' @return Invisibly \code{NULL}; called for its warning side effect.
#' @keywords internal
#' @noRd
warn_not_converged <- function(object) {
  conv <- object$converged
  if (is.null(conv)) conv <- object$convergence$converged
  if (isFALSE(conv)) {
    warning("Model did not converge; estimates may be unreliable.",
            call. = FALSE)
  }
  invisible(NULL)
}


#' Guarded nlminb optimization
#'
#' Wraps \code{stats::nlminb} so optimization failures (errors or a
#' non-finite objective) raise a single informative error instead of an
#' opaque one. Preserves nlminb's return shape on success.
#'
#' @param start,objective,gradient Passed to \code{stats::nlminb}.
#' @param ... Further arguments to \code{stats::nlminb} (e.g. \code{lower},
#'   \code{upper}, \code{control}).
#' @param context Short label describing the model, used in the error.
#' @return The \code{stats::nlminb} result list.
#' @keywords internal
#' @noRd
safe_nlminb <- function(start, objective, gradient = NULL, ...,
                        context = "model") {
  opt <- tryCatch(
    stats::nlminb(start = start, objective = objective,
                  gradient = gradient, ...),
    error = function(e) {
      stop("Optimization failed for ", context, ": ", conditionMessage(e),
           ". Try different starting values or simplify the model.",
           call. = FALSE)
    }
  )
  ov <- opt$objective
  if (is.null(ov) || !is.finite(ov)) {
    stop("Optimization failed for ", context,
         ": the objective is not finite. Try different starting values ",
         "or simplify the model.", call. = FALSE)
  }
  opt
}


#' Guarded matrix inversion
#'
#' \code{tryCatch(solve(M))} that warns and returns \code{NULL} on failure
#' (e.g. a singular matrix) rather than erroring. Callers must handle the
#' \code{NULL} return.
#'
#' @param M A square matrix to invert.
#' @param context Short label describing the matrix, used in the warning.
#' @return The inverse of \code{M}, or \code{NULL} on failure.
#' @keywords internal
#' @noRd
safe_solve <- function(M, context = "matrix") {
  tryCatch(
    solve(M),
    error = function(e) {
      warning("Could not invert ", context, ": ", conditionMessage(e),
              ". Returning NULL.", call. = FALSE)
      NULL
    }
  )
}
