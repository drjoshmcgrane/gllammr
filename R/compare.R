#' Compare fitted gllammr models
#'
#' Builds a comparison table for any set of fitted gllammr models (or any
#' objects with \code{logLik}, \code{AIC}, \code{BIC}, and \code{n_obs}
#' components): log-likelihood, parameter count, AIC/BIC with deltas, and
#' Akaike weights. Works across model classes - GLMMs, IRT, latent class,
#' CDM, SEM, survival - because every fitter reports the same marginal
#' quantities.
#'
#' @param ... Named fitted model objects (names label the table rows;
#'   unnamed arguments are labelled by their expressions)
#' @param sort_by "none" (default; preserve input order), "AIC", or "BIC"
#'
#' @return An object of class \code{gllammr_model_comparison}: a data
#'   frame with one row per model.
#'
#' @details
#' Information criteria are only meaningful across models fitted to the
#' \emph{same response data}; the function checks that \code{n_obs}
#' agrees and warns otherwise, but cannot verify that the responses
#' themselves coincide - that is the analyst's responsibility. For
#' inequality-constrained models (ordered/partially ordered classes,
#' monotone CDMs) the parameter count is nominal and likelihood-ratio
#' comparisons have chi-bar-square null distributions; treat those
#' comparisons descriptively.
#'
#' @examples
#' \dontrun{
#' compare_models(rasch = fit_irt(Y, "Rasch"),
#'                twopl = fit_irt(Y, "2PL"),
#'                lca3  = fit_lca(Y, nclass = 3))
#' }
#'
#' @export
compare_models <- function(..., sort_by = c("none", "AIC", "BIC")) {
  sort_by <- match.arg(sort_by)
  fits <- list(...)
  if (length(fits) < 2) {
    stop("compare_models() needs at least two fitted models")
  }
  labels <- names(fits)
  exprs <- vapply(substitute(list(...))[-1], deparse, character(1))
  if (is.null(labels)) labels <- exprs
  labels[labels == ""] <- exprs[labels == ""]

  get_num <- function(f, what) {
    v <- f[[what]]
    if (is.null(v) || !is.numeric(v)) NA_real_ else as.numeric(v[1])
  }
  ll <- vapply(fits, get_num, 0, what = "logLik")
  aic <- vapply(fits, get_num, 0, what = "AIC")
  bic <- vapply(fits, get_num, 0, what = "BIC")
  nobs <- vapply(fits, function(f) {
    v <- f[["n_obs"]] %||% f[["n_persons"]]
    if (is.null(v)) NA_real_ else as.numeric(v[1])
  }, 0)
  if (anyNA(ll) || anyNA(aic)) {
    stop("Every model must provide logLik and AIC components")
  }
  if (length(unique(stats::na.omit(nobs))) > 1) {
    warning("Models were fitted to different numbers of observations (",
            paste(nobs, collapse = ", "),
            "); information criteria are not comparable")
  }

  k <- (aic + 2 * ll) / 2          # nominal parameter count
  d_aic <- aic - min(aic)
  d_bic <- bic - min(bic)
  w <- exp(-d_aic / 2)
  w <- w / sum(w)

  tab <- data.frame(
    model = labels,
    class = vapply(fits, function(f) class(f)[1], character(1)),
    logLik = ll,
    n_params = k,
    AIC = aic,
    dAIC = d_aic,
    BIC = bic,
    dBIC = d_bic,
    akaike_weight = w,
    stringsAsFactors = FALSE
  )
  if (sort_by != "none") {
    tab <- tab[order(tab[[sort_by]]), , drop = FALSE]
  }
  rownames(tab) <- NULL
  structure(tab, fits = setNames(fits, labels),
            class = c("gllammr_model_comparison", "data.frame"))
}


#' @export
print.gllammr_model_comparison <- function(x, digits = 2, ...) {
  cat("Model comparison (", nrow(x), " models)\n\n", sep = "")
  tab <- as.data.frame(x)
  num <- vapply(tab, is.numeric, TRUE)
  tab[num] <- lapply(tab[num], round, digits)
  print(tab, row.names = FALSE)
  cat("\nBest by AIC:", x$model[which.min(x$AIC)],
      "| by BIC:", x$model[which.min(x$BIC)], "\n")
  invisible(x)
}
