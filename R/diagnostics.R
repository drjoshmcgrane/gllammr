#' Diagnostic Methods for GLLAMM Models
#'
#' Comprehensive diagnostics for model checking
#'
#' @name diagnostics
NULL


#' Plot diagnostics for GLLAMM models
#'
#' @param x A gllamm object
#' @param which Which plots to produce (1-6)
#' @param ... Additional arguments
#'
#' @export
plot.gllamm <- function(x, which = c(1, 2, 3, 5), ...) {

  # Dispatch to model-specific plotting functions
  if (inherits(x, "gllamm_irt")) {
    return(plot.gllamm_irt(x, which = which, ...))
  }
  if (inherits(x, "gllamm_lca")) {
    return(plot.gllamm_lca(x, which = which, ...))
  }
  if (inherits(x, "gllamm_ordinal")) {
    return(plot.gllamm_ordinal(x, which = which, ...))
  }

  # Default GLMM diagnostics for standard models
  # Calculate residuals and fitted values
  resids <- residuals(x, type = "response")
  fitted_vals <- fitted(x)

  # Standardized residuals
  std_resids <- resids / sd(resids)

  # Set up plotting area
  if (length(which) > 1) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)

    n_plots <- length(which)
    n_rows <- ceiling(sqrt(n_plots))
    n_cols <- ceiling(n_plots / n_rows)
    par(mfrow = c(n_rows, n_cols))
  }

  for (w in which) {
    if (w == 1) {
      # Residuals vs Fitted
      plot(fitted_vals, resids,
           xlab = "Fitted values",
           ylab = "Residuals",
           main = "Residuals vs Fitted",
           pch = 20, col = "darkgray")
      abline(h = 0, lty = 2, col = "red")
      lines(lowess(fitted_vals, resids), col = "blue", lwd = 2)

    } else if (w == 2) {
      # Normal Q-Q plot
      qqnorm(std_resids, main = "Normal Q-Q Plot",
             pch = 20, col = "darkgray")
      qqline(std_resids, col = "red", lty = 2)

    } else if (w == 3) {
      # Scale-Location (sqrt standardized residuals vs fitted)
      sqrt_std_resids <- sqrt(abs(std_resids))
      plot(fitted_vals, sqrt_std_resids,
           xlab = "Fitted values",
           ylab = expression(sqrt("|Standardized residuals|")),
           main = "Scale-Location",
           pch = 20, col = "darkgray")
      lines(lowess(fitted_vals, sqrt_std_resids), col = "blue", lwd = 2)

    } else if (w == 4) {
      # Cook's distance at the group level (case deletion per cluster)
      D <- tryCatch(cooks.distance(x), error = function(e) {
        warning("Cook's distance unavailable: ", conditionMessage(e))
        NULL
      })
      if (!is.null(D)) {
        barplot(D,
                xlab = "Group",
                ylab = "Cook's distance",
                main = "Group-Level Cook's Distance",
                col = "steelblue",
                border = NA)
        abline(h = 4 / length(D), lty = 2, col = "red")
      }

    } else if (w == 5) {
      # Residuals vs Leverage (not standard for mixed models)
      # Instead: Residuals by group (first grouping factor)
      if (length(x$random_terms) > 0 && !is.null(x$data)) {
        rt <- x$random_terms[[1]]
        gvars <- rt$grouping_vars
        if (is.null(gvars)) gvars <- rt$grouping
        if (length(gvars) > 1) {
          group_factor <- interaction(x$data[, gvars], drop = TRUE)
        } else {
          group_factor <- factor(x$data[[gvars]])
        }

        boxplot(resids ~ group_factor,
                xlab = "Group",
                ylab = "Residuals",
                main = "Residuals by Group",
                col = "lightblue",
                outline = FALSE)
        abline(h = 0, lty = 2, col = "red")
      }

    } else if (w == 6) {
      # Random effects distribution
      if (length(x$random_effects) > 0) {
        re_vals <- unlist(x$random_effects)
        hist(re_vals,
             main = "Distribution of Random Effects",
             xlab = "Random Effects",
             col = "lightgreen",
             breaks = 20)
        curve(dnorm(x, 0, sd(re_vals)) * length(re_vals) * diff(range(re_vals))/20,
              add = TRUE, col = "red", lwd = 2)
      }
    }
  }

  invisible(x)
}


#' Influence diagnostics for GLLAMM
#'
#' @param object A gllamm object
#' @param ... Additional arguments
#'
#' @return Data frame with influence measures
#'
#' @export
influence.gllamm <- function(object, ...) {

  n <- object$n_obs
  resids <- residuals(object, type = "response")
  fitted_vals <- fitted(object)

  # Hat values (approximate for mixed models)
  # Use diagonal of projection matrix approximation
  hat_values <- rep(NA, n)

  # Standardized residuals
  std_resids <- resids / sd(resids)

  # Create influence data frame
  infl <- data.frame(
    obs = 1:n,
    residual = resids,
    std_residual = std_resids,
    fitted = fitted_vals,
    hat = hat_values
  )

  class(infl) <- c("influence.gllamm", "data.frame")

  return(infl)
}


#' Outlier detection for GLLAMM
#'
#' @param object A gllamm object
#' @param threshold Threshold for standardized residuals (default: 3)
#' @param ... Additional arguments
#'
#' @return Indices of potential outliers
#'
#' @export
find_outliers <- function(object, threshold = 3, ...) {

  resids <- residuals(object, type = "response")
  std_resids <- resids / sd(resids)

  outliers <- which(abs(std_resids) > threshold)

  if (length(outliers) > 0) {
    message("Found ", length(outliers), " potential outliers (|std residual| > ", threshold, ")")

    outlier_info <- data.frame(
      index = outliers,
      residual = resids[outliers],
      std_residual = std_resids[outliers],
      fitted = fitted(object)[outliers]
    )

    return(outlier_info)
  } else {
    message("No outliers detected")
    return(NULL)
  }
}


#' Goodness-of-fit tests for GLLAMM
#'
#' @param object A gllamm object
#' @param ... Additional arguments
#'
#' @export
gof.gllamm <- function(object, ...) {

  cat("Goodness of Fit\n")
  cat("===============\n\n")

  # Log-likelihood and information criteria
  cat("Log-likelihood:", round(object$logLik, 2), "\n")
  cat("AIC:", round(object$AIC, 2), "\n")
  cat("BIC:", round(object$BIC, 2), "\n\n")

  # Residual diagnostics
  resids <- residuals(object, type = "response")

  cat("Residual Summary:\n")
  cat("  Min:", round(min(resids), 3), "\n")
  cat("  Q1:", round(quantile(resids, 0.25), 3), "\n")
  cat("  Median:", round(median(resids), 3), "\n")
  cat("  Q3:", round(quantile(resids, 0.75), 3), "\n")
  cat("  Max:", round(max(resids), 3), "\n")
  cat("  SD:", round(sd(resids), 3), "\n\n")

  # Normality test (Shapiro-Wilk if n < 5000)
  if (length(resids) < 5000) {
    sw_test <- shapiro.test(resids)
    cat("Shapiro-Wilk normality test:\n")
    cat("  W =", round(sw_test$statistic, 4), "\n")
    cat("  p-value =", format.pval(sw_test$p.value), "\n")
    if (sw_test$p.value < 0.05) {
      cat("  (Residuals may deviate from normality)\n")
    }
  }

  invisible(object)
}


#' Variance decomposition (ICC)
#'
#' @param x A gllamm object
#' @param quiet Suppress printed output and messages (default: FALSE)
#' @param ... Additional arguments
#'
#' @export
icc.gllamm <- function(x, quiet = FALSE, ...) {
  cat <- if (quiet) function(...) invisible(NULL) else base::cat
  message <- if (quiet) function(...) invisible(NULL) else base::message
  print <- if (quiet) function(...) invisible(NULL) else base::print

  if (length(x$coefficients$random_var) == 0) {
    stop("icc is only available for multi-level models. ",
         "Model does not contain random effects.")
  }

  # Random-effect variances: random_var is a list of per-term covariance
  # matrices; take the diagonal (variances) of each
  .re_variances <- function(rv) {
    unlist(lapply(rv, function(m) {
      if (is.matrix(m)) diag(m) else as.numeric(m)
    }))
  }

  fam <- x$family$family
  if (is.null(fam)) {
    fam <- if (inherits(x, "gllamm_binomial")) "binomial" else "gaussian"
  }

  # For Gaussian models
  if (fam == "gaussian") {
    # Residual variance
    resid_var <- var(residuals(x))

    # Random effects variance
    re_vars <- .re_variances(x$coefficients$random_var)

    # Total variance
    total_var <- sum(re_vars) + resid_var

    # ICC for each level
    iccs <- re_vars / total_var

    names(iccs) <- paste0("Level", seq_along(iccs))

    cat("Intraclass Correlation Coefficients (ICC):\n")
    for (i in seq_along(iccs)) {
      cat(" ", names(iccs)[i], ":", round(iccs[i], 4), "\n")
    }

    cat("\nProportion of variance:\n")
    cat("  Random effects:", round(sum(re_vars) / total_var, 4), "\n")
    cat("  Residual:", round(resid_var / total_var, 4), "\n")

    return(invisible(iccs))

  } else {
    message("ICC computation for non-Gaussian families uses approximation")

    # Approximate ICC for GLMMs
    # Use latent variable formulation
    re_vars <- .re_variances(x$coefficients$random_var)

    if (fam == "binomial") {
      # Logistic: residual variance = pi^2/3
      resid_var <- pi^2/3
    } else if (fam == "poisson") {
      # Use overdispersion approximation
      resid_var <- 1  # Approximate
    }

    total_var <- sum(re_vars) + resid_var
    iccs <- re_vars / total_var

    names(iccs) <- paste0("Level", seq_along(iccs))

    cat("Approximate ICC (", fam, "):\n", sep = "")
    print(round(iccs, 4))

    return(invisible(iccs))
  }
}


#' Group-level Cook's distance for GLLAMM models
#'
#' Influence of each cluster on the fixed effects, computed by refitting the
#' model with the cluster deleted: D_j = (beta - beta_(-j))' V^(-1)
#' (beta - beta_(-j)) / p, with V the estimated covariance of the fixed
#' effects. Case deletion at the cluster level is the standard influence
#' measure for mixed models; observation-level deletion would break the
#' random-effects structure.
#'
#' @param model A fitted \code{gllamm} object (from \code{gllamm()})
#' @param max_groups Refuse to run for more clusters than this (each cluster
#'   costs one model refit); raise the limit explicitly for large data
#' @param ... Additional arguments (currently unused)
#'
#' @return Named vector of Cook's distances, one per cluster. Clusters whose
#'   deletion refit fails get \code{NA}.
#'
#' @export
cooks.distance.gllamm <- function(model, max_groups = 50, ...) {
  if (is.null(model$formula) || is.null(model$data) || is.null(model$family)) {
    stop("Cook's distance requires a model fitted via gllamm() ",
         "(formula, data, and family must be stored on the object)")
  }

  parsed <- parse_formula(model$formula, model$data)
  model_data <- make_model_matrices(parsed, model$data)
  groups <- model_data$groups[[1]]      # 0-indexed
  n_groups <- model_data$n_groups[1]

  if (n_groups > max_groups) {
    stop("Model has ", n_groups, " clusters; each one costs a refit. ",
         "Raise max_groups (currently ", max_groups, ") to proceed.")
  }

  beta_full <- model$coefficients$fixed
  p <- length(beta_full)
  V <- model$vcov$fixed
  if (is.null(V) || anyNA(V)) {
    stop("Fixed-effects covariance unavailable; cannot scale Cook's distance")
  }
  V_inv <- solve(V)

  # Tight iteration caps: refits start near the full-data optimum
  refit_control <- list(iter.max = 200, eval.max = 400)

  D <- rep(NA_real_, n_groups)
  for (j in seq_len(n_groups)) {
    keep <- groups != (j - 1)
    refit <- try(
      gllamm(model$formula, data = model$data[keep, , drop = FALSE],
             family = model$family, control = refit_control),
      silent = TRUE
    )
    if (inherits(refit, "try-error")) next
    delta <- beta_full - refit$coefficients$fixed
    D[j] <- as.numeric(t(delta) %*% V_inv %*% delta) / p
  }

  names(D) <- paste0("Group", seq_len(n_groups))
  D
}
