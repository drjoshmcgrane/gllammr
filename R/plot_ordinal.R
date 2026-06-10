#' Plotting Functions for Ordinal Regression Models
#'
#' Visualize ordinal regression model results
#'
#' @name plot_ordinal
NULL


#' Plot Ordinal Regression Model Diagnostics
#'
#' Create diagnostic plots for ordinal models including cumulative probabilities,
#' category probabilities, threshold parameters, and covariate effects
#'
#' @param x A gllamm_ordinal object
#' @param which Which plots to produce (1=cumulative, 2=category, 3=thresholds, 4=effects)
#' @param covariate Name of covariate to plot (default: first non-intercept covariate)
#' @param covariate_values Optional vector of covariate values to plot (default: -2 to 2)
#' @param ... Additional arguments passed to plotting functions
#'
#' @details
#' Plot types:
#' \itemize{
#'   \item \code{which = 1}: Cumulative Probabilities P(Y <= k) vs covariate
#'   \item \code{which = 2}: Category Probabilities P(Y = k) vs covariate
#'   \item \code{which = 3}: Threshold Parameters on latent scale
#'   \item \code{which = 4}: Covariate Effects (shows non-proportional effects for PPO)
#' }
#'
#' @examples
#' \dontrun{
#' # Fit ordinal model
#' fit <- fit_ordinal(rating ~ temp + contact + (1 | judge),
#'                    data = wine, link = "logit")
#'
#' # Plot all diagnostics for 'temp' covariate
#' plot(fit, which = 1:4, covariate = "temp")
#'
#' # Plot only cumulative probabilities
#' plot(fit, which = 1, covariate = "contact")
#' }
#'
#' @export
plot.gllamm_ordinal <- function(x, which = 1:3, covariate = NULL,
                                covariate_values = NULL, ...) {

  # Determine which covariate to plot
  if (is.null(covariate)) {
    # Get covariate names (different structure for PPO vs other links)
    if (x$link == "ppo") {
      covar_names <- colnames(x$coefficients$beta_ppo)
    } else {
      covar_names <- names(x$coefficients$fixed)
    }

    covariate <- covar_names[covar_names != "(Intercept)"][1]
    if (length(covariate) == 0 || is.na(covariate)) {
      stop("No covariates found in model. Specify covariate argument.")
    }
  }

  # Check that covariate exists
  if (x$link == "ppo") {
    if (!covariate %in% colnames(x$coefficients$beta_ppo)) {
      stop("Covariate '", covariate, "' not found in model.")
    }
  } else {
    if (!covariate %in% names(x$coefficients$fixed)) {
      stop("Covariate '", covariate, "' not found in model.")
    }
  }

  # Set up plotting area if multiple plots
  if (length(which) > 1) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)

    if (length(which) == 2) {
      par(mfrow = c(1, 2))
    } else if (length(which) == 3) {
      par(mfrow = c(2, 2))
    } else {
      par(mfrow = c(2, 2))
    }
  }

  # Generate each requested plot
  for (w in which) {
    if (w == 1) plot_cumulative_probs_ordinal(x, covariate, covariate_values, ...)
    if (w == 2) plot_category_probs_ordinal(x, covariate, covariate_values, ...)
    if (w == 3) plot_thresholds_ordinal(x, ...)
    if (w == 4) plot_covariate_effects_ordinal(x, covariate, ...)
  }

  invisible(x)
}


#' Plot Cumulative Probabilities
#' @keywords internal
plot_cumulative_probs_ordinal <- function(x, covariate, covariate_values, ...) {

  # Get covariate values
  if (is.null(covariate_values)) {
    covariate_values <- seq(-2, 2, length.out = 100)
  }

  # Get parameters
  thresholds <- x$coefficients$thresholds
  n_categories <- x$n_categories

  # Get coefficients (different for PPO)
  if (x$link == "ppo") {
    beta_ppo <- x$coefficients$beta_ppo[, covariate]  # Vector of length n_categories - 1
  } else {
    beta <- x$coefficients$fixed[covariate]
  }

  # Compute cumulative probabilities for each threshold
  cum_probs <- matrix(NA, length(covariate_values), n_categories - 1)

  for (i in seq_along(covariate_values)) {
    covar_val <- covariate_values[i]

    for (k in 1:(n_categories - 1)) {
      # Get eta for this threshold
      if (x$link == "ppo") {
        eta <- beta_ppo[k] * covar_val
      } else {
        eta <- beta * covar_val
      }

      # Compute cumulative probability
      if (x$link %in% c("logit", "ppo")) {
        cum_probs[i, k] <- plogis(thresholds[k] - eta)
      } else if (x$link == "probit") {
        cum_probs[i, k] <- pnorm(thresholds[k] - eta)
      } else {
        # For other links, use logit as approximation
        cum_probs[i, k] <- plogis(thresholds[k] - eta)
      }
    }
  }

  # Plot
  colors <- rainbow(n_categories - 1)

  plot(covariate_values, rep(0, length(covariate_values)), type = "n",
       xlab = covariate,
       ylab = "Cumulative Probability",
       main = paste("P(Y <= k) vs", covariate),
       ylim = c(0, 1),
       las = 1)

  grid(col = "gray90")

  # Add horizontal reference lines
  abline(h = c(0.25, 0.5, 0.75), lty = 2, col = "gray70")

  # Plot each cumulative probability
  for (k in 1:(n_categories - 1)) {
    lines(covariate_values, cum_probs[, k],
          col = colors[k], lwd = 2)
  }

  # Add legend
  legend("right",
         legend = paste("P(Y <=", 1:(n_categories - 1), ")"),
         col = colors,
         lwd = 2,
         bty = "n")
}


#' Plot Category Probabilities
#' @keywords internal
plot_category_probs_ordinal <- function(x, covariate, covariate_values, ...) {

  # Get covariate values
  if (is.null(covariate_values)) {
    covariate_values <- seq(-2, 2, length.out = 100)
  }

  # Get parameters
  beta <- x$coefficients$fixed[covariate]
  thresholds <- x$coefficients$thresholds
  n_categories <- x$n_categories

  # Compute category probabilities
  cat_probs <- matrix(NA, length(covariate_values), n_categories)

  for (i in seq_along(covariate_values)) {
    covar_val <- covariate_values[i]
    eta <- beta * covar_val

    # Compute cumulative probabilities
    cum_probs <- numeric(n_categories - 1)
    for (k in 1:(n_categories - 1)) {
      if (x$link == "logit") {
        cum_probs[k] <- plogis(thresholds[k] - eta)
      } else if (x$link == "probit") {
        cum_probs[k] <- pnorm(thresholds[k] - eta)
      } else {
        cum_probs[k] <- plogis(thresholds[k] - eta)
      }
    }

    # Convert to category probabilities
    cat_probs[i, 1] <- cum_probs[1]
    for (k in 2:(n_categories - 1)) {
      cat_probs[i, k] <- cum_probs[k] - cum_probs[k - 1]
    }
    cat_probs[i, n_categories] <- 1 - cum_probs[n_categories - 1]
  }

  # Plot
  colors <- rainbow(n_categories)

  plot(covariate_values, rep(0, length(covariate_values)), type = "n",
       xlab = covariate,
       ylab = "Probability",
       main = paste("P(Y = k) vs", covariate),
       ylim = c(0, 1),
       las = 1)

  grid(col = "gray90")

  # Plot each category probability
  for (k in 1:n_categories) {
    lines(covariate_values, cat_probs[, k],
          col = colors[k], lwd = 2)
  }

  # Add legend
  if (!is.null(x$category_labels)) {
    legend_labels <- paste("Y =", x$category_labels)
  } else {
    legend_labels <- paste("Y =", 1:n_categories)
  }

  legend("topright",
         legend = legend_labels,
         col = colors,
         lwd = 2,
         bty = "n",
         ncol = if (n_categories > 5) 2 else 1)
}


#' Plot Threshold Parameters
#' @keywords internal
plot_thresholds_ordinal <- function(x, ...) {

  thresholds <- x$coefficients$thresholds
  n_thresholds <- length(thresholds)

  # Create plot
  plot(1:n_thresholds, thresholds,
       type = "b", pch = 19, cex = 1.5, col = "blue", lwd = 2,
       xlab = "Threshold (k)",
       ylab = "Latent Position",
       main = "Threshold Parameters",
       las = 1,
       xaxt = "n")

  # Add x-axis
  axis(1, at = 1:n_thresholds, labels = 1:n_thresholds)

  # Add grid
  grid(col = "gray90")

  # Add horizontal line at 0
  abline(h = 0, lty = 2, col = "red", lwd = 1.5)

  # Add labels with values
  text(1:n_thresholds, thresholds,
       labels = sprintf("%.2f", thresholds),
       pos = 4, offset = 0.5, cex = 0.9, col = "darkblue")

  # Add interpretation text
  mtext("Thresholds partition the latent scale",
        side = 3, line = 0.5, cex = 0.8, col = "gray40")
}


#' Plot Covariate Effects
#' @keywords internal
plot_covariate_effects_ordinal <- function(x, covariate, ...) {

  # Check if this is a PPO model (would have different effects per threshold)
  is_ppo <- (x$link == "ppo")

  if (is_ppo) {
    # Get PPO coefficients for this covariate (different for each threshold)
    beta_ppo <- x$coefficients$beta_ppo[, covariate]
    n_thresholds <- length(beta_ppo)

    plot(1:n_thresholds, beta_ppo, type = "b", pch = 19, cex = 1.5, col = "blue",
         lwd = 2,
         xlab = "Threshold",
         ylab = "Coefficient",
         main = paste("Non-Proportional Effect of", covariate),
         ylim = range(c(0, beta_ppo)) * 1.2,
         las = 1,
         xaxt = "n")

    # Add x-axis labels
    axis(1, at = 1:n_thresholds, labels = paste("τ", 1:n_thresholds, sep = ""))

    # Add grid
    grid(col = "gray90")

    # Add horizontal line at 0
    abline(h = 0, lty = 2, col = "red", lwd = 1.5)

    # Add value labels
    for (i in 1:n_thresholds) {
      text(i, beta_ppo[i], sprintf("%.3f", beta_ppo[i]),
           pos = 3, offset = 0.5, cex = 0.9, col = "darkblue")
    }

    # Add interpretation
    mtext("Different effects per threshold (partial proportional odds)",
          side = 3, line = 0.3, cex = 0.8, col = "gray40")

  } else {
    # Get single proportional coefficient
    beta <- x$coefficients$fixed[covariate]
    # Single proportional effect
    # Visualize as bar with confidence interval if available

    plot(1, beta, pch = 19, cex = 3, col = "blue",
         xlim = c(0.5, 1.5), ylim = range(c(0, beta)) * 1.2,
         xlab = "", ylab = "Coefficient",
         main = paste("Proportional Effect of", covariate),
         las = 1, xaxt = "n")

    # Add x-axis label
    axis(1, at = 1, labels = covariate)

    # Add grid
    grid(col = "gray90")

    # Add horizontal line at 0
    abline(h = 0, lty = 2, col = "red", lwd = 1.5)

    # Add value label
    text(1, beta, sprintf("%.3f", beta),
         pos = 3, offset = 1, cex = 1.2, font = 2, col = "darkblue")

    # Add interpretation
    direction <- if (beta > 0) "increases" else "decreases"
    mtext(paste("Higher", covariate, direction, "log-odds of higher categories"),
          side = 3, line = 0.5, cex = 0.8, col = "gray40")

    # Add confidence interval if available
    if (!is.null(x$tmb_sdr) && !inherits(x$tmb_sdr, "try-error")) {
      # Extract SE for this coefficient
      # This would require matching parameter names
      # For now, skip CI
    }
  }
}


#' Plot Ordinal Model Effects for Multiple Covariates
#'
#' Compare effects across multiple covariates in an ordinal model
#'
#' @param object A gllamm_ordinal object
#' @param covariates Vector of covariate names (default: all non-intercept)
#' @param sort_by Sort covariates by: "magnitude", "name", or "none"
#' @param ... Additional arguments
#'
#' @examples
#' \dontrun{
#' fit <- fit_ordinal(rating ~ temp + contact + (1 | judge), data = wine)
#' plot_ordinal_effects(fit)
#' }
#'
#' @export
plot_ordinal_effects <- function(object, covariates = NULL,
                                sort_by = c("magnitude", "name", "none"), ...) {

  sort_by <- match.arg(sort_by)

  # Get covariates
  if (is.null(covariates)) {
    covariates <- names(object$coefficients$fixed)
    covariates <- covariates[covariates != "(Intercept)"]
  }

  # Extract coefficients
  coefs <- object$coefficients$fixed[covariates]

  # Sort if requested
  if (sort_by == "magnitude") {
    coefs <- coefs[order(abs(coefs), decreasing = TRUE)]
  } else if (sort_by == "name") {
    coefs <- coefs[order(names(coefs))]
  }

  # Create barplot
  colors <- ifelse(coefs > 0, "steelblue", "coral")

  bp <- barplot(coefs,
                main = "Covariate Effects (Proportional Odds)",
                ylab = "Coefficient",
                las = 2,
                col = colors,
                border = "white")

  # Add horizontal line at 0
  abline(h = 0, lwd = 2)

  # Add grid
  grid(nx = NA, ny = NULL, col = "gray90")

  # Redraw bars on top of grid
  barplot(coefs, add = TRUE, col = colors, border = "white",
          axes = FALSE, las = 2)

  # Add value labels
  text(bp, coefs, sprintf("%.2f", coefs),
       pos = ifelse(coefs > 0, 3, 1),
       cex = 0.8, font = 2)
}
