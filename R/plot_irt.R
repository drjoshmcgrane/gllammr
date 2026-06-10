#' Plotting Functions for IRT Models
#'
#' Visualize item response theory models
#'
#' @name plot_irt
NULL


#' Plot IRT Model Diagnostics
#'
#' Create diagnostic plots for IRT models including item characteristic curves,
#' item information functions, test information, and ability distributions
#'
#' @param x A gllamm_irt object
#' @param which Which plots to produce (1=ICC, 2=IIF, 3=TIF, 4=Ability)
#' @param items Which items to plot (default: first 6 items)
#' @param ... Additional arguments passed to plotting functions
#'
#' @details
#' Plot types:
#' \itemize{
#'   \item \code{which = 1}: Item Characteristic Curves (ICC)
#'   \item \code{which = 2}: Item Information Functions (IIF)
#'   \item \code{which = 3}: Test Information Function (TIF)
#'   \item \code{which = 4}: Person Ability Distribution
#' }
#'
#' @examples
#' \dontrun{
#' # Fit 2PL model
#' fit <- fit_irt(responses, model = "2PL")
#'
#' # Plot all diagnostics for items 1-3
#' plot(fit, which = 1:4, items = 1:3)
#'
#' # Plot only ICCs
#' plot(fit, which = 1, items = 1:5)
#' }
#'
#' @export
plot.gllamm_irt <- function(x, which = 1:4, items = NULL, ...) {

  # Determine which items to plot
  if (is.null(items)) {
    items <- if (x$n_items <= 6) 1:x$n_items else 1:6
  }

  # Ensure items are valid
  items <- items[items <= x$n_items]

  # Set up plotting area if multiple plots
  if (length(which) > 1) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    n_plots <- length(which)
    if (n_plots <= 4) {
      par(mfrow = c(2, 2))
    } else {
      n_rows <- ceiling(sqrt(n_plots))
      n_cols <- ceiling(n_plots / n_rows)
      par(mfrow = c(n_rows, n_cols))
    }
  }

  # Generate each requested plot
  for (w in which) {
    if (w == 1) plot_icc_irt(x, items, ...)
    if (w == 2) plot_iif_irt(x, items, ...)
    if (w == 3) plot_tif_irt(x, ...)
    if (w == 4) plot_ability_distribution_irt(x, ...)
  }

  invisible(x)
}


#' Plot Item Characteristic Curves
#' @keywords internal
plot_icc_irt <- function(x, items, ...) {

  theta_seq <- seq(-4, 4, length.out = 200)

  # For dichotomous models
  if (x$model %in% c("Rasch", "2PL", "3PL")) {

    # Set up colors
    colors <- rainbow(length(items))

    # Initialize plot
    plot(theta_seq, rep(0, length(theta_seq)), type = "n",
         xlab = expression(theta ~ "(Ability)"),
         ylab = "P(Correct)",
         main = "Item Characteristic Curves",
         ylim = c(0, 1),
         las = 1)

    # Add grid
    grid(col = "gray90")

    # Plot each item
    for (idx in seq_along(items)) {
      item <- items[idx]

      if (x$model == "Rasch") {
        # Rasch: P(theta) = exp(theta - b) / (1 + exp(theta - b))
        b <- x$item_parameters$difficulty[item]
        prob <- plogis(theta_seq - b)

      } else if (x$model == "2PL") {
        # 2PL: P(theta) = 1 / (1 + exp(-a(theta - b)))
        a <- x$item_parameters$discrimination[item]
        b <- x$item_parameters$difficulty[item]
        prob <- plogis(a * (theta_seq - b))

      } else if (x$model == "3PL") {
        # 3PL: P(theta) = c + (1-c) / (1 + exp(-a(theta - b)))
        a <- x$item_parameters$discrimination[item]
        b <- x$item_parameters$difficulty[item]
        c <- x$item_parameters$guessing[item]
        prob <- c + (1 - c) * plogis(a * (theta_seq - b))
      }

      lines(theta_seq, prob, col = colors[idx], lwd = 2)

      # Mark difficulty parameter
      if (x$model != "3PL") {
        abline(v = b, col = colors[idx], lty = 2, lwd = 1)
      }
    }

    # Add legend
    legend("bottomright",
           legend = paste("Item", items),
           col = colors,
           lwd = 2,
           bty = "n")

  } else {
    # Polytomous models (GRM, PCM, etc.)
    # For simplicity, plot first item with all categories
    item <- items[1]

    plot(theta_seq, rep(0, length(theta_seq)), type = "n",
         xlab = expression(theta ~ "(Ability)"),
         ylab = "Probability",
         main = paste("Category Response Curves - Item", item),
         ylim = c(0, 1),
         las = 1)

    grid(col = "gray90")

    # Number of categories for this item
    if (x$model == "GRM") {
      # Graded Response Model
      n_cat <- length(x$item_parameters$thresholds[[item]]) + 1
      a <- x$item_parameters$discrimination[item]
      thresholds <- x$item_parameters$thresholds[[item]]

      # Compute cumulative probabilities
      cum_probs <- matrix(0, length(theta_seq), n_cat - 1)
      for (k in 1:(n_cat - 1)) {
        cum_probs[, k] <- plogis(a * (theta_seq - thresholds[k]))
      }

      # Convert to category probabilities
      cat_probs <- matrix(0, length(theta_seq), n_cat)
      cat_probs[, 1] <- cum_probs[, 1]
      for (k in 2:(n_cat - 1)) {
        cat_probs[, k] <- cum_probs[, k] - cum_probs[, k - 1]
      }
      cat_probs[, n_cat] <- 1 - cum_probs[, n_cat - 1]

      # Plot each category
      colors <- rainbow(n_cat)
      for (k in 1:n_cat) {
        lines(theta_seq, cat_probs[, k], col = colors[k], lwd = 2)
      }

      legend("topright",
             legend = paste("Category", 1:n_cat),
             col = colors,
             lwd = 2,
             bty = "n")
    } else {
      text(0, 0.5, "Polytomous ICC plot\nnot yet implemented\nfor this model",
           cex = 1.2, col = "gray50")
    }
  }
}


#' Plot Item Information Functions
#' @keywords internal
plot_iif_irt <- function(x, items, ...) {

  theta_seq <- seq(-4, 4, length.out = 200)
  colors <- rainbow(length(items))

  # Initialize plot
  plot(theta_seq, rep(0, length(theta_seq)), type = "n",
       xlab = expression(theta ~ "(Ability)"),
       ylab = "Information",
       main = "Item Information Functions",
       ylim = c(0, NA),
       las = 1)

  grid(col = "gray90")

  # Compute and plot information for each item
  for (idx in seq_along(items)) {
    item <- items[idx]

    if (x$model == "Rasch") {
      # I(theta) = P(theta)(1 - P(theta))
      b <- x$item_parameters$difficulty[item]
      p <- plogis(theta_seq - b)
      info <- p * (1 - p)

    } else if (x$model == "2PL") {
      # I(theta) = a² P(theta)(1 - P(theta))
      a <- x$item_parameters$discrimination[item]
      b <- x$item_parameters$difficulty[item]
      p <- plogis(a * (theta_seq - b))
      info <- a^2 * p * (1 - p)

    } else if (x$model == "3PL") {
      # I(theta) = a²(1-c)² P*(1-P*) / [(1-c+cP*)²P]
      a <- x$item_parameters$discrimination[item]
      b <- x$item_parameters$difficulty[item]
      c <- x$item_parameters$guessing[item]
      p <- plogis(a * (theta_seq - b))
      p_star <- c + (1 - c) * p
      info <- a^2 * (1 - c)^2 * p * (1 - p) / p_star^2

    } else {
      # Polytomous: use numerical approximation
      info <- rep(0.5, length(theta_seq))  # Placeholder
    }

    lines(theta_seq, info, col = colors[idx], lwd = 2)
  }

  # Add legend
  legend("topright",
         legend = paste("Item", items),
         col = colors,
         lwd = 2,
         bty = "n")
}


#' Plot Test Information Function
#' @keywords internal
plot_tif_irt <- function(x, ...) {

  theta_seq <- seq(-4, 4, length.out = 200)

  # Compute test information (sum of all item information)
  test_info <- rep(0, length(theta_seq))

  for (item in 1:x$n_items) {

    if (x$model == "Rasch") {
      b <- x$item_parameters$difficulty[item]
      p <- plogis(theta_seq - b)
      info <- p * (1 - p)

    } else if (x$model == "2PL") {
      a <- x$item_parameters$discrimination[item]
      b <- x$item_parameters$difficulty[item]
      p <- plogis(a * (theta_seq - b))
      info <- a^2 * p * (1 - p)

    } else if (x$model == "3PL") {
      a <- x$item_parameters$discrimination[item]
      b <- x$item_parameters$difficulty[item]
      c <- x$item_parameters$guessing[item]
      p <- plogis(a * (theta_seq - b))
      p_star <- c + (1 - c) * p
      info <- a^2 * (1 - c)^2 * p * (1 - p) / p_star^2

    } else {
      info <- rep(0, length(theta_seq))  # Placeholder for polytomous
    }

    test_info <- test_info + info
  }

  # Two-panel plot: Information and SE
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2))

  # Panel 1: Test Information
  plot(theta_seq, test_info, type = "l", lwd = 2, col = "blue",
       xlab = expression(theta ~ "(Ability)"),
       ylab = "Information",
       main = "Test Information Function",
       las = 1)
  grid(col = "gray90")

  # Mark maximum
  max_idx <- which.max(test_info)
  points(theta_seq[max_idx], test_info[max_idx], pch = 19, col = "red", cex = 1.5)
  text(theta_seq[max_idx], test_info[max_idx],
       paste0("\nMax = ", round(test_info[max_idx], 2)),
       pos = 3, col = "red")

  # Panel 2: Standard Error (SE = 1/sqrt(I))
  se <- 1 / sqrt(test_info + 1e-10)
  plot(theta_seq, se, type = "l", lwd = 2, col = "red",
       xlab = expression(theta ~ "(Ability)"),
       ylab = "Standard Error",
       main = "Measurement Error",
       las = 1)
  grid(col = "gray90")

  # Restore par
  par(old_par)
}


#' Plot Person Ability Distribution
#' @keywords internal
plot_ability_distribution_irt <- function(x, ...) {

  # Extract person abilities
  if (is.null(x$person_abilities)) {
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = "Person Ability Distribution")
    text(0, 0, "Person abilities not available", cex = 1.2, col = "gray50")
    return(invisible(NULL))
  }

  abilities <- x$person_abilities

  # Histogram with normal overlay
  hist(abilities,
       breaks = 30,
       col = "lightblue",
       border = "white",
       main = "Distribution of Person Abilities",
       xlab = expression(theta ~ "(Ability)"),
       ylab = "Frequency",
       probability = TRUE,
       las = 1)

  # Add normal density curve
  curve(dnorm(x, mean(abilities), sd(abilities)),
        add = TRUE, col = "darkred", lwd = 2)

  # Add rug plot
  rug(abilities, col = "blue", lwd = 0.5)

  # Add summary statistics
  legend("topright",
         legend = c(
           paste("Mean:", round(mean(abilities), 2)),
           paste("SD:", round(sd(abilities), 2)),
           paste("N:", length(abilities))
         ),
         bty = "n")
}
