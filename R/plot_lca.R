#' Plotting Functions for Latent Class Analysis
#'
#' Visualize latent class analysis results
#'
#' @name plot_lca
NULL


#' Plot Latent Class Analysis Results
#'
#' Create diagnostic plots for latent class models including class profiles,
#' item probability heatmaps, and classification summaries
#'
#' @param x A gllamm_lca object
#' @param which Which plots to produce (1=profiles, 2=heatmap, 3=classification)
#' @param ... Additional arguments passed to plotting functions
#'
#' @details
#' Plot types:
#' \itemize{
#'   \item \code{which = 1}: Class Profiles (line plot of item probabilities by class)
#'   \item \code{which = 2}: Item Probability Heatmap
#'   \item \code{which = 3}: Classification Summary (barplot of class assignments)
#' }
#'
#' @examples
#' \dontrun{
#' # Fit LCA model
#' fit <- fit_lca(indicators, nclass = 3)
#'
#' # Plot all diagnostics
#' plot(fit, which = 1:3)
#'
#' # Plot only class profiles
#' plot(fit, which = 1)
#' }
#'
#' @export
plot.gllamm_lca <- function(x, which = 1:3, ...) {

  # Set up plotting area if multiple plots
  if (length(which) > 1) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)

    if (length(which) == 2) {
      par(mfrow = c(1, 2))
    } else if (length(which) == 3) {
      par(mfrow = c(2, 2))
    }
  }

  # Generate each requested plot
  for (w in which) {
    if (w == 1) plot_class_profiles_lca(x, ...)
    if (w == 2) plot_item_probabilities_lca(x, ...)
    if (w == 3) plot_classification_lca(x, ...)
  }

  invisible(x)
}


#' Plot Class Profiles
#' @keywords internal
plot_class_profiles_lca <- function(x, ...) {

  # Item probabilities: stored as n_items x n_classes; transpose to
  # n_classes x n_items for plotting
  item_probs <- t(x$item_probs)
  n_classes <- x$nclass

  # Colors for each class
  colors <- rainbow(n_classes)

  # Plot setup
  n_items <- ncol(item_probs)
  plot(1:n_items, rep(0, n_items), type = "n",
       xlab = "Item",
       ylab = "P(Response = 1)",
       main = "Latent Class Profiles",
       ylim = c(0, 1),
       las = 1,
       xaxt = "n")

  # Add x-axis with item labels
  axis(1, at = 1:n_items, labels = 1:n_items)

  # Add grid
  grid(col = "gray90")

  # Add horizontal reference lines
  abline(h = c(0.25, 0.5, 0.75), lty = 2, col = "gray70")

  # Plot each class as a line
  for (k in 1:n_classes) {
    lines(1:n_items, item_probs[k, ],
          type = "b", pch = 19, col = colors[k], lwd = 2)
  }

  # Add legend
  legend("topright",
         legend = paste("Class", 1:n_classes),
         col = colors,
         lwd = 2,
         pch = 19,
         bty = "n")
}


#' Plot Item Probability Heatmap
#' @keywords internal
plot_item_probabilities_lca <- function(x, ...) {

  # Item probabilities: stored as n_items x n_classes; transpose to
  # n_classes x n_items so downstream code sees the expected orientation
  item_probs <- t(x$item_probs)

  # Transpose for better visualization (items on x-axis)
  item_probs_t <- t(item_probs)

  # Create color palette
  n_colors <- 50
  color_palette <- colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(n_colors)

  # Create heatmap using image()
  n_items <- ncol(item_probs)
  n_classes <- nrow(item_probs)

  # Set up plot
  image(1:n_items, 1:n_classes, item_probs_t,
        col = color_palette,
        xlab = "Item",
        ylab = "Class",
        main = "Item Response Probabilities",
        las = 1,
        axes = FALSE)

  # Add axes
  axis(1, at = 1:n_items, labels = 1:n_items)
  axis(2, at = 1:n_classes, labels = paste("Class", 1:n_classes), las = 1)
  box()

  # Add grid
  for (i in 0:n_items) {
    abline(v = i + 0.5, col = "gray80", lwd = 0.5)
  }
  for (k in 0:n_classes) {
    abline(h = k + 0.5, col = "gray80", lwd = 0.5)
  }

  # Add text values
  for (k in 1:n_classes) {
    for (i in 1:n_items) {
      prob_val <- item_probs[k, i]
      # Choose text color based on background
      text_col <- if (prob_val > 0.6) "white" else "black"
      text(i, k, sprintf("%.2f", prob_val),
           col = text_col, cex = 0.8)
    }
  }

  # Add color legend
  legend_x <- par("usr")[2] * 1.05
  legend_y_vals <- seq(par("usr")[3], par("usr")[4], length.out = 10)

  for (i in 1:9) {
    rect(legend_x, legend_y_vals[i], legend_x + 0.3, legend_y_vals[i + 1],
         col = color_palette[floor(i * n_colors / 10)],
         border = NA, xpd = TRUE)
  }
}


#' Plot Classification Summary
#' @keywords internal
plot_classification_lca <- function(x, ...) {

  # Get posterior probabilities
  posterior <- x$posterior
  n_classes <- x$nclass

  # Assign each person to most likely class
  modal_class <- apply(posterior, 1, which.max)

  # Count assignments
  class_counts <- table(factor(modal_class, levels = 1:n_classes))

  # Create barplot
  colors <- rainbow(n_classes)

  bp <- barplot(class_counts,
                main = "Class Assignments (Modal)",
                xlab = "Latent Class",
                ylab = "Frequency",
                col = colors,
                border = "white",
                las = 1,
                names.arg = paste("Class", 1:n_classes))

  # Add counts on top of bars
  text(bp, class_counts, labels = class_counts,
       pos = 3, cex = 1.1, font = 2)

  # Add proportions as text
  class_props <- class_counts / sum(class_counts)
  text(bp, class_counts / 2,
       labels = sprintf("(%.1f%%)", class_props * 100),
       cex = 0.9, col = "white", font = 2)

  # Add classification quality info
  if (!is.null(posterior)) {
    # Compute entropy
    entropy_raw <- -sum(posterior * log(posterior + 1e-10), na.rm = TRUE)
    max_entropy <- nrow(posterior) * log(n_classes)
    entropy <- 1 - entropy_raw / max_entropy

    # Add text box with entropy
    legend("topright",
           legend = c(
             paste("N =", nrow(posterior)),
             paste("Entropy =", round(entropy, 3)),
             if (entropy > 0.8) "(Excellent)" else
               if (entropy > 0.6) "(Good)" else
                 if (entropy > 0.4) "(Fair)" else "(Poor)"
           ),
           bty = "n",
           cex = 0.9)
  }
}


#' Plot Individual Classification Uncertainty
#'
#' Visualize posterior probabilities for individual cases
#'
#' @param x A gllamm_lca object
#' @param cases Which cases to plot (default: first 20)
#' @param sort_by Sort cases by: "entropy" (default), "modal", or "index"
#' @param ... Additional arguments
#'
#' @examples
#' \dontrun{
#' fit <- fit_lca(data, nclass = 3)
#' plot_classification_uncertainty(fit, cases = 1:30)
#' }
#'
#' @export
plot_classification_uncertainty <- function(x, cases = 1:min(20, nrow(x$posterior)),
                                           sort_by = c("entropy", "modal", "index"), ...) {

  sort_by <- match.arg(sort_by)
  posterior <- x$posterior

  # Ensure cases are valid
  cases <- cases[cases <= nrow(posterior)]

  # Sort cases if requested
  if (sort_by == "entropy") {
    # Sort by individual entropy (higher = more uncertain)
    ind_entropy <- -rowSums(posterior * log(posterior + 1e-10))
    cases <- order(ind_entropy, decreasing = TRUE)[1:length(cases)]
  } else if (sort_by == "modal") {
    # Sort by modal class
    modal_class <- apply(posterior, 1, which.max)
    cases <- cases[order(modal_class[cases])]
  }

  # Subset posterior
  post_subset <- posterior[cases, , drop = FALSE]

  # Create stacked barplot
  colors <- rainbow(x$nclass)

  barplot(t(post_subset),
          main = "Classification Probabilities",
          xlab = "Case",
          ylab = "Posterior Probability",
          col = colors,
          border = "white",
          las = 1,
          names.arg = cases,
          cex.names = 0.7)

  # Add legend
  legend("topright",
         legend = paste("Class", 1:x$nclass),
         fill = colors,
         bty = "n")

  # Add reference line at 0.5
  abline(h = 0.5, lty = 2, col = "gray50")
}
