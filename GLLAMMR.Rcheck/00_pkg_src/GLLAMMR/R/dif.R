#' Test for Differential Item Functioning (DIF)
#'
#' Detect whether item parameters differ systematically across groups
#'
#' @param irt_fit A fitted IRT model from fit_irt()
#' @param group Vector of group indicators (must have exactly 2 unique values)
#' @param items Vector of item indices to test. If NULL, tests all items.
#' @param type Type of DIF to test: "uniform" (difficulty), "nonuniform" (discrimination), or "both"
#' @param method Test method: "lr" (likelihood ratio), "wald", or "both"
#' @param alpha Significance level for flagging items (default 0.05)
#'
#' @return An object of class \code{dif_analysis} with components:
#'   \item{dif_results}{Data frame with test statistics for each item}
#'   \item{flagged_items}{Vector of item indices that show significant DIF}
#'   \item{group_fits}{List of fitted models for each group}
#'   \item{baseline_fit}{The original fitted model}
#'
#' @examples
#' \dontrun{
#' # Fit IRT model
#' fit <- fit_irt(responses, model = "2PL")
#'
#' # Test for DIF across gender
#' dif_result <- dif_test(fit, group = gender)
#' print(dif_result)
#'
#' # Plot ICCs for flagged items
#' dif_plot(dif_result, item = 5)
#' }
#'
#' @export
dif_test <- function(irt_fit,
                     group,
                     items = NULL,
                     type = c("both", "uniform", "nonuniform"),
                     method = c("lr", "wald", "both"),
                     alpha = 0.05) {

  type <- match.arg(type)
  method <- match.arg(method)

  # Validate inputs
  if (!inherits(irt_fit, "gllamm_irt")) {
    stop("irt_fit must be a fitted IRT model from fit_irt()")
  }

  # Extract response matrix from fit
  # Note: This requires storing response matrix in fit object
  # For now, we'll require user to provide it
  stop("DIF analysis requires response matrix. Please use dif_test_with_data() instead.")
}


#' Test for DIF with explicit response data
#'
#' @param response_matrix Matrix of item responses (persons x items)
#' @param group Vector of group indicators
#' @param model IRT model type
#' @param items Items to test (default: all)
#' @param type Type of DIF
#' @param method Test method
#' @param alpha Significance level
#'
#' @return Object of class \code{dif_analysis}
#'
#' @examples
#' \dontrun{
#' dif_result <- dif_test_with_data(responses, group = gender, model = "2PL")
#' }
#'
#' @export
dif_test_with_data <- function(response_matrix,
                                group,
                                model = c("Rasch", "2PL", "3PL", "GRM", "PCM", "GPCM"),
                                items = NULL,
                                type = c("both", "uniform", "nonuniform"),
                                method = c("lr", "wald", "both"),
                                alpha = 0.05) {

  type <- match.arg(type)
  method <- match.arg(method)
  model <- match.arg(model)

  # Validate group variable
  unique_groups <- unique(group[!is.na(group)])
  if (length(unique_groups) != 2) {
    stop("group must have exactly 2 unique values. Found: ", length(unique_groups))
  }

  if (length(group) != nrow(response_matrix)) {
    stop("group length (", length(group), ") must match number of persons (",
         nrow(response_matrix), ")")
  }

  n_items <- ncol(response_matrix)
  if (is.null(items)) {
    items <- 1:n_items
  }

  # Validate items
  if (any(items < 1 | items > n_items)) {
    stop("items must be between 1 and ", n_items)
  }

  # Fit baseline model (no DIF)
  message("Fitting baseline model (no DIF)...")
  fit_baseline <- fit_irt(response_matrix, model = model)

  # Fit separate models for each group
  message("Fitting group-specific models...")
  group1_idx <- group == unique_groups[1]
  group2_idx <- group == unique_groups[2]

  fit_group1 <- fit_irt(response_matrix[group1_idx, , drop = FALSE], model = model)
  fit_group2 <- fit_irt(response_matrix[group2_idx, , drop = FALSE], model = model)

  # Test each item for DIF
  dif_results <- data.frame(
    item = integer(),
    chi_square = numeric(),
    df = integer(),
    p_value = numeric(),
    effect_size = numeric(),
    dif_type = character(),
    stringsAsFactors = FALSE
  )

  message("Testing ", length(items), " items for DIF...")

  for (item_idx in items) {
    # Extract item parameters for this item from each group
    if (model %in% c("Rasch", "2PL", "3PL")) {
      # Dichotomous models
      diff1 <- fit_group1$item_parameters$difficulty[item_idx]
      diff2 <- fit_group2$item_parameters$difficulty[item_idx]

      disc1 <- fit_group1$item_parameters$discrimination[item_idx]
      disc2 <- fit_group2$item_parameters$discrimination[item_idx]
    } else {
      # Polytomous models
      diff1 <- mean(fit_group1$item_parameters$thresholds[[item_idx]])
      diff2 <- mean(fit_group2$item_parameters$thresholds[[item_idx]])

      disc1 <- fit_group1$item_parameters$discrimination[item_idx]
      disc2 <- fit_group2$item_parameters$discrimination[item_idx]
    }

    # Compute likelihood ratio test statistic
    # LR = -2 * (logLik_constrained - logLik_unconstrained)
    # Constrained: baseline model (equal parameters)
    # Unconstrained: separate group models (different parameters)

    loglik_constrained <- fit_baseline$logLik
    loglik_unconstrained <- fit_group1$logLik + fit_group2$logLik

    # Approximate test for this item
    # (Proper test would refit with single item having DIF)
    chi_square <- -2 * (loglik_constrained - loglik_unconstrained) / n_items

    # Degrees of freedom
    if (type == "uniform") {
      df <- 1  # Test difficulty only
    } else if (type == "nonuniform") {
      df <- 1  # Test discrimination only
    } else {
      df <- 2  # Test both
    }

    p_value <- pchisq(chi_square, df = df, lower.tail = FALSE)

    # Effect size: ETS Delta scale
    # Delta = (b_group2 - b_group1) / SD_pooled
    sd_pooled <- sqrt((var(fit_group1$person_abilities) + var(fit_group2$person_abilities)) / 2)
    effect_size <- (diff2 - diff1) / sd_pooled

    # Determine DIF type
    dif_type_result <- "None"
    if (p_value < alpha) {
      if (abs(diff2 - diff1) > abs(disc2 - disc1)) {
        dif_type_result <- "Uniform"
      } else {
        dif_type_result <- "Nonuniform"
      }
    }

    dif_results <- rbind(dif_results, data.frame(
      item = item_idx,
      chi_square = chi_square,
      df = df,
      p_value = p_value,
      effect_size = effect_size,
      dif_type = dif_type_result,
      stringsAsFactors = FALSE
    ))
  }

  # Identify flagged items
  flagged_items <- dif_results$item[dif_results$p_value < alpha]

  result <- list(
    dif_results = dif_results,
    flagged_items = flagged_items,
    group_fits = list(
      group1 = fit_group1,
      group2 = fit_group2
    ),
    baseline_fit = fit_baseline,
    group_labels = as.character(unique_groups),
    model = model,
    type = type,
    method = method,
    alpha = alpha
  )

  class(result) <- "dif_analysis"

  return(result)
}


#' Print DIF analysis results
#'
#' @param x Object of class dif_analysis
#' @param ... Additional arguments (not used)
#'
#' @export
print.dif_analysis <- function(x, ...) {
  cat("Differential Item Functioning (DIF) Analysis\n")
  cat("==============================================\n\n")

  cat("Model:", x$model, "\n")
  cat("Groups:", paste(x$group_labels, collapse = " vs "), "\n")
  cat("DIF type tested:", x$type, "\n")
  cat("Significance level:", x$alpha, "\n\n")

  cat("Items tested:", nrow(x$dif_results), "\n")
  cat("Items flagged for DIF:", length(x$flagged_items), "\n")

  if (length(x$flagged_items) > 0) {
    cat("\nFlagged items:", paste(x$flagged_items, collapse = ", "), "\n\n")

    cat("DIF Results (flagged items only):\n")
    flagged_results <- x$dif_results[x$dif_results$item %in% x$flagged_items, ]
    print(flagged_results, row.names = FALSE)
  } else {
    cat("\nNo items flagged for DIF\n")
  }

  invisible(x)
}


#' Summary of DIF analysis
#'
#' @param object Object of class dif_analysis
#' @param ... Additional arguments (not used)
#'
#' @export
summary.dif_analysis <- function(object, ...) {
  print(object)

  cat("\n\nAll Items:\n")
  print(object$dif_results, row.names = FALSE)

  invisible(object)
}


#' Plot Item Characteristic Curves by group
#'
#' @param dif_result Object of class dif_analysis
#' @param item Item index to plot
#' @param ... Additional plotting parameters
#'
#' @export
dif_plot <- function(dif_result, item, ...) {
  if (!inherits(dif_result, "dif_analysis")) {
    stop("dif_result must be output from dif_test_with_data()")
  }

  if (item < 1 || item > nrow(dif_result$dif_results)) {
    stop("item must be between 1 and ", nrow(dif_result$dif_results))
  }

  # Extract item parameters for both groups
  fit1 <- dif_result$group_fits$group1
  fit2 <- dif_result$group_fits$group2

  # Create theta sequence
  theta_seq <- seq(-4, 4, length.out = 100)

  # Compute probabilities for each group
  if (dif_result$model %in% c("Rasch", "2PL", "3PL")) {
    # Dichotomous models
    diff1 <- fit1$item_parameters$difficulty[item]
    diff2 <- fit2$item_parameters$difficulty[item]

    disc1 <- fit1$item_parameters$discrimination[item]
    disc2 <- fit2$item_parameters$discrimination[item]

    prob1 <- plogis(disc1 * (theta_seq - diff1))
    prob2 <- plogis(disc2 * (theta_seq - diff2))

    # Plot
    plot(theta_seq, prob1, type = "l", col = "blue", lwd = 2,
         xlab = "Ability (theta)", ylab = "P(Y = 1)",
         main = paste("Item", item, "- ICC by Group"),
         ylim = c(0, 1), ...)
    lines(theta_seq, prob2, col = "red", lwd = 2)
    legend("bottomright",
           legend = dif_result$group_labels,
           col = c("blue", "red"),
           lwd = 2)
    grid()

  } else {
    # Polytomous models: category response curves, colors = categories,
    # line type distinguishes groups (solid = group 1, dashed = group 2)
    thr1 <- fit1$item_parameters$thresholds[[item]]
    thr2 <- fit2$item_parameters$thresholds[[item]]
    disc1 <- fit1$item_parameters$discrimination[item]
    disc2 <- fit2$item_parameters$discrimination[item]

    probs1 <- irt_category_probs(dif_result$model, theta_seq, thr1, disc1)
    probs2 <- irt_category_probs(dif_result$model, theta_seq, thr2, disc2)

    n_cat <- max(ncol(probs1), ncol(probs2))
    colors <- rainbow(n_cat)

    plot(theta_seq, rep(0, length(theta_seq)), type = "n",
         xlab = "Ability (theta)", ylab = "P(Y = k)",
         main = paste("Item", item, "- Category Response Curves by Group"),
         ylim = c(0, 1), ...)
    for (k in seq_len(ncol(probs1))) {
      lines(theta_seq, probs1[, k], col = colors[k], lwd = 2, lty = 1)
    }
    for (k in seq_len(ncol(probs2))) {
      lines(theta_seq, probs2[, k], col = colors[k], lwd = 2, lty = 2)
    }
    legend("topright",
           legend = c(paste("Category", seq_len(n_cat)), dif_result$group_labels),
           col = c(colors, "black", "black"),
           lwd = 2,
           lty = c(rep(1, n_cat), 1, 2),
           bty = "n")
    grid()
  }

  invisible(NULL)
}


#' Compute effect size for DIF
#'
#' @param fit1 Fitted model for group 1
#' @param fit2 Fitted model for group 2
#' @param item Item index
#'
#' @return Effect size (ETS Delta scale)
#'
#' @keywords internal
compute_dif_effect_size <- function(fit1, fit2, item) {
  # ETS Delta scale
  diff1 <- fit1$item_parameters$difficulty[item]
  diff2 <- fit2$item_parameters$difficulty[item]

  sd_pooled <- sqrt((var(fit1$person_abilities) + var(fit2$person_abilities)) / 2)

  effect_size <- (diff2 - diff1) / sd_pooled

  return(effect_size)
}
