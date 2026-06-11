#' Test for Differential Item Functioning (DIF)
#'
#' Logistic-regression DIF (Swaminathan & Rogers 1990; Zumbo 1999) with
#' latent-trait or observed-score matching, generalized to multiple DIF
#' variables and their interactions through a formula interface, with
#' iterative purification of the matching criterion. Polytomous items are
#' tested with cumulative-logit (proportional odds) regression.
#'
#' @param response_matrix Item response matrix (persons x items): binary
#'   0/1 or polytomous 1..K (NA allowed)
#' @param dif The DIF specification: either a single grouping vector, or a
#'   one-sided formula over columns of \code{person_data} - e.g.
#'   \code{~ gender}, \code{~ gender + language}, or
#'   \code{~ gender * language} (the interaction tests whether DIF for one
#'   factor differs by the level of the other).
#' @param person_data Data frame with the DIF variables (required when
#'   \code{dif} is a formula)
#' @param model Matching (measurement) model for the latent criterion:
#'   "auto" (default: 2PL for dichotomous, GRM for polytomous), or any of
#'   "Rasch", "2PL", "GRM", "PCM", "GPCM"
#' @param match "theta" (default; EAP score from the anchor items under
#'   \code{model}) or "score" (observed anchor-set total score, the
#'   classical Swaminathan-Rogers criterion, comparable to
#'   \code{difR::difLogistic})
#' @param items Item indices to test (default: all)
#' @param type "both" (default; joint test of uniform + nonuniform DIF),
#'   "uniform" (group effects given the criterion), or "nonuniform"
#'   (criterion x group interactions given the uniform terms)
#' @param purify Iteratively purify the matching criterion (default TRUE):
#'   re-derive it from the currently unflagged items and re-test, until
#'   the flagged set stabilizes
#' @param anchors Optional item indices guaranteed DIF-free; they are
#'   always part of the matching criterion and never tested
#' @param alpha Significance level for flagging (default 0.05)
#' @param p_adjust Multiple-testing correction passed to
#'   \code{\link[stats]{p.adjust}} ("none" default; e.g. "BH", "holm")
#' @param max_iter Maximum purification iterations (default 10)
#'
#' @return An object of class \code{dif_analysis}: \code{dif_results}
#'   (per item: LR chi-square, df, p, adjusted p, Nagelkerke
#'   \eqn{\Delta R^2} effect size with the Jodoin-Gierl A/B/C
#'   classification, flag), \code{flagged_items}, \code{anchor_items},
#'   \code{purification} (iterations, history, converged), the matching
#'   scores, and per-item full-model coefficients for plotting.
#'
#' @details
#' For each studied item the nested models
#' \deqn{M_0: y ~ m,\quad M_1: y ~ m + Z,\quad M_2: y ~ m + Z + m:Z}
#' are fitted, where m is the matching criterion and Z the design matrix
#' of the DIF formula. Uniform DIF is \eqn{M_1} vs \eqn{M_0}, nonuniform
#' \eqn{M_2} vs \eqn{M_1}, and "both" \eqn{M_2} vs \eqn{M_0}. With
#' multiple DIF variables each test has as many degrees of freedom as Z
#' has columns, and an interaction in the formula (e.g.
#' \code{~ g1 * g2}) tests intersectional DIF beyond the additive
#' effects. Effect sizes are Nagelkerke \eqn{\Delta R^2} between the
#' compared models (A < 0.035, B < 0.07, C otherwise; Jodoin & Gierl
#' 2001).
#'
#' Purification (Lord 1980; Candell & Drasgow 1988): items flagged in one
#' round are removed from the matching criterion for the next, so DIF
#' items do not contaminate the score against which DIF is judged.
#'
#' @examples
#' \dontrun{
#' # Single factor, purified
#' res <- dif_test(resp, dif = gender)
#'
#' # Two factors plus their interaction, latent matching
#' res <- dif_test(resp, dif = ~ gender * language, person_data = persons)
#' summary(res)
#' dif_plot(res, item = 3, by = "gender")
#' }
#'
#' @export
dif_test <- function(response_matrix, dif, person_data = NULL,
                     model = c("auto", "Rasch", "2PL", "GRM", "PCM",
                               "GPCM"),
                     match = c("theta", "score"),
                     items = NULL,
                     type = c("both", "uniform", "nonuniform"),
                     purify = TRUE, anchors = NULL,
                     alpha = 0.05, p_adjust = "none",
                     max_iter = 10) {
  model <- match.arg(model)
  match <- match.arg(match)
  type <- match.arg(type)
  response_matrix <- as.matrix(response_matrix)
  n_persons <- nrow(response_matrix)
  n_items <- ncol(response_matrix)
  item_names <- colnames(response_matrix) %||% paste0("Item", 1:n_items)

  # ---- DIF design matrix from formula or vector ----
  if (inherits(dif, "formula")) {
    if (is.null(person_data)) {
      stop("person_data is required when dif is a formula")
    }
    if (nrow(person_data) != n_persons) {
      stop("person_data must have one row per person (",
           n_persons, "); found ", nrow(person_data))
    }
    mm <- model.matrix(dif, data = person_data)
    dif_formula <- dif
  } else {
    if (length(dif) != n_persons) {
      stop("dif length (", length(dif), ") must match number of persons (",
           n_persons, ")")
    }
    person_data <- data.frame(group = factor(dif))
    mm <- model.matrix(~ group, data = person_data)
    dif_formula <- ~ group
  }
  has_int <- "(Intercept)" %in% colnames(mm)
  Z <- mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
  if (ncol(Z) == 0) stop("The DIF specification has no terms to test")

  # ---- Item type bookkeeping ----
  n_cats <- apply(response_matrix, 2, function(v) {
    length(unique(v[!is.na(v)]))
  })
  is_poly <- any(n_cats > 2)
  if (model == "auto") model <- if (is_poly) "GRM" else "2PL"

  if (is.null(items)) items <- seq_len(n_items)
  if (any(items < 1 | items > n_items)) {
    stop("items must be between 1 and ", n_items)
  }
  if (!is.null(anchors)) {
    if (any(anchors < 1 | anchors > n_items)) {
      stop("anchors must be valid item indices")
    }
    items <- setdiff(items, anchors)
  }

  # ---- Matching criterion from a given anchor set ----
  matching <- function(anchor_set, studied = NULL) {
    if (match == "theta") {
      f <- fit_irt(response_matrix[, anchor_set, drop = FALSE],
                   model = model)
      unname(f$person_abilities)
    } else {
      # Observed score over the anchors plus the studied item (the
      # difR::difLogistic convention)
      cols <- union(anchor_set, studied)
      rowSums(response_matrix[, cols, drop = FALSE], na.rm = TRUE)
    }
  }

  # ---- One purification round: test every studied item ----
  test_round <- function(anchor_set) {
    th_global <- if (match == "theta") matching(anchor_set) else NULL
    rows <- lapply(items, function(j) {
      m <- if (match == "theta") th_global else matching(anchor_set, j)
      .dif_item_tests(response_matrix[, j], m, Z, type, n_cats[j])
    })
    res <- do.call(rbind, rows)
    res$item <- items
    res$name <- item_names[items]
    res$p_adj <- stats::p.adjust(res$p_value, method = p_adjust)
    res$flagged <- !is.na(res$p_adj) & res$p_adj < alpha
    list(res = res, th = th_global)
  }

  # ---- Purification loop ----
  all_idx <- seq_len(n_items)
  base_anchor <- if (is.null(anchors)) all_idx else
    union(anchors, setdiff(all_idx, items))
  if (is_poly && !requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for polytomous DIF tests")
  }

  flagged <- integer(0)
  history <- list()
  converged <- TRUE
  iter <- 1
  res <- NULL
  th_final <- NULL
  last_anchor <- base_anchor

  repeat {
    anchor_set <- setdiff(base_anchor, flagged)
    if (length(anchor_set) < 2) {
      if (is.null(res)) {
        stop("Fewer than 2 anchor items available; provide more items ",
             "or explicit anchors")
      }
      # Purification breakdown (typically: too large a fraction of the
      # test has same-direction DIF for the criterion to disentangle DIF
      # from impact). Keep the last valid round's results.
      warning("Purification removed (almost) all items from the anchor; ",
              "the DIF/impact decomposition is not identified from these ",
              "items alone. Reporting the last valid round - consider ",
              "supplying known DIF-free 'anchors'.")
      converged <- FALSE
      anchor_set <- last_anchor
      break
    }
    rt <- test_round(anchor_set)
    res <- rt$res
    th_final <- rt$th
    last_anchor <- anchor_set
    new_flagged <- res$item[res$flagged]
    history[[iter]] <- new_flagged
    if (!purify || setequal(new_flagged, flagged)) {
      flagged <- new_flagged
      break
    }
    flagged <- new_flagged
    if (iter >= max_iter) {
      warning("Purification did not stabilize in ", max_iter,
              " iterations")
      converged <- FALSE
      break
    }
    iter <- iter + 1
  }

  # ---- Per-item full-model coefficients for plotting ----
  # (matching criterion from the last valid round)
  item_models <- lapply(items, function(j) {
    m <- if (match == "theta") th_final else matching(anchor_set, j)
    .dif_item_fullfit(response_matrix[, j], m, Z, n_cats[j])
  })
  names(item_models) <- item_names[items]

  result <- list(
    dif_results = res[, c("item", "name", "chisq", "df", "p_value",
                          "p_adj", "delta_R2", "classification",
                          "flagged")],
    flagged_items = flagged,
    anchor_items = anchor_set,
    matching = if (match == "theta") th_final else NULL,
    match = match,
    dif_formula = dif_formula,
    dif_terms = colnames(Z),
    person_data = person_data,
    item_models = item_models,
    purification = list(purify = purify, iterations = iter,
                        history = history, converged = converged),
    model = model,
    type = type,
    alpha = alpha,
    p_adjust = p_adjust,
    n_items = n_items
  )
  class(result) <- "dif_analysis"
  result
}


#' Nested-model LR tests for one item
#' @keywords internal
.dif_item_tests <- function(y, m, Z, type, K) {
  ok <- !is.na(y) & !is.na(m) & stats::complete.cases(Z)
  d <- data.frame(y = y[ok], m = m[ok], Z[ok, , drop = FALSE],
                  check.names = FALSE)
  if (K > 2) d$y <- factor(d$y, ordered = TRUE)
  zn <- colnames(Z)
  bt <- function(v) paste0("`", v, "`")
  f0 <- stats::as.formula("y ~ m")
  f1 <- stats::as.formula(paste("y ~ m +", paste(bt(zn), collapse = "+")))
  f2 <- stats::as.formula(paste("y ~ m +", paste(bt(zn), collapse = "+"),
                                "+", paste(paste0("m:", bt(zn)),
                                           collapse = "+")))
  out <- data.frame(chisq = NA_real_, df = NA_integer_,
                    p_value = NA_real_, delta_R2 = NA_real_,
                    classification = NA_character_,
                    stringsAsFactors = FALSE)
  fit <- function(f) {
    if (K > 2) {
      g <- tryCatch(MASS::polr(f, data = d, Hess = FALSE),
                    error = function(e) NULL)
      if (is.null(g)) return(NULL)
      list(dev = g$deviance, np = length(g$coefficients) + length(g$zeta))
    } else {
      g <- tryCatch(stats::glm(f, data = d, family = stats::binomial()),
                    error = function(e) NULL)
      if (is.null(g)) return(NULL)
      list(dev = g$deviance, np = length(stats::coef(g)))
    }
  }
  m0 <- fit(f0); m1 <- fit(f1); m2 <- fit(f2)
  null_dev <- if (K > 2) {
    g <- tryCatch(MASS::polr(y ~ 1, data = d, Hess = FALSE),
                  error = function(e) NULL)
    if (is.null(g)) NA_real_ else g$deviance
  } else {
    g <- stats::glm(y ~ 1, data = d, family = stats::binomial())
    g$deviance
  }
  if (is.null(m0) || is.null(m1) || is.null(m2) || !is.finite(null_dev)) {
    return(out)
  }
  pair <- switch(type,
                 uniform = list(m0, m1),
                 nonuniform = list(m1, m2),
                 both = list(m0, m2))
  out$chisq <- max(pair[[1]]$dev - pair[[2]]$dev, 0)
  out$df <- pair[[2]]$np - pair[[1]]$np
  out$p_value <- stats::pchisq(out$chisq, out$df, lower.tail = FALSE)

  n <- nrow(d)
  nagelkerke <- function(dev) {
    r2_cs <- 1 - exp((dev - null_dev) / n)
    r2_cs / (1 - exp(-null_dev / n))
  }
  out$delta_R2 <- nagelkerke(pair[[2]]$dev) - nagelkerke(pair[[1]]$dev)
  out$classification <- if (!is.finite(out$delta_R2)) NA_character_
    else if (out$delta_R2 < 0.035) "A"
    else if (out$delta_R2 < 0.07) "B"
    else "C"
  out
}


#' Full-model (M2) fit for one item, kept for plotting
#' @keywords internal
.dif_item_fullfit <- function(y, m, Z, K) {
  ok <- !is.na(y) & !is.na(m) & stats::complete.cases(Z)
  d <- data.frame(y = y[ok], m = m[ok], Z[ok, , drop = FALSE],
                  check.names = FALSE)
  zn <- paste0("`", colnames(Z), "`")
  f2 <- stats::as.formula(paste("y ~ m +", paste(zn, collapse = "+"),
                                "+", paste(paste0("m:", zn),
                                           collapse = "+")))
  if (K > 2) {
    d$y <- factor(d$y, ordered = TRUE)
    tryCatch(MASS::polr(f2, data = d, Hess = FALSE),
             error = function(e) NULL)
  } else {
    tryCatch(stats::glm(f2, data = d, family = stats::binomial()),
             error = function(e) NULL)
  }
}


#' Test for DIF with explicit response data (deprecated)
#'
#' Deprecated single-factor wrapper kept for backward compatibility; use
#' \code{\link{dif_test}}, which supports multiple DIF variables,
#' interactions, and iterative purification.
#'
#' @param response_matrix Matrix of item responses (persons x items)
#' @param group Grouping vector
#' @param model Matching model
#' @param items Items to test (default: all)
#' @param type Type of DIF
#' @param method Ignored (kept for compatibility)
#' @param alpha Significance level
#'
#' @return Object of class \code{dif_analysis}
#' @export
dif_test_with_data <- function(response_matrix, group,
                               model = c("Rasch", "2PL", "3PL", "GRM",
                                         "PCM", "GPCM"),
                               items = NULL,
                               type = c("both", "uniform", "nonuniform"),
                               method = NULL, alpha = 0.05) {
  .Deprecated("dif_test")
  model <- match.arg(model)
  if (model == "3PL") model <- "2PL"
  dif_test(response_matrix, dif = group, model = model, items = items,
           type = match.arg(type), alpha = alpha, purify = TRUE)
}


#' @export
print.dif_analysis <- function(x, ...) {
  cat("Differential Item Functioning (DIF) Analysis\n")
  cat("==============================================\n\n")
  cat("DIF specification:", deparse(x$dif_formula),
      "(", length(x$dif_terms), "term(s) )\n")
  cat("Matching:", if (x$match == "theta")
    paste0("latent (", x$model, " EAP)") else "observed score",
    "| Type:", x$type, "\n")
  if (isTRUE(x$purification$purify)) {
    cat("Purification:", x$purification$iterations, "iteration(s),",
        if (x$purification$converged) "converged" else "NOT converged",
        "\n")
  }
  cat("Items tested:", nrow(x$dif_results),
      "| flagged:", length(x$flagged_items),
      "(alpha =", x$alpha,
      if (x$p_adjust != "none") paste0(", ", x$p_adjust, "-adjusted")
      else "", ")\n\n")
  if (length(x$flagged_items) > 0) {
    cat("Flagged items:\n")
    fr <- x$dif_results[x$dif_results$flagged, ]
    fr[, c("chisq", "p_value", "p_adj", "delta_R2")] <-
      round(fr[, c("chisq", "p_value", "p_adj", "delta_R2")], 4)
    print(fr, row.names = FALSE)
  } else {
    cat("No items flagged for DIF\n")
  }
  invisible(x)
}


#' @export
summary.dif_analysis <- function(object, ...) {
  print(object)
  cat("\nAll items:\n")
  fr <- object$dif_results
  fr[, c("chisq", "p_value", "p_adj", "delta_R2")] <-
    round(fr[, c("chisq", "p_value", "p_adj", "delta_R2")], 4)
  print(fr, row.names = FALSE)
  cat("\nEffect sizes: Nagelkerke delta-R2 between the compared models;\n")
  cat("A < 0.035 (negligible), B < 0.07 (moderate), C (large)\n")
  invisible(object)
}


#' Plot item response curves by DIF group
#'
#' Plots the model-implied probability of a correct/positive response (or
#' the expected score, for polytomous items) against the matching
#' criterion, for each level of one DIF variable, from the per-item
#' full DIF model. Other DIF variables are held at their reference level.
#'
#' @param dif_result Object from \code{\link{dif_test}}
#' @param item Item index (position among the tested items)
#' @param by Name of the DIF variable to display (default: the first one)
#' @param ... Additional graphical parameters
#'
#' @export
dif_plot <- function(dif_result, item, by = NULL, ...) {
  if (!inherits(dif_result, "dif_analysis")) {
    stop("dif_result must be output from dif_test()")
  }
  ridx <- match(item, dif_result$dif_results$item)
  if (is.na(ridx)) stop("item ", item, " was not among the tested items")
  fit <- dif_result$item_models[[ridx]]
  if (is.null(fit)) stop("no stored model for item ", item)

  pd <- dif_result$person_data
  vars <- all.vars(dif_result$dif_formula)
  if (is.null(by)) by <- vars[1]
  if (!by %in% vars) stop("'by' must be one of: ", paste(vars, collapse = ", "))

  by_levels <- if (is.factor(pd[[by]])) levels(pd[[by]])
               else sort(unique(pd[[by]]))
  if (length(by_levels) > 8) {
    by_levels <- stats::quantile(pd[[by]], c(0.1, 0.5, 0.9))
  }

  m_seq <- if (!is.null(dif_result$matching)) {
    seq(min(dif_result$matching, na.rm = TRUE),
        max(dif_result$matching, na.rm = TRUE), length.out = 100)
  } else seq(0, dif_result$n_items, length.out = 100)

  ref_row <- pd[1, vars, drop = FALSE]
  for (v in vars) {
    ref_row[[v]] <- if (is.factor(pd[[v]])) factor(levels(pd[[v]])[1],
                                                   levels = levels(pd[[v]]))
                    else 0
  }

  curves <- lapply(by_levels, function(lv) {
    nd_vars <- ref_row[rep(1, length(m_seq)), , drop = FALSE]
    nd_vars[[by]] <- if (is.factor(pd[[by]]))
      factor(lv, levels = levels(pd[[by]])) else lv
    Znew <- model.matrix(dif_result$dif_formula, data = nd_vars)
    Znew <- Znew[, setdiff(colnames(Znew), "(Intercept)"), drop = FALSE]
    nd <- data.frame(m = m_seq, Znew, check.names = FALSE)
    if (inherits(fit, "polr")) {
      pr <- stats::predict(fit, newdata = nd, type = "probs")
      as.numeric(pr %*% seq_len(ncol(pr)))      # expected score
    } else {
      stats::predict(fit, newdata = nd, type = "response")
    }
  })

  ylim <- range(unlist(curves))
  cols <- grDevices::hcl.colors(length(by_levels), "Dark 2")
  plot(m_seq, curves[[1]], type = "l", lwd = 2, col = cols[1],
       xlab = if (dif_result$match == "theta") "Ability (theta)"
              else "Matching score",
       ylab = if (inherits(fit, "polr")) "Expected score" else "P(Y = 1)",
       main = paste0(dif_result$dif_results$name[ridx],
                     " - response curves by ", by),
       ylim = ylim, ...)
  if (length(by_levels) > 1) {
    for (k in 2:length(by_levels)) {
      graphics::lines(m_seq, curves[[k]], lwd = 2, col = cols[k])
    }
  }
  graphics::legend("bottomright", legend = as.character(by_levels),
                   col = cols, lwd = 2, title = by)
  graphics::grid()
  invisible(dif_result)
}
