#' Latent structure comparison: categorization, ordering or quantification
#'
#' Fits the six latent structure models of Torres Irribarra & Diakow's
#' model selection framework to a binary response matrix and returns a
#' comparison table. The models progressively constrain the latent
#' structure, decomposing the order and scale assumptions:
#'
#' \describe{
#'   \item{UN}{Unconstrained latent class model - qualitative structure
#'     (differences of kind).}
#'   \item{MON}{Ordered classes with class (person) monotonicity
#'     (Croon 1990).}
#'   \item{IIO}{Ordered classes with invariant item ordering (item
#'     monotonicity).}
#'   \item{DM}{Double monotonicity: both restrictions.}
#'   \item{LCR}{Located latent classes (latent class Rasch; Lindsay,
#'     Clogg & Grego 1991) - a discrete quantitative structure.}
#'   \item{RM}{The Rasch model with a normal latent distribution - a
#'     continuous quantitative structure.}
#' }
#'
#' Successive comparisons carry the framework's logic: UN vs MON/IIO asks
#' whether an ordering is tenable at all; the single-monotonicity models
#' vs DM asks whether persons and items share one proficiency
#' progression; DM vs LCR isolates the interval-scale (parameter
#' separability) assumption; LCR vs RM asks whether a continuous latent
#' variable adds anything beyond located classes. UN, MON, IIO, and DM
#' share the same nominal parameter count (inequality constraints do not
#' reduce it), so their information criteria differ only through fit;
#' treat those comparisons descriptively (chi-bar-square caveat).
#'
#' @param Y Binary response matrix (persons x items)
#' @param nclass Number of latent classes for the class-based models.
#'   The default \code{ceiling((ncol(Y) + 1) / 2)} is the smallest number
#'   at which the located class model is fit-equivalent to the
#'   semiparametric Rasch model (Lindsay et al. 1991).
#' @param item_order How to order items for the IIO/DM constraints:
#'   "auto" (default; by marginal proportion correct) or "columns" (the
#'   column order of \code{Y})
#' @param n_starts Random starts for each class-based model (default 5)
#' @param control Control list passed to the fitters
#'
#' @return An object of class \code{lca_structure_comparison}: a data
#'   frame with one row per model (structure type, logLik, nominal
#'   parameter count, AIC, BIC) plus the fitted models in
#'   \code{attr(, "fits")}.
#'
#' @examples
#' \dontrun{
#' cmp <- latent_structure_comparison(resp, nclass = 4)
#' print(cmp)
#' }
#'
#' @references
#' Torres Irribarra, D., & Diakow, R. Categorization, ordering or
#' quantification: selecting a latent variable model by comparing latent
#' structures.
#'
#' @export
latent_structure_comparison <- function(Y, nclass = NULL,
                                        item_order = c("auto", "columns"),
                                        n_starts = 5, control = list()) {
  item_order <- match.arg(item_order)
  Y <- as.matrix(Y)
  vals <- Y[!is.na(Y)]
  if (!all(vals %in% c(0, 1))) {
    stop("latent_structure_comparison requires binary (0/1) responses")
  }
  n_items <- ncol(Y)
  if (is.null(nclass)) nclass <- ceiling((n_items + 1) / 2)

  # Item order for the IIO/DM constraints: easiest-to-hardest chain
  ord <- if (item_order == "auto") {
    order(colMeans(Y, na.rm = TRUE))
  } else {
    seq_len(n_items)
  }
  item_pairs <- lapply(seq_len(n_items - 1),
                       function(k) c(ord[k], ord[k + 1]))

  ctl <- modifyList(list(n_starts = n_starts), control)

  fits <- list(
    UN = fit_lca(Y, nclass = nclass, control = ctl),
    MON = fit_lca(Y, nclass = nclass, ordering = "increasing",
                  control = ctl),
    IIO = fit_lca(Y, nclass = nclass, item_ordering = item_pairs,
                  control = ctl),
    DM = fit_lca(Y, nclass = nclass, ordering = "increasing",
                 item_ordering = item_pairs, control = ctl),
    LCR = fit_lca(Y, nclass = nclass, structure = "rasch",
                  control = ctl),
    RM = fit_irt(Y, model = "Rasch", se = FALSE)
  )

  # Delegate the comparison table to the general machinery, then add the
  # framework's structure labels
  cmp <- compare_models(UN = fits$UN, MON = fits$MON, IIO = fits$IIO,
                        DM = fits$DM, LCR = fits$LCR, RM = fits$RM)
  tab <- as.data.frame(cmp)
  tab$structure <- c("qualitative", "ordinal", "ordinal", "ordinal",
                     "quantitative (discrete)",
                     "quantitative (continuous)")
  tab <- tab[, c("model", "structure", "logLik", "n_params", "AIC",
                 "dAIC", "BIC", "dBIC")]
  rownames(tab) <- NULL

  structure(tab, fits = fits, nclass = nclass, item_order = ord,
            class = c("lca_structure_comparison",
                      "gllammr_model_comparison", "data.frame"))
}


#' @export
print.lca_structure_comparison <- function(x, ...) {
  cat("Latent structure comparison (Torres Irribarra & Diakow",
      "framework)\n")
  cat("Classes:", attr(x, "nclass"),
      "| IIO/DM item order:", paste(attr(x, "item_order"),
                                    collapse = " < "), "\n\n")
  tab <- as.data.frame(x)
  num <- vapply(tab, is.numeric, TRUE)
  tab[num] <- lapply(tab[num], round, 2)
  print(tab, row.names = FALSE)
  best <- tab$model[which.min(tab$BIC)]
  cat("\nLowest BIC:", best, "\n")
  cat("Read successively: UN vs MON/IIO (is ordering tenable?),\n")
  cat("single vs double monotonicity (one shared progression?),\n")
  cat("DM vs LCR (interval scale?), LCR vs RM (continuum vs grain).\n")
  invisible(x)
}
