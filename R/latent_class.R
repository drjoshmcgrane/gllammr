#' Fit Latent Class Analysis Models
#'
#' Fit latent class models to categorical data
#'
#' @param formula Formula specifying manifest variables (can be a matrix)
#' @param data Data frame containing variables
#' @param nclass Number of latent classes
#' @param start Optional starting values
#' @param control Control parameters
#'
#' @return An object of class \code{gllamm_lca}
#'
#' @examples
#' \dontrun{
#' # Simulate 2-class data
#' set.seed(123)
#' n <- 500
#'
#' # Class 1: high probability of yes
#' class1_probs <- c(0.8, 0.7, 0.9, 0.75)
#' # Class 2: low probability of yes
#' class2_probs <- c(0.2, 0.3, 0.1, 0.25)
#'
#' # Generate data
#' true_class <- sample(1:2, n, replace = TRUE, prob = c(0.6, 0.4))
#' data <- matrix(NA, n, 4)
#' for (i in 1:n) {
#'   probs <- if (true_class[i] == 1) class1_probs else class2_probs
#'   data[i, ] <- rbinom(4, 1, probs)
#' }
#' colnames(data) <- paste0("Item", 1:4)
#'
#' # Fit 2-class model
#' fit <- fit_lca(data, nclass = 2)
#' summary(fit)
#' }
#'
#' @export
fit_lca <- function(formula, data = NULL, nclass = 2,
                    weights = NULL,
                    start = NULL, control = list()) {

  # Handle formula or matrix input
  if (is.matrix(formula)) {
    Y <- formula
  } else if (is.data.frame(formula)) {
    Y <- as.matrix(formula)
  } else {
    stop("Input must be a matrix or data frame of manifest variables")
  }

  n_obs <- nrow(Y)
  n_items <- ncol(Y)

  # Validate weights if provided
  if (!is.null(weights)) {
    if (length(weights) != n_obs) {
      stop("Length of weights (", length(weights), ") must match number of observations (", n_obs, ")")
    }
    if (any(weights < 0, na.rm = TRUE)) {
      stop("All weights must be non-negative")
    }
    if (any(is.na(weights))) {
      stop("weights cannot contain missing values")
    }
  }

  # Check binary responses
  if (!all(Y %in% c(0, 1, NA))) {
    stop("All manifest variables must be binary (0/1)")
  }

  # Remove missing (for now - could implement EM with missing data)
  complete_rows <- complete.cases(Y)
  if (sum(complete_rows) < n_obs) {
    warning("Removing ", n_obs - sum(complete_rows), " rows with missing data")
    Y <- Y[complete_rows, ]
    if (!is.null(weights)) {
      weights <- weights[complete_rows]
    }
    n_obs <- nrow(Y)
  }

  # Prepare weights vector (default to 1.0 if NULL)
  if (is.null(weights)) {
    weights_vec <- rep(1.0, n_obs)
  } else {
    weights_vec <- as.numeric(weights)
  }

  # Prepare TMB data
  tmb_data <- list(
    Y = Y,
    n_obs = as.integer(n_obs),
    n_items = as.integer(n_items),
    n_classes = as.integer(nclass),
    weights = weights_vec
  )

  # Initialize parameters
  if (is.null(start)) {
    # Initialize item probabilities with jitter around marginal proportions
    item_means <- colMeans(Y)

    item_probs_init <- matrix(NA, n_items, nclass)
    for (k in 1:nclass) {
      jitter_factor <- runif(n_items, 0.7, 1.3)
      item_probs_init[, k] <- pmin(pmax(item_means * jitter_factor, 0.05), 0.95)
    }

    # Initialize class probabilities (equal)
    class_logits_init <- rep(0, nclass - 1)

    tmb_params <- list(
      item_logits = qlogis(item_probs_init),
      class_logits = class_logits_init
    )
  } else {
    tmb_params <- start
  }

  # Create TMB object
  tmb_data$model_name <- "latent_class"
  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    DLL = "GLLAMMR",
    silent = TRUE
  )

  # Optimize with multiple restarts (LCA is prone to local optima)
  best_opt <- NULL
  best_obj <- NULL
  best_obj_val <- Inf

  n_starts <- control$n_starts %||% 3

  for (restart in 1:n_starts) {
    if (restart > 1) {
      # Random restart
      tmb_params$item_logits <- matrix(qlogis(runif(n_items * nclass, 0.2, 0.8)),
                                       n_items, nclass)
      tmb_params$class_logits <- rnorm(nclass - 1, 0, 0.5)
      obj <- TMB::MakeADFun(
        data = tmb_data,
        parameters = tmb_params,
        DLL = "GLLAMMR",
        silent = TRUE
      )
    }

    opt <- try(nlminb(
      start = obj$par,
      objective = obj$fn,
      gradient = obj$gr,
      control = list(eval.max = 2000, iter.max = 1000, trace = 0)
    ), silent = TRUE)

    if (!inherits(opt, "try-error") && opt$objective < best_obj_val) {
      best_opt <- opt
      best_obj <- obj   # keep the TMB object whose restart won
      best_obj_val <- opt$objective
    }
  }

  if (is.null(best_opt)) {
    stop("All optimization attempts failed")
  }

  opt <- best_opt
  obj <- best_obj

  # Get standard errors
  sdr <- try(TMB::sdreport(obj), silent = TRUE)

  # Extract parameters
  par_full <- obj$env$last.par.best

  # Item probabilities (logit-parameterized in the template)
  item_probs_matrix <- matrix(
    plogis(par_full[names(par_full) == "item_logits"]),
    nrow = n_items,
    ncol = nclass
  )
  rownames(item_probs_matrix) <- colnames(Y) %||% paste0("Item", 1:n_items)
  colnames(item_probs_matrix) <- paste0("Class", 1:nclass)

  # Class probabilities
  class_logits <- par_full[names(par_full) == "class_logits"]
  sum_exp <- 1 + sum(exp(class_logits))
  class_probs <- c(exp(class_logits) / sum_exp, 1 / sum_exp)
  names(class_probs) <- paste0("Class", 1:nclass)

  # Posterior class membership
  posterior <- matrix(NA, n_obs, nclass)
  for (i in 1:n_obs) {
    for (k in 1:nclass) {
      likelihood <- 1
      for (j in 1:n_items) {
        p <- item_probs_matrix[j, k]
        likelihood <- likelihood * (p^Y[i,j]) * ((1-p)^(1-Y[i,j]))
      }
      posterior[i, k] <- class_probs[k] * likelihood
    }
    posterior[i, ] <- posterior[i, ] / sum(posterior[i, ])
  }

  # Modal class assignment
  modal_class <- apply(posterior, 1, which.max)

  # Construct result
  result <- list(
    nclass = nclass,
    class_probs = class_probs,
    item_probs = item_probs_matrix,
    posterior = posterior,
    modal_class = modal_class,
    logLik = -opt$objective,
    AIC = 2 * opt$objective + 2 * length(obj$par),
    BIC = 2 * opt$objective + log(n_obs) * length(obj$par),
    convergence = list(
      converged = (opt$convergence == 0),
      message = opt$message
    ),
    n_obs = n_obs,
    n_items = n_items,
    tmb_obj = obj,
    tmb_opt = opt,
    tmb_sdr = sdr
  )

  class(result) <- c("gllamm_lca", "gllamm")

  return(result)
}


#' @export
print.gllamm_lca <- function(x, ...) {
  cat("Latent Class Analysis\n\n")
  cat("Number of classes:", x$nclass, "\n")
  cat("Number of observations:", x$n_obs, "\n")
  cat("Number of items:", x$n_items, "\n\n")

  cat("Class probabilities:\n")
  print(round(x$class_probs, 3))

  cat("\nItem response probabilities by class:\n")
  print(round(x$item_probs, 3))

  cat("\nModel fit:\n")
  cat("Log-likelihood:", round(x$logLik, 2), "\n")
  cat("AIC:", round(x$AIC, 2), "\n")
  cat("BIC:", round(x$BIC, 2), "\n")

  invisible(x)
}


#' @export
summary.gllamm_lca <- function(object, ...) {
  print(object)

  cat("\nClass sizes (modal assignment):\n")
  print(table(object$modal_class))

  cat("\nMean posterior probabilities by modal class:\n")
  for (k in 1:object$nclass) {
    modal_k <- object$modal_class == k
    if (sum(modal_k) > 0) {
      cat("Class", k, ":", round(mean(object$posterior[modal_k, k]), 3), "\n")
    }
  }

  invisible(object)
}


# Helper function for NULL coalescing
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
