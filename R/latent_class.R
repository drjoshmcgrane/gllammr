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

  # Detect indicator types: binary (0/1), categorical (integers 1..K, K > 2),
  # continuous (anything else numeric)
  item_type <- integer(n_items)
  n_cats <- integer(n_items)
  for (j in seq_len(n_items)) {
    v <- Y[!is.na(Y[, j]), j]
    if (all(v %in% c(0, 1))) {
      item_type[j] <- 0L
    } else if (all(v == round(v)) && min(v) == 1 && length(unique(v)) > 2 &&
               max(v) == length(unique(v))) {
      item_type[j] <- 1L
      n_cats[j] <- as.integer(max(v))
    } else {
      item_type[j] <- 2L
    }
  }
  max_cats <- max(c(n_cats, 2L))

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
    item_type = item_type,
    n_cats = n_cats,
    max_cats = as.integer(max_cats),
    weights = weights_vec
  )

  # Initialize parameters for every block (unused blocks are mapped off)
  init_params <- function() {
    marg <- colMeans(Y)
    item_logits_init <- matrix(0, n_items, nclass)
    cat_logits_init <- array(0, dim = c(n_items, nclass, max_cats - 1))
    means_init <- matrix(0, n_items, nclass)
    log_sds_init <- matrix(0, n_items, nclass)

    for (j in seq_len(n_items)) {
      if (item_type[j] == 0L) {
        p0 <- pmin(pmax(marg[j], 0.05), 0.95)
        item_logits_init[j, ] <- qlogis(pmin(pmax(
          p0 * runif(nclass, 0.7, 1.3), 0.05), 0.95))
      } else if (item_type[j] == 1L) {
        cat_logits_init[j, , ] <- rnorm(nclass * (max_cats - 1), 0, 0.3)
      } else {
        v <- Y[, j]
        means_init[j, ] <- mean(v) + stats::sd(v) * rnorm(nclass, 0, 0.7)
        log_sds_init[j, ] <- log(max(stats::sd(v), 0.1))
      }
    }
    list(item_logits = item_logits_init,
         cat_logits = cat_logits_init,
         item_means = means_init,
         item_log_sds = log_sds_init,
         class_logits = rnorm(nclass - 1, 0, 0.2))
  }

  if (is.null(start)) {
    tmb_params <- init_params()
  } else {
    tmb_params <- start
  }

  # Map off parameter blocks for absent indicator types (per item row)
  block_map <- function(active_rows, dims) {
    m <- array(seq_len(prod(dims)), dim = dims)
    flat <- as.vector(m)
    keep <- as.vector(slice.index(m, 1) %in% which(active_rows))
    flat[!keep] <- NA
    factor(flat)
  }
  tmb_map <- list(
    item_logits = block_map(item_type == 0L, c(n_items, nclass)),
    cat_logits = block_map(item_type == 1L, c(n_items, nclass, max_cats - 1)),
    item_means = block_map(item_type == 2L, c(n_items, nclass)),
    item_log_sds = block_map(item_type == 2L, c(n_items, nclass))
  )

  # Create TMB object
  tmb_data$model_name <- "latent_class"
  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    map = tmb_map,
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
      # Random restart across all parameter blocks
      tmb_params <- init_params()
      obj <- TMB::MakeADFun(
        data = tmb_data,
        parameters = tmb_params,
        map = tmb_map,
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

  # Extract parameters: parList() reassembles full-shaped blocks (maps respected)
  pl <- obj$env$parList()
  item_names <- colnames(Y) %||% paste0("Item", 1:n_items)

  # Binary item probabilities (rows for non-binary items are reported as NA)
  item_probs_matrix <- plogis(pl$item_logits)
  item_probs_matrix[item_type != 0L, ] <- NA
  dimnames(item_probs_matrix) <- list(item_names, paste0("Class", 1:nclass))

  # Categorical item probabilities: per-item list of K x nclass matrices
  cat_probs <- NULL
  if (any(item_type == 1L)) {
    cat_probs <- list()
    for (j in which(item_type == 1L)) {
      K <- n_cats[j]
      eta <- rbind(0, matrix(pl$cat_logits[j, , seq_len(K - 1)],
                             ncol = nclass, byrow = TRUE))
      pj <- apply(eta, 2, function(e) exp(e - max(e)) / sum(exp(e - max(e))))
      dimnames(pj) <- list(paste0("Cat", 1:K), paste0("Class", 1:nclass))
      cat_probs[[item_names[j]]] <- pj
    }
  }

  # Gaussian indicator parameters
  gaussian_params <- NULL
  if (any(item_type == 2L)) {
    idx <- which(item_type == 2L)
    gaussian_params <- list(
      means = matrix(pl$item_means[idx, ], nrow = length(idx),
                     dimnames = list(item_names[idx], paste0("Class", 1:nclass))),
      sds = matrix(exp(pl$item_log_sds[idx, ]), nrow = length(idx),
                   dimnames = list(item_names[idx], paste0("Class", 1:nclass)))
    )
  }

  # Class probabilities
  class_logits <- pl$class_logits
  sum_exp <- 1 + sum(exp(class_logits))
  class_probs <- c(exp(class_logits) / sum_exp, 1 / sum_exp)
  names(class_probs) <- paste0("Class", 1:nclass)

  # Posterior class membership (log-space, mixed indicator types)
  log_lik_class <- matrix(0, n_obs, nclass)
  for (k in 1:nclass) {
    ll <- rep(log(class_probs[k]), n_obs)
    for (j in 1:n_items) {
      if (item_type[j] == 0L) {
        x <- pl$item_logits[j, k]
        ll <- ll + Y[, j] * x - log1p(exp(x))
      } else if (item_type[j] == 1L) {
        K <- n_cats[j]
        eta <- c(0, pl$cat_logits[j, k, seq_len(K - 1)])
        logp <- eta - log(sum(exp(eta - max(eta)))) - max(eta)
        ll <- ll + logp[Y[, j]]
      } else {
        ll <- ll + dnorm(Y[, j], pl$item_means[j, k],
                         exp(pl$item_log_sds[j, k]), log = TRUE)
      }
    }
    log_lik_class[, k] <- ll
  }
  m_row <- apply(log_lik_class, 1, max)
  posterior <- exp(log_lik_class - m_row)
  posterior <- posterior / rowSums(posterior)

  # Modal class assignment
  modal_class <- apply(posterior, 1, which.max)

  # Construct result
  result <- list(
    nclass = nclass,
    class_probs = class_probs,
    item_probs = item_probs_matrix,
    cat_probs = cat_probs,
    gaussian_params = gaussian_params,
    item_type = item_type,
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
