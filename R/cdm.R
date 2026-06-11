#' Fit Cognitive Diagnosis Models
#'
#' Fits cognitive diagnosis models (CDMs) for binary responses: latent
#' classes are attribute profiles (combinations of A binary skills), and a
#' Q-matrix declares which attributes each item measures. Item response
#' probabilities depend only on the item's measured attributes (its
#' "reduced profile"), giving the saturated G-DINA / LCDM family, with
#' DINA and DINO as two-group special cases.
#'
#' @param Y Binary response matrix (persons x items; NA allowed)
#' @param Q Binary Q-matrix (items x attributes): \code{Q[j, a] = 1} if
#'   item j measures attribute a. Every row must have at least one 1.
#' @param model Item model:
#'   \itemize{
#'     \item \code{"gdina"} (default): saturated G-DINA - one free
#'       response probability per distinct reduced profile of the item's
#'       measured attributes.
#'     \item \code{"dina"}: conjunctive - two probabilities per item
#'       (guessing for profiles missing any required attribute, 1 - slip
#'       for profiles mastering all of them).
#'     \item \code{"dino"}: disjunctive - mastery of any measured
#'       attribute suffices.
#'   }
#' @param hierarchy Optional attribute hierarchy: a list of attribute
#'   index pairs (or two-column matrix), each pair \code{c(a, b)} meaning
#'   attribute \code{a} is a prerequisite of attribute \code{b}; profiles
#'   with \code{b} mastered but not \code{a} are removed from the latent
#'   space. Attribute names (colnames of Q) may be used instead of
#'   indices.
#' @param monotone Constrain response probabilities to be nondecreasing
#'   in the attributes (default TRUE): mastering more of the measured
#'   attributes can never lower the success probability. Implemented as
#'   weighted isotonic regression over the reduced-profile lattice in the
#'   M-step, so estimation remains closed-form EM.
#' @param weights Optional vector of case weights (one per person)
#' @param control List: \code{n_starts} (default 5), \code{max_iter}
#'   (default 1000), \code{tol} (default 1e-8)
#'
#' @return An object of class \code{gllamm_cdm} with components including
#'   \code{item_params} (per item, P(Y = 1) by reduced profile; for
#'   DINA/DINO also \code{guess}/\code{slip}), \code{profile_probs},
#'   \code{posterior} (persons x profiles), \code{attribute_posteriors}
#'   (persons x attributes marginal mastery probabilities),
#'   \code{modal_profile}, and \code{logLik}/\code{AIC}/\code{BIC}.
#'
#' @details
#' Estimation is marginal maximum likelihood via EM with closed-form
#' M-steps (weighted proportions pooled over the item's reduced-profile
#' groups, isotonically projected when \code{monotone = TRUE}) and
#' safeguarded Ramsay acceleration. AIC/BIC use the nominal parameter
#' count; likelihood-ratio tests against less constrained models have
#' non-standard null distributions when monotonicity binds.
#'
#' @examples
#' \dontrun{
#' Q <- rbind(c(1, 0), c(0, 1), c(1, 1), c(1, 0), c(0, 1))
#' fit <- fit_cdm(Y, Q, model = "dina")
#' summary(fit)
#' fit$attribute_posteriors  # P(mastery) per person and attribute
#' }
#'
#' @export
fit_cdm <- function(Y, Q, model = c("gdina", "dina", "dino"),
                    hierarchy = NULL, monotone = TRUE,
                    weights = NULL, control = list()) {
  model <- match.arg(model)

  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  n_obs <- nrow(Y)
  n_items <- ncol(Y)
  if (nrow(Q) != n_items) {
    stop("Q must have one row per item (", n_items, "); found ", nrow(Q))
  }
  if (!all(Q %in% c(0, 1))) stop("Q must be a binary (0/1) matrix")
  storage.mode(Q) <- "integer"
  if (any(rowSums(Q) == 0)) {
    stop("Every item must measure at least one attribute; Q rows ",
         paste(which(rowSums(Q) == 0), collapse = ", "), " are all zero")
  }
  vals <- Y[!is.na(Y)]
  if (!all(vals %in% c(0, 1))) {
    stop("Y must be binary (0/1, NA allowed)")
  }
  if (!is.null(weights)) {
    if (length(weights) != n_obs) {
      stop("weights length (", length(weights),
           ") must equal number of persons (", n_obs, ")")
    }
    if (anyNA(weights) || any(weights < 0)) {
      stop("weights must be non-negative and complete")
    }
  }

  n_attr <- ncol(Q)
  attr_names <- colnames(Q) %||% paste0("A", seq_len(n_attr))
  colnames(Q) <- attr_names
  item_names <- colnames(Y) %||% paste0("Item", seq_len(n_items))

  # ---- Latent space: attribute profiles, optionally hierarchy-pruned ----
  profiles <- as.matrix(expand.grid(rep(list(0:1), n_attr)))
  dimnames(profiles) <- NULL
  storage.mode(profiles) <- "integer"

  if (!is.null(hierarchy)) {
    hedges <- if (is.list(hierarchy)) {
      if (!all(vapply(hierarchy, length, 0L) == 2)) {
        stop("hierarchy pairs must each have exactly two attributes")
      }
      do.call(rbind, lapply(hierarchy, function(e) {
        if (is.character(e)) match(e, attr_names) else as.integer(e)
      }))
    } else if (is.matrix(hierarchy) && ncol(hierarchy) == 2) {
      if (is.character(hierarchy)) {
        matrix(match(hierarchy, attr_names), ncol = 2)
      } else {
        hierarchy
      }
    } else {
      stop("hierarchy must be a list or two-column matrix of attribute pairs")
    }
    storage.mode(hedges) <- "integer"
    if (anyNA(hedges) || any(hedges < 1) || any(hedges > n_attr)) {
      stop("hierarchy attributes must be valid indices or names ",
           "(prerequisite first)")
    }
    if (any(hedges[, 1] == hedges[, 2])) {
      stop("hierarchy pairs must relate two distinct attributes")
    }
    .topological_order(n_attr, hedges)   # errors on cycles
    keep <- rep(TRUE, nrow(profiles))
    for (e in seq_len(nrow(hedges))) {
      keep <- keep & !(profiles[, hedges[e, 2]] == 1L &
                         profiles[, hedges[e, 1]] == 0L)
    }
    profiles <- profiles[keep, , drop = FALSE]
  }

  K <- nrow(profiles)
  profile_labels <- apply(profiles, 1, paste, collapse = "")

  # ---- Per-item reduced-profile groups and monotonicity edges ----
  group_id <- matrix(0L, n_items, K)
  group_patterns <- vector("list", n_items)   # G_j x k_j pattern matrices
  group_edges <- vector("list", n_items)
  for (j in seq_len(n_items)) {
    measured <- which(Q[j, ] == 1L)
    red <- profiles[, measured, drop = FALSE]
    if (model == "dina") {
      eta <- as.integer(rowSums(red) == length(measured))
      pat <- matrix(c(rep(0L, length(measured)), rep(1L, length(measured))),
                    nrow = 2, byrow = TRUE)
      present <- sort(unique(eta))
      group_id[j, ] <- match(eta, present)
      group_patterns[[j]] <- pat[present + 1L, , drop = FALSE]
      group_edges[[j]] <- if (length(present) == 2) {
        cbind(1L, 2L)
      } else {
        matrix(integer(0), 0, 2)
      }
    } else if (model == "dino") {
      eta <- as.integer(rowSums(red) > 0L)
      pat <- matrix(c(rep(0L, length(measured)), rep(1L, length(measured))),
                    nrow = 2, byrow = TRUE)
      present <- sort(unique(eta))
      group_id[j, ] <- match(eta, present)
      group_patterns[[j]] <- pat[present + 1L, , drop = FALSE]
      group_edges[[j]] <- if (length(present) == 2) {
        cbind(1L, 2L)
      } else {
        matrix(integer(0), 0, 2)
      }
    } else {
      key <- apply(red, 1, paste, collapse = "")
      ukey <- unique(key)
      upat <- red[match(ukey, key), , drop = FALSE]
      ord <- order(rowSums(upat), ukey)   # display order: by level
      upat <- upat[ord, , drop = FALSE]
      ukey <- ukey[ord]
      group_id[j, ] <- match(key, ukey)
      group_patterns[[j]] <- upat
      G <- length(ukey)
      ed <- which(outer(seq_len(G), seq_len(G), Vectorize(function(a, b) {
        a != b && all(upat[a, ] <= upat[b, ]) && any(upat[a, ] < upat[b, ])
      })), arr.ind = TRUE)
      group_edges[[j]] <- if (nrow(ed)) {
        cbind(as.integer(ed[, 1]), as.integer(ed[, 2]))
      } else {
        matrix(integer(0), 0, 2)
      }
    }
  }
  n_groups <- vapply(group_patterns, nrow, 0L)
  group_levels <- lapply(group_patterns, rowSums)

  em <- fit_cdm_em(Y, profiles, group_id, n_groups, group_edges,
                   group_levels = group_levels,
                   monotone = monotone, weights = weights,
                   n_starts = control$n_starts %||% 3,
                   max_iter = control$max_iter %||% 2000,
                   tol = control$tol %||% 1e-7)

  # ---- Assemble results ----
  item_params <- vector("list", n_items)
  names(item_params) <- item_names
  for (j in seq_len(n_items)) {
    pr <- em$params$items[[j]]
    names(pr) <- apply(group_patterns[[j]], 1, paste, collapse = "")
    entry <- list(measured = attr_names[Q[j, ] == 1L], prob = pr)
    if (model %in% c("dina", "dino") && n_groups[j] == 2) {
      entry$guess <- unname(pr[1])
      entry$slip <- unname(1 - pr[2])
    }
    item_params[[j]] <- entry
  }

  posterior <- em$posterior
  colnames(posterior) <- profile_labels
  attribute_posteriors <- posterior %*% profiles
  colnames(attribute_posteriors) <- attr_names
  modal_idx <- apply(posterior, 1, which.max)

  profile_probs <- em$params$pi
  names(profile_probs) <- profile_labels

  n_par <- (K - 1) + sum(n_groups)

  result <- list(
    model = model,
    nclass = K,
    n_attributes = n_attr,
    Q = Q,
    profiles = profiles,
    profile_labels = profile_labels,
    hierarchy = hierarchy,
    monotone = monotone,
    item_params = item_params,
    profile_probs = profile_probs,
    posterior = posterior,
    attribute_posteriors = attribute_posteriors,
    modal_profile = profile_labels[modal_idx],
    modal_attributes = profiles[modal_idx, , drop = FALSE],
    logLik = em$loglik,
    AIC = -2 * em$loglik + 2 * n_par,
    BIC = -2 * em$loglik + log(n_obs) * n_par,
    n_par = n_par,
    convergence = list(
      converged = em$converged,
      message = if (em$converged) "EM converged"
                else "Maximum EM iterations reached",
      iterations = em$iterations
    ),
    n_obs = n_obs,
    n_items = n_items
  )
  class(result) <- c("gllamm_cdm", "gllamm")
  result
}


#' EM estimation for cognitive diagnosis models
#'
#' Marginal ML for CDMs: the E-step is the standard finite-mixture
#' posterior over attribute profiles (BLAS matrix products); the M-step
#' pools expected counts over each item's reduced-profile groups (the
#' closed-form weighted proportion) and, when monotone, projects with
#' isotonic regression over the reduced-profile lattice. Safeguarded
#' Ramsay acceleration as in \code{fit_lca_em}.
#'
#' @keywords internal
fit_cdm_em <- function(Y, profiles, group_id, n_groups, group_edges,
                       group_levels, monotone = TRUE, weights = NULL,
                       n_starts = 5, max_iter = 1000, tol = 1e-8) {
  n_obs <- nrow(Y)
  n_items <- ncol(Y)
  K <- nrow(profiles)
  obs_mask <- !is.na(Y)
  w <- if (is.null(weights)) rep(1, n_obs) else as.numeric(weights)

  Yb <- ifelse(obs_mask, Y, 0)
  Mb <- obs_mask * 1
  complete <- all(obs_mask)

  # Class membership lists per item group (M-step pooling)
  grp_members <- lapply(seq_len(n_items), function(j) {
    lapply(seq_len(n_groups[j]), function(g) which(group_id[j, ] == g))
  })

  init_params <- function() {
    base <- pmin(pmax(colSums(Yb) / pmax(colSums(Mb), 1), 0.1), 0.9)
    items <- vector("list", n_items)
    for (j in seq_len(n_items)) {
      # Monotone-friendly starts: spread the groups by their lattice level
      # (number of mastered measured attributes) around the item base rate
      G <- n_groups[j]
      lev <- group_levels[[j]]
      span <- runif(1, 1.5, 3.5)
      eta <- qlogis(base[j]) + span * (lev / max(lev, 1) - 0.5) +
        rnorm(G, 0, 0.25)
      items[[j]] <- pmin(pmax(plogis(eta), 0.05), 0.95)
    }
    pi0 <- runif(K, 0.5, 1.5)
    list(pi = pi0 / sum(pi0), items = items)
  }

  expand_probs <- function(p) {
    Pf <- matrix(0, n_items, K)
    for (j in seq_len(n_items)) {
      Pf[j, ] <- p$items[[j]][group_id[j, ]]
    }
    Pf
  }

  class_loglik <- function(p) {
    Pf <- expand_probs(p)
    lq <- log1p(-Pf)
    if (complete) {
      # One BLAS product: y*log(p) + (1-y)*log(1-p) = y*logit(p) + log(1-p)
      sweep(Yb %*% (log(Pf) - lq), 2, colSums(lq), "+")
    } else {
      Yb %*% log(Pf) + (Mb - Yb) %*% lq
    }
  }

  row_max <- function(M) {
    m <- M[, 1]
    for (k in 2:ncol(M)) m <- pmax(m, M[, k])
    m
  }

  flatten <- function(p) c(log(p$pi), qlogis(unlist(p$items)))
  unflatten <- function(v, p) {
    pi_raw <- exp(v[1:K]); p$pi <- pi_raw / sum(pi_raw)
    k <- K
    for (j in seq_len(n_items)) {
      G <- n_groups[j]
      pj <- plogis(v[k + seq_len(G)])
      if (monotone && nrow(group_edges[[j]])) {
        # Project extrapolations back into the constraint set
        pj <- .isotonic_poset(pj, rep(1, G), group_edges[[j]])
      }
      p$items[[j]] <- pmin(pmax(pj, 1e-6), 1 - 1e-6)
      k <- k + G
    }
    p
  }

  run_em <- function(p) {
    loglik_old <- -Inf
    loglik <- NA_real_
    converged <- FALSE
    iter <- 0
    th_before <- NULL
    step_prev_norm <- -1
    accelerated <- FALSE
    p_revert <- NULL

    for (iter in seq_len(max_iter)) {
      # ---- E-step ----
      Lw <- sweep(class_loglik(p), 2, log(p$pi), "+")
      m <- row_max(Lw)
      W <- exp(Lw - m)
      rs <- rowSums(W)
      loglik <- sum(w * (m + log(rs)))
      W <- (W / rs) * w

      if (accelerated && loglik < loglik_old) {
        p <- p_revert
        accelerated <- FALSE
        next
      }
      accelerated <- FALSE

      if (abs(loglik - loglik_old) < tol) {
        converged <- TRUE
        break
      }
      loglik_old <- loglik
      th_before <- flatten(p)

      # ---- closed-form M-steps ----
      # Floor pi: with 2^A profiles many prevalences are legitimately ~0,
      # and an exact zero makes log(pi) infinite (breaks acceleration)
      Nk <- colSums(W)
      pi_new <- Nk / sum(Nk)
      pi_new[pi_new < 1e-12] <- 1e-12
      p$pi <- pi_new / sum(pi_new)
      num <- crossprod(Yb, W)                      # n_items x K
      # With complete data the per-class denominator is just Nk for every
      # item; the second crossprod is only needed under missingness
      den <- if (complete) NULL else crossprod(Mb, W)
      for (j in seq_len(n_items)) {
        dj <- if (complete) Nk else den[j, ]
        ng <- vapply(grp_members[[j]], function(idx) sum(num[j, idx]), 0)
        dg <- vapply(grp_members[[j]], function(idx) sum(dj[idx]), 0)
        dg[dg < 1e-12] <- 1e-12
        pg <- ng / dg
        if (monotone && nrow(group_edges[[j]])) {
          pg <- .isotonic_poset(pg, dg, group_edges[[j]])
        }
        pg[pg < 1e-6] <- 1e-6
        pg[pg > 1 - 1e-6] <- 1 - 1e-6
        p$items[[j]] <- pg
      }

      # ---- Ramsay acceleration (safeguarded; unflatten projects) ----
      th_after <- flatten(p)
      step_now <- th_after - th_before
      step_norm <- sqrt(sum(step_now^2))
      if (step_prev_norm > 0 && step_norm > 0) {
        r <- step_norm / step_prev_norm
        if (r > 0.1 && r < 0.98) {
          gain <- min(r / (1 - r), 20)
          p_revert <- p
          p <- unflatten(th_after + gain * step_now, p)
          accelerated <- TRUE
        }
      }
      step_prev_norm <- step_norm
    }
    list(p = p, loglik = loglik, converged = converged, iter = iter)
  }

  best <- NULL
  for (s in seq_len(n_starts)) {
    fit <- run_em(init_params())
    if (is.null(best) || fit$loglik > best$loglik) best <- fit
  }

  p <- best$p
  Lw <- sweep(class_loglik(p), 2, log(p$pi), "+")
  m <- row_max(Lw)
  W <- exp(Lw - m)
  posterior <- W / rowSums(W)

  list(params = p, loglik = best$loglik, posterior = posterior,
       converged = best$converged, iterations = best$iter)
}


#' @export
print.gllamm_cdm <- function(x, ...) {
  cat("Cognitive Diagnosis Model (",
      switch(x$model, gdina = "saturated G-DINA", dina = "DINA",
             dino = "DINO"), ")\n\n", sep = "")
  cat("Attributes:", x$n_attributes,
      "| Latent profiles:", x$nclass,
      if (!is.null(x$hierarchy)) "(hierarchy-pruned)" else "", "\n")
  cat("Persons:", x$n_obs, "| Items:", x$n_items,
      "| Monotone:", x$monotone, "\n\n")

  cat("Profile prevalences:\n")
  print(round(x$profile_probs, 3))

  cat("\nAttribute mastery prevalences:\n")
  prev <- as.numeric(x$profile_probs %*% x$profiles)
  names(prev) <- colnames(x$Q)
  print(round(prev, 3))

  if (x$model %in% c("dina", "dino")) {
    gs <- t(vapply(x$item_params, function(e)
      c(guess = e$guess %||% NA_real_, slip = e$slip %||% NA_real_),
      c(guess = 0, slip = 0)))
    cat("\nItem guess/slip:\n")
    print(round(gs, 3))
  }

  cat("\nLog-likelihood:", round(x$logLik, 2),
      "| AIC:", round(x$AIC, 2), "| BIC:", round(x$BIC, 2), "\n")
  invisible(x)
}


#' @export
summary.gllamm_cdm <- function(object, ...) {
  print(object)

  if (object$model == "gdina") {
    cat("\nItem success probabilities by reduced profile:\n")
    for (j in seq_along(object$item_params)) {
      e <- object$item_params[[j]]
      cat("  ", names(object$item_params)[j],
          " (", paste(e$measured, collapse = ","), "): ",
          paste(sprintf("%s=%.3f", names(e$prob), e$prob), collapse = "  "),
          "\n", sep = "")
    }
  }

  cat("\nClassification certainty (mean max posterior):",
      round(mean(apply(object$posterior, 1, max)), 3), "\n")
  cat("Modal profile distribution:\n")
  print(table(object$modal_profile))
  invisible(object)
}
