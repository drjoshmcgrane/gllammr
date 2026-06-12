#' EM estimation for latent class models
#'
#' Classic EM for finite mixtures with conditionally independent indicators
#' (the poLCA algorithm). All M-steps are closed-form (weighted proportions
#' for binary/categorical indicators, weighted moments for gaussian ones),
#' so each iteration is one N x K posterior computation plus a handful of
#' cross-products - the E-step runs on BLAS matrix products.
#'
#' @keywords internal
fit_lca_em <- function(Y, nclass, item_type, n_cats, weights = NULL,
                       n_starts = 3, max_iter = 1000, tol = 1e-8) {
  n_obs <- nrow(Y)
  n_items <- ncol(Y)
  K <- nclass
  obs_mask <- !is.na(Y)
  w <- if (is.null(weights)) rep(1, n_obs) else as.numeric(weights)

  bin_idx <- which(item_type == 0L)
  cat_idx <- which(item_type == 1L)
  gau_idx <- which(item_type == 2L)

  # Pre-extracted blocks (missing set to 0 with masks)
  Yb <- NULL; Mb <- NULL
  if (length(bin_idx)) {
    Yb <- ifelse(obs_mask[, bin_idx, drop = FALSE],
                 Y[, bin_idx, drop = FALSE], 0)
    Mb <- obs_mask[, bin_idx, drop = FALSE] * 1
  }

  init_params <- function() {
    p <- list(pi = rep(1 / K, K))
    if (length(bin_idx)) {
      marg <- colSums(Yb) / pmax(colSums(Mb), 1)
      p$bin <- matrix(pmin(pmax(marg * runif(length(bin_idx) * K, 0.6, 1.4),
                                0.05), 0.95),
                      length(bin_idx), K)
    }
    if (length(cat_idx)) {
      p$cat <- lapply(cat_idx, function(j) {
        Kj <- n_cats[j]
        m <- matrix(runif(Kj * K, 0.5, 1.5), Kj, K)
        sweep(m, 2, colSums(m), "/")
      })
    }
    if (length(gau_idx)) {
      mu0 <- vapply(gau_idx, function(j) mean(Y[obs_mask[, j], j]), 0)
      sd0 <- vapply(gau_idx, function(j) sd(Y[obs_mask[, j], j]), 0)
      p$mu <- outer(mu0, rep(1, K)) +
        outer(sd0, rnorm(K, 0, 0.7))
      p$sd <- matrix(pmax(sd0, 0.1), length(gau_idx), K)
    }
    p
  }

  # N x K log-likelihood under each class
  class_loglik <- function(p) {
    L <- matrix(0, n_obs, K)
    if (length(bin_idx)) {
      lp <- log(p$bin); lq <- log1p(-p$bin)        # (J_b x K)
      L <- L + Yb %*% lp + (Mb - Yb) %*% lq
    }
    for (jj in seq_along(cat_idx)) {
      j <- cat_idx[jj]
      logp <- log(pmax(p$cat[[jj]], 1e-12))        # Kj x K
      obs <- obs_mask[, j]
      L[obs, ] <- L[obs, ] + logp[Y[obs, j], , drop = FALSE]
    }
    for (jj in seq_along(gau_idx)) {
      j <- gau_idx[jj]
      obs <- obs_mask[, j]
      y <- Y[obs, j]
      for (k in 1:K) {
        L[obs, k] <- L[obs, k] + dnorm(y, p$mu[jj, k], p$sd[jj, k], log = TRUE)
      }
    }
    L
  }

  row_max <- function(M) {
    m <- M[, 1]
    for (k in 2:ncol(M)) m <- pmax(m, M[, k])
    m
  }

  # Unconstrained flatten/unflatten for Ramsay extrapolation
  flatten <- function(p) {
    v <- log(p$pi)
    if (length(bin_idx)) v <- c(v, qlogis(p$bin))
    for (m in p$cat) v <- c(v, log(m))
    if (length(gau_idx)) v <- c(v, p$mu, log(p$sd))
    v
  }
  unflatten <- function(v, p) {
    k <- K
    pi_raw <- exp(v[1:K]); p$pi <- pi_raw / sum(pi_raw)
    if (length(bin_idx)) {
      nb <- length(bin_idx) * K
      p$bin <- matrix(plogis(v[k + 1:nb]), length(bin_idx), K)
      p$bin <- pmin(pmax(p$bin, 1e-6), 1 - 1e-6)
      k <- k + nb
    }
    for (jj in seq_along(cat_idx)) {
      Kj <- n_cats[cat_idx[jj]]
      m <- matrix(exp(v[k + 1:(Kj * K)]), Kj, K)
      p$cat[[jj]] <- sweep(m, 2, colSums(m), "/")
      k <- k + Kj * K
    }
    if (length(gau_idx)) {
      ng <- length(gau_idx) * K
      p$mu <- matrix(v[k + 1:ng], length(gau_idx), K); k <- k + ng
      p$sd <- matrix(pmax(exp(v[k + 1:ng]), 1e-4), length(gau_idx), K)
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
        # Safeguard: extrapolation lost likelihood - revert to plain EM point
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
      Nk <- colSums(W)
      p$pi <- Nk / sum(Nk)
      if (length(bin_idx)) {
        num <- crossprod(Yb, W)                    # J_b x K
        den <- crossprod(Mb, W)
        p$bin <- pmin(pmax(num / pmax(den, 1e-12), 1e-6), 1 - 1e-6)
      }
      for (jj in seq_along(cat_idx)) {
        j <- cat_idx[jj]
        Kj <- n_cats[j]
        obs <- obs_mask[, j]
        num <- matrix(0, Kj, K)
        for (cc in 1:Kj) {
          sel <- obs & Y[, j] == cc
          if (any(sel)) num[cc, ] <- colSums(W[sel, , drop = FALSE])
        }
        den <- colSums(W[obs, , drop = FALSE])
        p$cat[[jj]] <- pmax(sweep(num, 2, pmax(den, 1e-12), "/"), 1e-9)
      }
      for (jj in seq_along(gau_idx)) {
        j <- gau_idx[jj]
        obs <- obs_mask[, j]
        Wo <- W[obs, , drop = FALSE]
        y <- Y[obs, j]
        den <- pmax(colSums(Wo), 1e-12)
        mu <- as.numeric(crossprod(Wo, y)) / den
        v <- as.numeric(crossprod(Wo, y^2)) / den - mu^2
        p$mu[jj, ] <- mu
        p$sd[jj, ] <- sqrt(pmax(v, 1e-6))
      }

      # ---- Ramsay acceleration on successive EM steps (safeguarded) ----
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

  # Final posteriors at the best solution
  p <- best$p
  Lw <- sweep(class_loglik(p), 2, log(p$pi), "+")
  m <- row_max(Lw)
  W <- exp(Lw - m)
  posterior <- W / rowSums(W)

  list(params = p, loglik = best$loglik, posterior = posterior,
       converged = best$converged, iterations = best$iter,
       bin_idx = bin_idx, cat_idx = cat_idx, gau_idx = gau_idx)
}
