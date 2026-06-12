#' Weighted isotonic regression (pool-adjacent-violators)
#'
#' Returns the nondecreasing vector minimizing sum(w * (y - x)^2), which is
#' the constrained M-step maximizer for binomial proportions and normal
#' means under a monotone-classes restriction (Croon 1990).
#'
#' @keywords internal
.pava_weighted <- function(y, w) {
  n <- length(y)
  if (n <= 1) return(y)
  vals <- numeric(n); wts <- numeric(n); sz <- integer(n)
  nb <- 0L
  for (i in seq_len(n)) {
    nb <- nb + 1L
    vals[nb] <- y[i]; wts[nb] <- w[i]; sz[nb] <- 1L
    while (nb > 1L && vals[nb - 1L] > vals[nb]) {
      wt <- wts[nb - 1L] + wts[nb]
      vals[nb - 1L] <- (wts[nb - 1L] * vals[nb - 1L] + wts[nb] * vals[nb]) / wt
      wts[nb - 1L] <- wt
      sz[nb - 1L] <- sz[nb - 1L] + sz[nb]
      nb <- nb - 1L
    }
  }
  rep(vals[seq_len(nb)], sz[seq_len(nb)])
}


#' Weighted isotonic regression over a partial order
#'
#' Minimizes sum(w * (y - x)^2) subject to x[a] <= x[b] for every row
#' (a, b) of \code{edges} - the constrained M-step maximizer for binomial
#' proportions and normal means under a partially ordered classes
#' restriction. A chain is solved exactly by pool-adjacent-violators;
#' a general DAG by Dykstra's cyclic projection algorithm (each constraint
#' set is a half-space whose weighted projection is a two-point pool),
#' which converges to the exact projection. Problems here are tiny
#' (length(y) = number of classes), so iteration cost is negligible.
#'
#' @keywords internal
.isotonic_poset <- function(y, w, edges, tol = 1e-12, max_iter = 10000) {
  K <- length(y)
  E <- nrow(edges)
  if (E == 0) return(y)
  # Chain in canonical form: PAVA is exact and direct
  if (E == K - 1 && all(edges[, 1] == seq_len(K - 1)) &&
      all(edges[, 2] == 2:K)) {
    return(.pava_weighted(y, w))
  }
  x <- y
  # Dykstra increments, one pair per constraint
  inc_a <- numeric(E)
  inc_b <- numeric(E)
  for (it in seq_len(max_iter)) {
    delta <- 0
    for (e in seq_len(E)) {
      a <- edges[e, 1]; b <- edges[e, 2]
      xa <- x[a] + inc_a[e]
      xb <- x[b] + inc_b[e]
      if (xa > xb) {
        m <- (w[a] * xa + w[b] * xb) / (w[a] + w[b])
        new_a <- m; new_b <- m
      } else {
        new_a <- xa; new_b <- xb
      }
      inc_a[e] <- xa - new_a
      inc_b[e] <- xb - new_b
      delta <- delta + abs(new_a - x[a]) + abs(new_b - x[b])
      x[a] <- new_a
      x[b] <- new_b
    }
    if (delta < tol) break
  }
  # Guarantee feasibility to numerical precision
  for (e in seq_len(E)) {
    a <- edges[e, 1]; b <- edges[e, 2]
    if (x[a] > x[b]) {
      m <- (w[a] * x[a] + w[b] * x[b]) / (w[a] + w[b])
      x[a] <- m; x[b] <- m
    }
  }
  x
}


#' Topological order of classes under a partial order (Kahn's algorithm)
#'
#' Returns a permutation of 1:K such that every edge (a, b) has a before b.
#' Errors if the order specification contains a cycle.
#'
#' @keywords internal
.topological_order <- function(K, edges) {
  indeg <- integer(K)
  for (e in seq_len(nrow(edges))) {
    indeg[edges[e, 2]] <- indeg[edges[e, 2]] + 1L
  }
  out <- integer(0)
  queue <- which(indeg == 0L)
  while (length(queue)) {
    v <- queue[1]; queue <- queue[-1]
    out <- c(out, v)
    children <- edges[edges[, 1] == v, 2]
    for (ch in children) {
      indeg[ch] <- indeg[ch] - 1L
      if (indeg[ch] == 0L) queue <- c(queue, ch)
    }
  }
  if (length(out) < K) {
    stop("The class ordering constraints contain a cycle; a partial ",
         "order must be acyclic")
  }
  out
}


#' EM estimation for latent class models
#'
#' Classic EM for finite mixtures with conditionally independent indicators
#' (the poLCA algorithm). All M-steps are closed-form (weighted proportions
#' for binary/categorical indicators, weighted moments for gaussian ones),
#' so each iteration is one N x K posterior computation plus a handful of
#' cross-products - the E-step runs on BLAS matrix products. When
#' \code{order_edges} is supplied, binary probabilities and gaussian means
#' are additionally constrained to be nondecreasing along the given class
#' partial order via isotonic regression in the M-step.
#'
#' @keywords internal
fit_lca_em <- function(Y, nclass, item_type, n_cats, weights = NULL,
                       n_starts = 3, max_iter = 1000, tol = 1e-8,
                       order_edges = NULL, item_edges = NULL,
                       structure = "free") {
  ordering <- !is.null(order_edges)
  item_ordering <- !is.null(item_edges)
  rasch <- identical(structure, "rasch")
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

  # Item-monotonicity (and double-monotonicity) constraints couple the
  # binary items: build one combined edge set over the J_b x K probability
  # grid, column-major node index (c - 1) * J_b + j
  grid_edges <- NULL
  if (item_ordering && length(bin_idx)) {
    J_b <- length(bin_idx)
    ge <- list()
    if (ordering) {
      for (e in seq_len(nrow(order_edges))) {
        a <- order_edges[e, 1]; b <- order_edges[e, 2]
        ge[[length(ge) + 1]] <- cbind((a - 1) * J_b + seq_len(J_b),
                                      (b - 1) * J_b + seq_len(J_b))
      }
    }
    for (e in seq_len(nrow(item_edges))) {
      ge[[length(ge) + 1]] <- cbind((seq_len(K) - 1) * J_b +
                                      item_edges[e, 1],
                                    (seq_len(K) - 1) * J_b +
                                      item_edges[e, 2])
    }
    grid_edges <- do.call(rbind, ge)
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
    if (ordering) {
      # Start from a class permutation consistent with the constraint so
      # the first isotonic M-steps don't collapse classes: classes in a
      # topological order of the DAG get increasing initial levels
      topo <- .topological_order(K, order_edges)
      crit <- rep(0, K)
      if (length(bin_idx)) crit <- crit + colMeans(p$bin)
      if (length(gau_idx)) crit <- crit + colMeans(p$mu)
      ord <- integer(K)
      ord[topo] <- order(crit)   # topo position i gets i-th smallest start
      if (length(bin_idx)) p$bin <- p$bin[, ord, drop = FALSE]
      if (length(gau_idx)) {
        p$mu <- p$mu[, ord, drop = FALSE]
        p$sd <- p$sd[, ord, drop = FALSE]
      }
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
        den <- pmax(crossprod(Mb, W), 1e-12)
        praw <- num / den
        if (rasch) {
          # Located latent classes (latent class Rasch): the M-step fits
          # logit pi_jc = theta_c - delta_j by weighted logistic regression
          # on the expected counts
          praw <- .rasch_mstep(num, den)
        } else if (!is.null(grid_edges)) {
          # Item (or double) monotonicity: one isotonic projection of the
          # whole probability grid over the combined partial order
          v <- .isotonic_poset(as.vector(praw), as.vector(den), grid_edges)
          praw <- matrix(v, nrow(praw), K)
        } else if (ordering) {
          # Ordered/partially ordered LCM: the constrained binomial M-step
          # is the weighted isotonic regression of the raw proportions
          # over the class partial order
          for (jj in seq_len(nrow(praw))) {
            praw[jj, ] <- .isotonic_poset(praw[jj, ], den[jj, ], order_edges)
          }
        }
        p$bin <- pmin(pmax(praw, 1e-6), 1 - 1e-6)
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
        m1 <- as.numeric(crossprod(Wo, y)) / den
        m2 <- as.numeric(crossprod(Wo, y^2)) / den
        mu <- if (ordering) .isotonic_poset(m1, den, order_edges) else m1
        # E[(y - mu)^2] at the (possibly constrained) mean
        v <- m2 - 2 * mu * m1 + mu^2
        p$mu[jj, ] <- mu
        p$sd[jj, ] <- sqrt(pmax(v, 1e-6))
      }

      # ---- Ramsay acceleration on successive EM steps (safeguarded) ----
      # (Acceleration is skipped under the Rasch structure: extrapolating
      # cellwise logits would leave the additive model space)
      if (!rasch) {
      th_after <- flatten(p)
      step_now <- th_after - th_before
      step_norm <- sqrt(sum(step_now^2))
      if (step_prev_norm > 0 && step_norm > 0) {
        r <- step_norm / step_prev_norm
        if (r > 0.1 && r < 0.98) {
          gain <- min(r / (1 - r), 20)
          p_revert <- p
          p <- unflatten(th_after + gain * step_now, p)
          if (ordering || item_ordering) {
            # Project the extrapolated point back into the constraint set
            # (isotonic projection); the revert safeguard already protects
            # monotone likelihood
            unit_w <- rep(1, K)
            if (length(bin_idx)) {
              if (!is.null(grid_edges)) {
                v <- .isotonic_poset(as.vector(p$bin),
                                     rep(1, length(p$bin)), grid_edges)
                p$bin <- matrix(v, nrow(p$bin), K)
              } else {
                for (jj in seq_len(nrow(p$bin))) {
                  p$bin[jj, ] <- .isotonic_poset(p$bin[jj, ], unit_w,
                                                 order_edges)
                }
              }
              p$bin <- pmin(pmax(p$bin, 1e-6), 1 - 1e-6)
            }
            if (length(gau_idx) && ordering) {
              for (jj in seq_len(nrow(p$mu))) {
                p$mu[jj, ] <- .isotonic_poset(p$mu[jj, ], unit_w,
                                              order_edges)
              }
            }
          }
          accelerated <- TRUE
        }
      }
      step_prev_norm <- step_norm
      }
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


#' Weighted logistic M-step for the located latent class (Rasch) model
#'
#' Given expected success counts and totals per item x class cell, fits
#' logit pi_jc = theta_c - delta_j by weighted binomial GLM (the
#' complete-data M-step of the latent class Rasch model of Lindsay,
#' Clogg & Grego 1991) and returns the fitted probability matrix.
#'
#' @keywords internal
.rasch_mstep <- function(num, den) {
  J <- nrow(num)
  K <- ncol(num)
  yprop <- as.vector(num / den)
  w <- as.vector(den)
  item_f <- factor(rep(seq_len(J), K))
  class_f <- factor(rep(seq_len(K), each = J))
  X <- stats::model.matrix(~ item_f + class_f)
  fit <- suppressWarnings(
    stats::glm.fit(X, yprop, weights = w, family = stats::binomial()))
  matrix(stats::plogis(as.numeric(X %*% fit$coefficients)), J, K)
}
