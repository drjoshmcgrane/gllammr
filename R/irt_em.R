#' MML-EM estimation for IRT models (Bock-Aitkin)
#'
#' Marginal maximum likelihood via the EM algorithm with fixed Gauss-Hermite
#' quadrature over the ability distribution - the algorithm class used by
#' mirt and TAM. The E-step computes person-by-node posterior weights in one
#' matrix product; M-steps are independent small optimizations per item.
#'
#' Identification matches the Laplace path: 2PL/3PL/GRM/GPCM/NRM fix the
#' ability SD at 1; Rasch and PCM keep a free ability SD via the
#' common-slope equivalence (theta ~ N(0, sigma) with unit slope is
#' theta ~ N(0,1) with shared slope sigma; difficulties back-transform as
#' b = d * sigma).
#'
#' @keywords internal
fit_irt_em <- function(response_matrix, model, weights = NULL,
                       mc_items = NULL, quad_points = 61,
                       max_iter = 500, tol = 1e-4, control = list()) {

  n_persons <- nrow(response_matrix)
  n_items <- ncol(response_matrix)
  Y <- as.matrix(response_matrix)
  obs_mask <- !is.na(Y)

  is_poly <- model %in% c("GRM", "PCM", "GPCM", "NRM")
  free_sigma <- model %in% c("Rasch", "PCM")   # common-slope equivalence
  # On the internal N(0,1) scale, PCM with a free ability SD is GPCM with a
  # shared slope (irt_category_probs fixes PCM's slope at 1 by design)
  prob_model <- if (model == "PCM") "GPCM" else model

  w_person <- if (is.null(weights)) rep(1, n_persons) else as.numeric(weights)

  # ---- Quadrature against the standard normal ----
  # Equally-spaced rectangular grid with normal-density weights (the
  # classical Bock-Aitkin / mirt scheme). For long tests the person
  # posteriors are very sharp; a uniform grid keeps resolution where the
  # data lives, whereas Gauss-Hermite nodes spread far into the tails.
  nodes <- seq(-6, 6, length.out = quad_points)
  log_A <- dnorm(nodes, log = TRUE)
  log_A <- log_A - log(sum(exp(log_A)))        # normalize: sum(A) = 1
  Q <- quad_points

  if (is_poly) {
    vp <- validate_poly_responses(Y)
    Y <- vp$response_matrix
    obs_mask <- !is.na(Y)
    K_per_item <- vp$n_categories_per_item
    max_K <- vp$max_categories

    # The entire polytomous EM loop runs in compiled C++ (E-step posteriors,
    # expected counts, per-item damped-Newton M-steps, common-slope update)
    model_code <- switch(model, GRM = 1L, PCM = 3L, GPCM = 3L, NRM = 4L)
    res <- .Call(C_em_poly,
                 matrix(as.integer(Y), n_persons, n_items),
                 as.integer(K_per_item),
                 as.numeric(w_person),
                 as.numeric(nodes),
                 as.numeric(log_A),
                 model_code,
                 as.integer(free_sigma),
                 as.integer(max_iter),
                 as.numeric(tol))

    sigma_hat <- if (free_sigma) res$a[1] else 1
    thresholds <- if (free_sigma) {
      lapply(res$thresholds, function(tr) tr * sigma_hat)
    } else res$thresholds
    names(thresholds) <- paste0("Item", 1:n_items)
    discrimination <- if (free_sigma) rep(1, n_items) else res$a
    names(discrimination) <- paste0("Item", 1:n_items)

    theta_hat <- res$theta * sigma_hat
    names(theta_hat) <- paste0("Person", seq_len(n_persons))

    n_par <- sum(K_per_item - 1) +
      (if (model %in% c("GRM", "GPCM", "NRM")) n_items else 0) +
      as.integer(free_sigma)

    result <- list(
      model = model,
      method = "EM",
      item_parameters = list(
        thresholds = thresholds,
        discrimination = discrimination
      ),
      person_abilities = theta_hat,
      ability_sd = sigma_hat,
      logLik = res$logLik,
      AIC = -2 * res$logLik + 2 * n_par,
      BIC = -2 * res$logLik + log(n_persons) * n_par,
      convergence = list(
        converged = res$converged,
        message = if (res$converged) "EM converged" else "Maximum EM iterations reached",
        iterations = res$iterations
      ),
      n_persons = n_persons,
      n_items = n_items,
      n_categories_per_item = K_per_item,
      max_categories = max_K,
      quad_points = quad_points
    )
    class(result) <- c("gllamm_irt_poly", "gllamm_irt", "gllamm")
    return(result)
  }

  # ---- Item probability matrices at the quadrature nodes ----
  # Dichotomous: returns J x Q matrix of P(Y=1 | z_q).
  # Polytomous: list over items of Q x K category-probability matrices.
  item_probs_at_nodes <- function(par) {
    if (!is_poly) {
      P <- matrix(0, n_items, Q)
      for (j in seq_len(n_items)) {
        p <- plogis(par$a[j] * (nodes - par$d[j]))
        if (model == "3PL" && par$mc[j] == 1) {
          p <- par$c[j] + (1 - par$c[j]) * p
        }
        P[j, ] <- p
      }
      pmin(pmax(P, 1e-10), 1 - 1e-10)
    } else {
      lapply(seq_len(n_items), function(j) {
        irt_category_probs(prob_model, nodes, par$thresholds[[j]], par$a[j])
      })
    }
  }

  # ---- Person-by-node log-likelihood (vectorized) ----
  loglik_matrix <- function(P) {
    if (!is_poly) {
      Y0 <- ifelse(obs_mask, Y, 0)
      M <- obs_mask * 1
      # N x Q: sum_j [ y log p + (1-y) log(1-p) ] over observed items
      Y0 %*% log(P) + (M - Y0) %*% log(1 - P)
    } else {
      L <- matrix(0, n_persons, Q)
      for (j in seq_len(n_items)) {
        logPj <- log(pmax(P[[j]], 1e-12))       # Q x K
        obs_j <- obs_mask[, j]
        L[obs_j, ] <- L[obs_j, ] + t(logPj[, Y[obs_j, j], drop = FALSE])
      }
      L
    }
  }

  # ---- Initial values ----
  par <- list(a = rep(1, n_items), d = numeric(n_items),
              c = rep(0, n_items), mc = rep(0L, n_items),
              thresholds = NULL)
  if (!is_poly) {
    pbar <- colMeans(Y, na.rm = TRUE)
    par$d <- -qlogis(pmin(pmax(pbar, 0.02), 0.98))
    if (model == "3PL") {
      par$mc <- if (is.null(mc_items)) rep(1L, n_items) else {
        m <- rep(0L, n_items)
        if (is.logical(mc_items)) m <- as.integer(mc_items) else m[mc_items] <- 1L
        m
      }
      par$c <- ifelse(par$mc == 1, 0.15, 0)
    }
  } else {
    par$thresholds <- lapply(seq_len(n_items), function(j) {
      K <- K_per_item[j]
      v <- Y[obs_mask[, j], j]
      cum_p <- cumsum(table(factor(v, levels = 1:K)))[-K] / length(v)
      qlogis(pmin(pmax(cum_p, 0.05), 0.95))
    })
  }

  # ---- Per-item M-step objective (expected complete-data loglik) ----
  # Dichotomous sufficient stats: r_j (Q), n_j (Q)
  mstep_dichot <- function(j, r_j, n_j) {
    common_a <- free_sigma
    f <- function(th) {
      a <- if (common_a) par$a[j] else exp(th[1])
      d <- if (common_a) th[1] else th[2]
      p <- plogis(a * (nodes - d))
      if (model == "3PL" && par$mc[j] == 1) {
        cc <- plogis(th[length(th)])
        p <- cc + (1 - cc) * p
      }
      p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
      -sum(r_j * log(p) + (n_j - r_j) * log(1 - p))
    }
    if (common_a) {
      th0 <- par$d[j]
    } else if (model == "3PL" && par$mc[j] == 1) {
      th0 <- c(log(par$a[j]), par$d[j], qlogis(max(par$c[j], 0.02)))
    } else {
      th0 <- c(log(par$a[j]), par$d[j])
    }
    o <- optim(th0, f, method = "BFGS",
               control = list(maxit = 50, reltol = 1e-10))
    th <- o$par
    if (common_a) {
      par$d[j] <<- th[1]
    } else {
      par$a[j] <<- min(exp(th[1]), 10)
      par$d[j] <<- th[2]
      if (model == "3PL" && par$mc[j] == 1) {
        par$c[j] <<- plogis(th[length(th)])
      }
    }
  }

  # Polytomous sufficient stats: R_j (Q x K expected counts)
  mstep_poly <- function(j, R_j) {
    K <- K_per_item[j]
    common_a <- free_sigma                     # PCM
    uses_a <- model %in% c("GRM", "GPCM", "NRM")
    # unconstrained parameterization of thresholds:
    # GRM: tau_1, log-spacings; PCM/GPCM/NRM: free values
    to_thresh <- function(v) {
      if (model == "GRM" && K > 2) {
        cumsum(c(v[1], exp(v[-1])))
      } else v
    }
    from_thresh <- function(tr) {
      if (model == "GRM" && K > 2) {
        c(tr[1], log(pmax(diff(tr), 1e-3)))
      } else tr
    }
    f <- function(th) {
      a <- if (common_a) par$a[j]
           else if (uses_a) exp(th[K])
           else 1
      tr <- to_thresh(th[seq_len(K - 1)])
      P <- pmax(irt_category_probs(prob_model, nodes, tr, a), 1e-12)
      -sum(R_j[, seq_len(K)] * log(P))
    }
    th0 <- from_thresh(par$thresholds[[j]])
    if (uses_a && !common_a) th0 <- c(th0, log(par$a[j]))
    o <- optim(th0, f, method = "BFGS",
               control = list(maxit = 60, reltol = 1e-10))
    par$thresholds[[j]] <<- to_thresh(o$par[seq_len(K - 1)])
    if (uses_a && !common_a) par$a[j] <<- min(exp(o$par[K]), 10)
  }

  # ---- Common-slope update (Rasch/PCM ability SD) ----
  mstep_common_slope <- function(stats) {
    f <- function(log_a) {
      a <- exp(log_a)
      total <- 0
      if (!is_poly) {
        for (j in seq_len(n_items)) {
          p <- pmin(pmax(plogis(a * (nodes - par$d[j])), 1e-10), 1 - 1e-10)
          total <- total - sum(stats$r[j, ] * log(p) +
                                 (stats$n[j, ] - stats$r[j, ]) * log(1 - p))
        }
      } else {
        for (j in seq_len(n_items)) {
          P <- pmax(irt_category_probs(prob_model, nodes,
                                       par$thresholds[[j]], a), 1e-12)
          total <- total - sum(stats$R[[j]] * log(P))
        }
      }
      total
    }
    o <- optimize(f, interval = log(c(0.05, 10)))
    par$a[] <<- exp(o$minimum)
  }

  # ---- EM loop ----
  loglik_old <- -Inf
  loglik <- NA_real_
  converged <- FALSE
  iter <- 0

  for (iter in seq_len(max_iter)) {
    P <- item_probs_at_nodes(par)
    L <- loglik_matrix(P)                       # N x Q

    # E-step: posterior over nodes (log-sum-exp by row)
    Lw <- sweep(L, 2, log_A, "+")
    m_row <- apply(Lw, 1, max)
    W <- exp(Lw - m_row)
    row_sums <- rowSums(W)
    loglik <- sum(w_person * (m_row + log(row_sums)))
    W <- (W / row_sums) * w_person              # weighted posteriors

    if (abs(loglik - loglik_old) < tol) {   # absolute, mirt-style
      converged <- TRUE
      break
    }
    loglik_old <- loglik

    # Expected counts and M-step
    if (!is_poly) {
      Y0 <- ifelse(obs_mask, Y, 0)
      r <- t(Y0) %*% W                          # J x Q expected correct
      n_exp <- t(obs_mask * 1) %*% W            # J x Q expected exposure
      for (j in seq_len(n_items)) mstep_dichot(j, r[j, ], n_exp[j, ])
      if (free_sigma) mstep_common_slope(list(r = r, n = n_exp))
    } else {
      R <- lapply(seq_len(n_items), function(j) {
        K <- K_per_item[j]
        Rj <- matrix(0, Q, K)
        obs_j <- obs_mask[, j]
        for (k in seq_len(K)) {
          sel <- obs_j & Y[, j] == k
          if (any(sel)) Rj[, k] <- colSums(W[sel, , drop = FALSE])
        }
        Rj
      })
      for (j in seq_len(n_items)) mstep_poly(j, R[[j]])
      if (free_sigma) mstep_common_slope(list(R = R))
    }
  }

  # ---- Output on the external scale ----
  sigma_hat <- if (free_sigma) par$a[1] else 1

  # EAP abilities (posterior means), back-transformed for free-sigma models
  P <- item_probs_at_nodes(par)
  L <- loglik_matrix(P)
  Lw <- sweep(L, 2, log_A, "+")
  m_row <- apply(Lw, 1, max)
  W <- exp(Lw - m_row); W <- W / rowSums(W)
  theta_hat <- as.numeric(W %*% nodes) * sigma_hat
  names(theta_hat) <- paste0("Person", seq_len(n_persons))

  if (!is_poly) {
    difficulty <- if (free_sigma) par$d * sigma_hat else par$d
    discrimination <- if (free_sigma) rep(1, n_items) else par$a
    names(difficulty) <- names(discrimination) <- paste0("Item", 1:n_items)
    item_parameters <- data.frame(
      difficulty = difficulty,
      discrimination = discrimination,
      guessing = if (model == "3PL") ifelse(par$mc == 1, par$c, NA) else NA
    )
    n_par <- n_items * switch(model, Rasch = 1, "2PL" = 2, "3PL" = 2) +
      (if (model == "3PL") sum(par$mc) else 0) + as.integer(free_sigma)
  } else {
    thresholds <- if (free_sigma) {
      lapply(par$thresholds, function(tr) tr * sigma_hat)
    } else par$thresholds
    names(thresholds) <- paste0("Item", 1:n_items)
    discrimination <- if (free_sigma) rep(1, n_items) else par$a
    names(discrimination) <- paste0("Item", 1:n_items)
    item_parameters <- list(
      thresholds = thresholds,
      discrimination = discrimination
    )
    n_par <- sum(K_per_item - 1) +
      (if (model %in% c("GRM", "GPCM", "NRM")) n_items else 0) +
      as.integer(free_sigma)
  }

  result <- list(
    model = model,
    method = "EM",
    item_parameters = item_parameters,
    person_abilities = theta_hat,
    ability_sd = sigma_hat,
    mc_items = if (model == "3PL") par$mc else NULL,
    logLik = loglik,
    AIC = -2 * loglik + 2 * n_par,
    BIC = -2 * loglik + log(n_persons) * n_par,
    convergence = list(
      converged = converged,
      message = if (converged) "EM converged" else "Maximum EM iterations reached",
      iterations = iter
    ),
    n_persons = n_persons,
    n_items = n_items,
    quad_points = quad_points
  )
  if (is_poly) {
    result$n_categories_per_item <- K_per_item
    result$max_categories <- max_K
  }

  class(result) <- if (is_poly) {
    c("gllamm_irt_poly", "gllamm_irt", "gllamm")
  } else {
    c("gllamm_irt", "gllamm")
  }
  result
}
