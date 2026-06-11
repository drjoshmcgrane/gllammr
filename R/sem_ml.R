#' Covariance-based ML estimation for SEM (Wishart / FIML likelihood)
#'
#' Fits the measurement + recursive structural model. With complete data
#' the data reduce to the sample covariance matrix and mean vector
#' (Wishart ML, the lavaan/LISREL approach - fitting cost independent of
#' N); with \code{missing = "fiml"} the casewise (missing-pattern) normal
#' likelihood is maximized directly, with means as free parameters.
#'
#' Model structure: y = nu + Lambda eta + epsilon; eta = B eta + zeta.
#' Exogenous latent variables (no incoming B paths - including the
#' pseudo-latents that carry observed structural covariates) have a free
#' covariance matrix parameterized by its Cholesky factor; endogenous
#' disturbances are uncorrelated.
#'
#' @param Y Observed-variable matrix (indicators first, then any observed
#'   structural covariates appended by the caller)
#' @param lambda_pattern p x q loading pattern (0 zero / 1 free / 2 fixed 1)
#' @param beta_pattern q x q structural pattern (1 = free path row ~ col)
#' @param theta_zero Logical p-vector: TRUE for rows whose residual
#'   variance is fixed at 0 (covariate pseudo-indicators)
#' @param missing "listwise" or "fiml"
#' @param se Compute standard errors (numerical Hessian; default TRUE)
#' @param control Optimization control list
#'
#' @keywords internal
fit_sem_ml <- function(Y, lambda_pattern, beta_pattern, theta_zero = NULL,
                       missing = c("listwise", "fiml"), se = TRUE,
                       control = list()) {
  missing <- match.arg(missing)
  p <- ncol(Y)
  q <- ncol(lambda_pattern)
  if (is.null(theta_zero)) theta_zero <- rep(FALSE, p)

  has_na <- anyNA(Y)
  if (has_na && missing == "listwise") {
    keep <- stats::complete.cases(Y)
    Y <- Y[keep, , drop = FALSE]
    has_na <- FALSE
  }
  n_obs <- nrow(Y)
  fiml <- missing == "fiml" && has_na

  # ---- Structure bookkeeping ----
  lambda_names <- rownames(lambda_pattern)
  latent_names <- colnames(lambda_pattern)
  exo <- which(rowSums(beta_pattern) == 0)        # no incoming paths
  endo <- setdiff(seq_len(q), exo)
  n_exo <- length(exo)
  n_endo <- length(endo)
  free_theta <- which(!theta_zero)
  n_theta_free <- length(free_theta)

  n_lambda <- sum(lambda_pattern == 1L)
  n_b <- sum(beta_pattern == 1L)
  n_chol <- n_exo * (n_exo + 1) / 2
  n_cov_params <- n_lambda + n_b + n_chol + n_endo + n_theta_free
  n_par <- n_cov_params + if (fiml) p else 0L

  I_q <- diag(q)

  # th layout: [lambda_free | b_free | exo chol (diag log, lower free) |
  #             log psi_endo | log theta_free | (fiml) mu]
  implied <- function(th) {
    k <- 0
    Lambda <- matrix(0, p, q)
    Lambda[lambda_pattern == 2L] <- 1
    if (n_lambda) {
      Lambda[lambda_pattern == 1L] <- th[k + seq_len(n_lambda)]
      k <- k + n_lambda
    }
    B <- matrix(0, q, q)
    if (n_b) {
      B[beta_pattern == 1L] <- th[k + seq_len(n_b)]
      k <- k + n_b
    }
    Psi <- matrix(0, q, q)
    if (n_exo) {
      L <- matrix(0, n_exo, n_exo)
      diag(L) <- exp(th[k + seq_len(n_exo)])
      k <- k + n_exo
      if (n_exo > 1) {
        L[lower.tri(L)] <- th[k + seq_len(n_exo * (n_exo - 1) / 2)]
        k <- k + n_exo * (n_exo - 1) / 2
      }
      Psi[exo, exo] <- L %*% t(L)
    }
    if (n_endo) {
      diag(Psi)[endo] <- exp(th[k + seq_len(n_endo)])^2
      k <- k + n_endo
    }
    theta_var <- rep(0, p)
    if (n_theta_free) {
      theta_var[free_theta] <- exp(th[k + seq_len(n_theta_free)])^2
      k <- k + n_theta_free
    }
    mu <- if (fiml) th[k + seq_len(p)] else NULL

    IB_inv <- solve(I_q - B)
    V_eta <- IB_inv %*% Psi %*% t(IB_inv)
    Sigma <- Lambda %*% V_eta %*% t(Lambda) + diag(theta_var, p)
    list(Lambda = Lambda, B = B, Psi = Psi, theta_var = theta_var,
         V_eta = V_eta, Sigma = Sigma, mu = mu)
  }

  # ---- Objectives ----
  if (!fiml) {
    mu_hat <- colMeans(Y)
    S <- crossprod(sweep(Y, 2, mu_hat)) / n_obs
    logdet_S <- as.numeric(determinant(S, logarithm = TRUE)$modulus)

    nll <- function(th) {
      Sg <- implied(th)$Sigma
      ld <- determinant(Sg, logarithm = TRUE)
      if (ld$sign <= 0) return(1e10)
      tr <- sum(diag(solve(Sg, S)))
      if (!is.finite(tr)) return(1e10)
      0.5 * n_obs * (p * log(2 * pi) + as.numeric(ld$modulus) + tr)
    }
  } else {
    # Pattern-based FIML
    pat_key <- apply(!is.na(Y), 1, paste, collapse = "")
    pat_levels <- unique(pat_key)
    pat_rows <- split(seq_len(n_obs), factor(pat_key, levels = pat_levels))
    pat_obs <- lapply(pat_rows, function(r) which(!is.na(Y[r[1], ])))
    if (any(vapply(pat_obs, length, 0L) == 0)) {
      stop("Some rows have no observed variables; remove them first")
    }

    nll <- function(th) {
      mom <- implied(th)
      total <- 0
      for (g in seq_along(pat_rows)) {
        ob <- pat_obs[[g]]
        rows <- pat_rows[[g]]
        Sg <- mom$Sigma[ob, ob, drop = FALSE]
        ch <- tryCatch(chol(Sg), error = function(e) NULL)
        if (is.null(ch)) return(1e10)
        ctr <- sweep(Y[rows, ob, drop = FALSE], 2, mom$mu[ob])
        z <- forwardsolve(t(ch), t(ctr))
        total <- total + length(rows) *
          (length(ob) * log(2 * pi) + 2 * sum(log(diag(ch)))) + sum(z^2)
      }
      0.5 * total
    }
  }

  # ---- Starting values ----
  sd0 <- apply(Y, 2, stats::sd, na.rm = TRUE)
  sd0[!is.finite(sd0) | sd0 <= 0] <- 1
  th0 <- c(rep(1, n_lambda), rep(0, n_b),
           rep(log(stats::median(sd0)), n_exo),
           rep(0, n_exo * (n_exo - 1) / 2),
           rep(log(stats::median(sd0)), n_endo),
           log(pmax(sd0[free_theta] / 2, 1e-3)))
  if (fiml) th0 <- c(th0, colMeans(Y, na.rm = TRUE))

  control_defaults <- list(eval.max = 10000, iter.max = 5000, trace = 0)
  ctl <- modifyList(control_defaults, control)
  opt <- nlminb(th0, nll, control = ctl)

  est <- implied(opt$par)
  loglik <- -opt$objective
  mu_out <- if (fiml) est$mu else colMeans(Y)

  # ---- Saturated and baseline log-likelihoods, fit indices ----
  if (!fiml) {
    loglik_sat <- -0.5 * n_obs * (p * log(2 * pi) + logdet_S + p)
    loglik_base <- -0.5 * n_obs *
      (p * log(2 * pi) + sum(log(diag(S))) + p)
  } else {
    sat <- .mvn_em_saturated(Y)
    loglik_sat <- sat$loglik
    loglik_base <- sum(vapply(seq_len(p), function(j) {
      yj <- Y[!is.na(Y[, j]), j]
      v <- mean((yj - mean(yj))^2)
      sum(stats::dnorm(yj, mean(yj), sqrt(v), log = TRUE))
    }, 0))
  }

  chisq <- max(2 * (loglik_sat - loglik), 0)
  df <- p * (p + 1) / 2 - n_cov_params
  chisq_base <- max(2 * (loglik_sat - loglik_base), 0)
  df_base <- p * (p - 1) / 2

  cfi <- 1 - max(chisq - df, 0) / max(chisq - df, chisq_base - df_base, 1e-10)
  tli <- if (df > 0 && df_base > 0) {
    (chisq_base / df_base - chisq / df) / (chisq_base / df_base - 1)
  } else NA_real_
  rmsea <- if (df > 0) sqrt(max(chisq - df, 0) / (df * n_obs)) else NA_real_
  rmsea_ci <- if (df > 0) .rmsea_ci(chisq, df, n_obs) else c(NA, NA)

  srmr <- {
    Sg <- est$Sigma
    Sd <- if (!fiml) S else .mvn_em_saturated(Y)$Sigma
    sd_d <- sqrt(diag(Sd))
    R <- (Sd - Sg) / (sd_d %o% sd_d)
    sqrt(mean(R[lower.tri(R, diag = TRUE)]^2))
  }

  fit_measures <- c(chisq = chisq, df = df,
                    pvalue = if (df > 0) stats::pchisq(chisq, df,
                                                       lower.tail = FALSE)
                             else NA_real_,
                    cfi = cfi, tli = tli, rmsea = rmsea,
                    rmsea_ci_lower = rmsea_ci[1],
                    rmsea_ci_upper = rmsea_ci[2], srmr = srmr)

  # ---- Natural-parameter table with delta-method SEs ----
  nat <- function(th) {
    m <- implied(th)
    out <- numeric(0)
    if (n_lambda) out <- c(out, m$Lambda[lambda_pattern == 1L])
    if (n_b) out <- c(out, m$B[beta_pattern == 1L])
    if (n_exo) {
      Pe <- m$Psi[exo, exo, drop = FALSE]
      out <- c(out, Pe[lower.tri(Pe, diag = TRUE)])
    }
    if (n_endo) out <- c(out, diag(m$Psi)[endo])
    if (n_theta_free) out <- c(out, m$theta_var[free_theta])
    out
  }
  nat_names <- character(0)
  if (n_lambda) {
    idx <- which(lambda_pattern == 1L, arr.ind = TRUE)
    nat_names <- c(nat_names, paste0(latent_names[idx[, 2]], "=~",
                                     lambda_names[idx[, 1]]))
  }
  if (n_b) {
    idx <- which(beta_pattern == 1L, arr.ind = TRUE)
    nat_names <- c(nat_names, paste0(latent_names[idx[, 1]], "~",
                                     latent_names[idx[, 2]]))
  }
  if (n_exo) {
    en <- latent_names[exo]
    idx <- which(lower.tri(diag(n_exo), diag = TRUE), arr.ind = TRUE)
    nat_names <- c(nat_names, paste0(en[idx[, 2]], "~~", en[idx[, 1]]))
  }
  if (n_endo) {
    nat_names <- c(nat_names, paste0(latent_names[endo], "~~",
                                     latent_names[endo]))
  }
  if (n_theta_free) {
    nat_names <- c(nat_names, paste0(lambda_names[free_theta], "~~",
                                     lambda_names[free_theta]))
  }

  est_nat <- nat(opt$par)
  se_nat <- rep(NA_real_, length(est_nat))
  vcov_nat <- NULL
  if (se) {
    H <- tryCatch(stats::optimHess(opt$par, nll), error = function(e) NULL)
    V <- if (!is.null(H)) tryCatch(solve(H), error = function(e) NULL)
         else NULL
    if (!is.null(V)) {
      # Numerical Jacobian of the natural parameters
      m0 <- est_nat
      J <- matrix(0, length(m0), length(opt$par))
      h <- pmax(abs(opt$par), 1) * 1e-5
      for (j in seq_along(opt$par)) {
        tp <- opt$par; tp[j] <- tp[j] + h[j]
        tm <- opt$par; tm[j] <- tm[j] - h[j]
        J[, j] <- (nat(tp) - nat(tm)) / (2 * h[j])
      }
      vcov_nat <- J %*% V[seq_len(ncol(J)), seq_len(ncol(J))] %*% t(J)
      d <- diag(vcov_nat)
      se_nat <- ifelse(d > 0, sqrt(d), NA_real_)
      dimnames(vcov_nat) <- list(nat_names, nat_names)
    }
  }
  param_table <- data.frame(
    label = nat_names, est = est_nat, se = se_nat,
    z = est_nat / se_nat,
    pvalue = 2 * stats::pnorm(-abs(est_nat / se_nat)),
    stringsAsFactors = FALSE)

  # ---- Factor scores (regression method; casewise under FIML) ----
  CrossCov <- est$V_eta %*% t(est$Lambda)        # q x p: cov(eta, y)
  if (!fiml) {
    scores <- sweep(Y, 2, mu_out) %*% t(CrossCov %*% solve(est$Sigma))
  } else {
    scores <- matrix(NA_real_, n_obs, q)
    for (g in seq_along(pat_rows)) {
      ob <- pat_obs[[g]]
      rows <- pat_rows[[g]]
      W <- CrossCov[, ob, drop = FALSE] %*%
        solve(est$Sigma[ob, ob, drop = FALSE])
      scores[rows, ] <- sweep(Y[rows, ob, drop = FALSE], 2,
                              mu_out[ob]) %*% t(W)
    }
  }

  list(Lambda = est$Lambda, B = est$B, Psi = est$Psi,
       theta_var = est$theta_var, V_eta = est$V_eta, Sigma = est$Sigma,
       intercepts = mu_out, factor_scores = scores,
       logLik = loglik, logLik_saturated = loglik_sat,
       logLik_baseline = loglik_base,
       fit_measures = fit_measures,
       param_table = param_table, vcov = vcov_nat,
       exo = exo, endo = endo,
       missing = if (fiml) "fiml" else "listwise",
       n_obs = n_obs,
       converged = (opt$convergence == 0), message = opt$message,
       n_params = n_par + if (fiml) 0L else p)   # + p saturated means
}


#' EM for the saturated multivariate-normal model under missingness
#'
#' Standard EM on the expected sufficient statistics: the E-step fills in
#' conditional means (and adds the conditional covariance) per missing
#' pattern; the M-step is the sample mean / covariance of the completed
#' statistics. Used for the FIML saturated log-likelihood (fit indices).
#'
#' @keywords internal
.mvn_em_saturated <- function(Y, max_iter = 500, tol = 1e-8) {
  n <- nrow(Y); p <- ncol(Y)
  pat_key <- apply(!is.na(Y), 1, paste, collapse = "")
  pat_levels <- unique(pat_key)
  pat_rows <- split(seq_len(n), factor(pat_key, levels = pat_levels))
  pat_obs <- lapply(pat_rows, function(r) which(!is.na(Y[r[1], ])))

  mu <- colMeans(Y, na.rm = TRUE)
  Sg <- diag(apply(Y, 2, stats::var, na.rm = TRUE), p)
  ll_old <- -Inf
  ll <- NA_real_

  for (it in seq_len(max_iter)) {
    T1 <- numeric(p)
    T2 <- matrix(0, p, p)
    ll <- 0
    for (g in seq_along(pat_rows)) {
      ob <- pat_obs[[g]]
      mi <- setdiff(seq_len(p), ob)
      rows <- pat_rows[[g]]
      Yo <- Y[rows, ob, drop = FALSE]
      So <- Sg[ob, ob, drop = FALSE]
      ch <- chol(So)
      ctr <- sweep(Yo, 2, mu[ob])
      z <- forwardsolve(t(ch), t(ctr))
      ll <- ll - 0.5 * (length(rows) *
        (length(ob) * log(2 * pi) + 2 * sum(log(diag(ch)))) + sum(z^2))

      Yc <- matrix(0, length(rows), p)
      Yc[, ob] <- Yo
      Cadd <- matrix(0, p, p)
      if (length(mi)) {
        W <- Sg[mi, ob, drop = FALSE] %*% chol2inv(ch)
        Yc[, mi] <- matrix(mu[mi], length(rows), length(mi),
                           byrow = TRUE) + ctr %*% t(W)
        Cmiss <- Sg[mi, mi, drop = FALSE] -
          W %*% Sg[ob, mi, drop = FALSE]
        Cadd[mi, mi] <- length(rows) * Cmiss
      }
      T1 <- T1 + colSums(Yc)
      T2 <- T2 + crossprod(Yc) + Cadd
    }
    mu <- T1 / n
    Sg <- T2 / n - tcrossprod(mu)
    if (abs(ll - ll_old) < tol) break
    ll_old <- ll
  }
  list(mu = mu, Sigma = Sg, loglik = ll)
}


#' 90% RMSEA confidence interval (noncentral chi-square inversion)
#' @keywords internal
.rmsea_ci <- function(chisq, df, n, level = 0.90) {
  alpha <- (1 - level) / 2
  ncp_for <- function(target) {
    if (stats::pchisq(chisq, df, ncp = 0) < target) return(0)
    upper <- max(chisq * 2, df * 10, 10)
    while (stats::pchisq(chisq, df, ncp = upper) > target) upper <- upper * 2
    stats::uniroot(function(nc) stats::pchisq(chisq, df, ncp = nc) - target,
                   c(0, upper))$root
  }
  lo <- tryCatch(ncp_for(1 - alpha), error = function(e) 0)
  hi <- tryCatch(ncp_for(alpha), error = function(e) NA_real_)
  c(sqrt(max(lo, 0) / (df * n)), sqrt(max(hi, 0) / (df * n)))
}
