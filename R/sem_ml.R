#' Covariance-based ML estimation for SEM (Wishart likelihood)
#'
#' Fits the measurement + recursive structural model from the sample
#' covariance matrix by minimizing the ML discrepancy
#' F = log|Sigma(theta)| + tr(S Sigma^-1) - log|S| - p,
#' the lavaan/LISREL approach. The data reduce to S and the mean vector, so
#' fitting cost is independent of N. For complete data this yields the same
#' ML estimates as the full-data marginal likelihood.
#'
#' @keywords internal
fit_sem_ml <- function(Y, lambda_pattern, beta_pattern, control = list()) {
  n_obs <- nrow(Y)
  p <- ncol(Y)
  q <- ncol(lambda_pattern)

  mu_hat <- colMeans(Y)
  S <- crossprod(sweep(Y, 2, mu_hat)) / n_obs        # ML covariance (divisor N)
  logdet_S <- determinant(S, logarithm = TRUE)$modulus

  n_lambda <- sum(lambda_pattern == 1L)
  n_beta <- sum(beta_pattern == 1L)
  I_q <- diag(q)

  implied_sigma <- function(th) {
    Lambda <- matrix(0, p, q)
    Lambda[lambda_pattern == 2L] <- 1
    if (n_lambda) Lambda[lambda_pattern == 1L] <- th[seq_len(n_lambda)]
    k <- n_lambda
    B <- matrix(0, q, q)
    if (n_beta) {
      B[beta_pattern == 1L] <- th[k + seq_len(n_beta)]
      k <- k + n_beta
    }
    psi <- exp(th[k + seq_len(q)])                   # latent residual SDs
    k <- k + q
    theta_sd <- exp(th[k + seq_len(p)])              # indicator residual SDs

    IB_inv <- solve(I_q - B)
    V_eta <- IB_inv %*% diag(psi^2, q) %*% t(IB_inv)
    list(Sigma = Lambda %*% V_eta %*% t(Lambda) + diag(theta_sd^2, p),
         Lambda = Lambda, B = B, psi = psi, theta_sd = theta_sd,
         V_eta = V_eta)
  }

  discrepancy <- function(th) {
    Sg <- implied_sigma(th)$Sigma
    ld <- determinant(Sg, logarithm = TRUE)
    if (ld$sign <= 0) return(1e10)
    as.numeric(ld$modulus + sum(diag(solve(Sg, S))) - logdet_S - p)
  }

  th0 <- c(rep(1, n_lambda), rep(0, n_beta),
           rep(log(stats::sd(Y[, 1])), q),
           log(apply(Y, 2, stats::sd) / 2))

  control_defaults <- list(eval.max = 5000, iter.max = 2000, trace = 0)
  ctl <- modifyList(control_defaults, control)
  opt <- nlminb(th0, discrepancy, control = ctl)

  est <- implied_sigma(opt$par)
  # Full-data Gaussian log-likelihood at the ML solution
  loglik <- -0.5 * n_obs *
    (p * log(2 * pi) +
       as.numeric(determinant(est$Sigma, logarithm = TRUE)$modulus) +
       sum(diag(solve(est$Sigma, S))))

  # Regression-method factor scores: E[eta | y] = V_eta Lambda' Sigma^-1 (y - mu)
  scores <- sweep(Y, 2, mu_hat) %*%
    t(est$V_eta %*% t(est$Lambda) %*% solve(est$Sigma))

  list(Lambda = est$Lambda, B = est$B, psi = est$psi,
       theta_sd = est$theta_sd, intercepts = mu_hat,
       factor_scores = scores, logLik = loglik,
       discrepancy = opt$objective,
       converged = (opt$convergence == 0), message = opt$message,
       n_params = n_lambda + n_beta + q + p + p)   # + p intercepts (saturated)
}
