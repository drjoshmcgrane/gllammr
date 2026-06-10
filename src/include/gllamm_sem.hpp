// Structural Equation Model (SEM) with latent variables
// Measurement model with continuous indicators + recursive structural model.
//
// lambda_pattern codes: 0 = loading fixed at zero, 1 = free loading,
// 2 = loading fixed at one (marker-variable identification).
// Variances are log-parameterized for unconstrained optimization.

#ifndef GLLAMM_SEM_HPP
#define GLLAMM_SEM_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_sem(objective_function<Type>* obj)
{
  // Data inputs
  DATA_MATRIX(Y);               // Observed indicators matrix (n x p)
  DATA_INTEGER(n_obs);          // Number of observations
  DATA_INTEGER(n_indicators);   // Number of observed indicators
  DATA_INTEGER(n_latent);       // Number of latent variables
  DATA_IMATRIX(lambda_pattern); // (p x q): 0 = zero, 1 = free, 2 = fixed at 1
  DATA_IMATRIX(beta_pattern);   // Structural paths (q x q): 1 = free path j <- k

  // Parameters
  PARAMETER_VECTOR(nu);          // Indicator intercepts (p)
  PARAMETER_VECTOR(lambda_free); // Free factor loadings (concatenated)
  PARAMETER_VECTOR(beta_free);   // Free structural coefficients (concatenated)
  PARAMETER_VECTOR(log_psi);     // Log latent residual SDs (q)
  PARAMETER_VECTOR(log_theta);   // Log indicator residual SDs (p)
  PARAMETER_MATRIX(eta);         // Latent variable scores (n x q)

  // Build loading matrix Lambda (p x q)
  matrix<Type> Lambda(n_indicators, n_latent);
  Lambda.setZero();
  int lambda_idx = 0;
  for (int i = 0; i < n_indicators; i++) {
    for (int j = 0; j < n_latent; j++) {
      if (lambda_pattern(i, j) == 1) {
        Lambda(i, j) = lambda_free(lambda_idx);
        lambda_idx++;
      } else if (lambda_pattern(i, j) == 2) {
        Lambda(i, j) = Type(1.0);   // marker variable
      }
    }
  }

  // Build structural coefficient matrix Beta (q x q), recursive
  matrix<Type> Beta(n_latent, n_latent);
  Beta.setZero();
  int beta_idx = 0;
  for (int i = 0; i < n_latent; i++) {
    for (int j = 0; j < n_latent; j++) {
      if (beta_pattern(i, j) == 1) {
        Beta(i, j) = beta_free(beta_idx);
        beta_idx++;
      }
    }
  }

  vector<Type> psi_sd = exp(log_psi.array());
  vector<Type> theta_sd = exp(log_theta.array());

  Type nll = 0.0;

  for (int i = 0; i < n_obs; i++) {
    vector<Type> eta_i = eta.row(i);

    // Structural model (recursive): eta_j ~ N(sum_k Beta_jk eta_k, psi_j^2)
    for (int j = 0; j < n_latent; j++) {
      Type eta_expected = 0.0;
      for (int k = 0; k < n_latent; k++) {
        if (k != j) {
          eta_expected += Beta(j, k) * eta_i(k);
        }
      }
      nll -= dnorm(eta_i(j), eta_expected, psi_sd(j), true);
    }

    // Measurement model: Y = nu + Lambda * eta + epsilon
    for (int p = 0; p < n_indicators; p++) {
      Type mu_ip = nu(p);
      for (int q = 0; q < n_latent; q++) {
        mu_ip += Lambda(p, q) * eta_i(q);
      }
      nll -= dnorm(Y(i, p), mu_ip, theta_sd(p), true);
    }
  }

  ADREPORT(Lambda);
  ADREPORT(Beta);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_SEM_HPP
