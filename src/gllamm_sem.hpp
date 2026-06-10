// Structural Equation Model (SEM) with latent variables
// Basic implementation with measurement model and structural model

#ifndef GLLAMM_SEM_HPP
#define GLLAMM_SEM_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_sem(objective_function<Type>* obj)
{
  // Data inputs
  DATA_MATRIX(Y);              // Observed indicators matrix (n x p)
  DATA_INTEGER(n_obs);         // Number of observations
  DATA_INTEGER(n_indicators);  // Number of observed indicators
  DATA_INTEGER(n_latent);      // Number of latent variables
  DATA_IMATRIX(lambda_pattern); // Loading pattern matrix (p x q): 1=free, 0=fixed
  DATA_IMATRIX(beta_pattern);  // Structural paths pattern (q x q): 1=free, 0=fixed

  // Parameters
  PARAMETER_VECTOR(lambda_free); // Free factor loadings (concatenated)
  PARAMETER_VECTOR(beta_free);   // Free structural coefficients (concatenated)
  PARAMETER_VECTOR(psi);         // Latent variable residual variances
  PARAMETER_VECTOR(theta);       // Indicator residual variances
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
      }
    }
  }

  // Build structural coefficient matrix Beta (q x q)
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

  // Initialize negative log-likelihood
  Type nll = 0.0;

  // For each observation
  for (int i = 0; i < n_obs; i++) {
    vector<Type> eta_i(n_latent);
    for (int j = 0; j < n_latent; j++) {
      eta_i(j) = eta(i, j);
    }

    // Structural model: eta = Beta * eta + zeta
    // Solve for eta given structural model
    // For simplicity, assume recursive structure (no simultaneity)
    // eta = (I - Beta)^{-1} * zeta

    // Prior for latent variables (given structural model)
    // Simplified: assume diagonal Psi for now
    for (int j = 0; j < n_latent; j++) {
      Type eta_expected = 0.0;
      for (int k = 0; k < n_latent; k++) {
        if (k != j) {
          eta_expected += Beta(j, k) * eta_i(k);
        }
      }
      nll -= dnorm(eta_i(j), eta_expected, sqrt(psi(j)), true);
    }

    // Measurement model: Y = Lambda * eta + epsilon
    vector<Type> y_i(n_indicators);
    vector<Type> mu_i(n_indicators);

    for (int p = 0; p < n_indicators; p++) {
      y_i(p) = Y(i, p);
      mu_i(p) = 0.0;
      for (int q = 0; q < n_latent; q++) {
        mu_i(p) += Lambda(p, q) * eta_i(q);
      }
    }

    // Likelihood for indicators
    for (int p = 0; p < n_indicators; p++) {
      nll -= dnorm(y_i(p), mu_i(p), sqrt(theta(p)), true);
    }
  }

  // Report
  ADREPORT(Lambda);
  ADREPORT(Beta);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_SEM_HPP
