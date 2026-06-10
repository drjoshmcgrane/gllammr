// Gaussian GLMM with random intercepts
// TMB template for basic 2-level model

#ifndef GLLAMM_GAUSSIAN_HPP
#define GLLAMM_GAUSSIAN_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_gaussian(objective_function<Type>* obj)
{
  // Data inputs
  DATA_VECTOR(y);              // Response vector
  DATA_MATRIX(X);              // Fixed effects design matrix
  DATA_MATRIX(Z);              // Random effects design matrix
  DATA_IVECTOR(groups);        // Group indices (0-indexed)
  DATA_INTEGER(n_groups);      // Number of groups
  DATA_INTEGER(n_obs);         // Number of observations
  DATA_INTEGER(n_fixed);       // Number of fixed effects
  DATA_INTEGER(n_random);      // Number of random effects per group
  DATA_VECTOR(weights);        // Case weights (fweights or pweights)

  // Parameters
  PARAMETER_VECTOR(beta);      // Fixed effects coefficients
  PARAMETER_VECTOR(u);         // Random effects (all groups concatenated)
  PARAMETER(log_sigma);        // Log residual standard deviation
  PARAMETER(log_sigma_u);      // Log random effects standard deviation

  // Transform parameters
  Type sigma = exp(log_sigma);
  Type sigma_u = exp(log_sigma_u);

  // Initialize negative log-likelihood
  // parallel_accumulator splits the likelihood across OpenMP threads
  // when available (no-op on single-threaded builds)
  parallel_accumulator<Type> nll(obj);

  // Prior for random effects: u ~ N(0, sigma_u^2)
  for (int j = 0; j < n_groups; j++) {
    for (int k = 0; k < n_random; k++) {
      int idx = j * n_random + k;
      nll -= dnorm(u[idx], Type(0.0), sigma_u, true);
    }
  }

  // Likelihood for observations
  for (int i = 0; i < n_obs; i++) {
    // Linear predictor: fixed effects + random effects
    Type eta = 0.0;

    // Add fixed effects
    for (int j = 0; j < n_fixed; j++) {
      eta += X(i, j) * beta[j];
    }

    // Add random effects
    int g = groups[i];  // Group for this observation
    for (int k = 0; k < n_random; k++) {
      int u_idx = g * n_random + k;
      eta += Z(i, k) * u[u_idx];
    }

    // Gaussian likelihood: y ~ N(eta, sigma^2) (weighted)
    Type w_i = weights[i];
    nll -= w_i * dnorm(y[i], eta, sigma, true);
  }

  // Report fitted values
  vector<Type> fitted(n_obs);
  for (int i = 0; i < n_obs; i++) {
    Type eta = 0.0;
    for (int j = 0; j < n_fixed; j++) {
      eta += X(i, j) * beta[j];
    }
    int g = groups[i];
    for (int k = 0; k < n_random; k++) {
      int u_idx = g * n_random + k;
      eta += Z(i, k) * u[u_idx];
    }
    fitted[i] = eta;
  }

  REPORT(fitted);
  ADREPORT(sigma);
  ADREPORT(sigma_u);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_GAUSSIAN_HPP
