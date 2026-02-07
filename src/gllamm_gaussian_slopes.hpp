// Gaussian GLMM with random intercepts and slopes
// Supports multiple random effects per group with variance-covariance structure

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data inputs
  DATA_VECTOR(y);              // Response vector
  DATA_MATRIX(X);              // Fixed effects design matrix
  DATA_SPARSE_MATRIX(Z);       // Random effects design matrix (sparse)
  DATA_IVECTOR(groups);        // Group indices (0-indexed)
  DATA_INTEGER(n_groups);      // Number of groups
  DATA_INTEGER(n_obs);         // Number of observations
  DATA_INTEGER(n_fixed);       // Number of fixed effects
  DATA_INTEGER(n_random);      // Number of random effects per group
  DATA_INTEGER(correlated);    // 1 if correlated, 0 if uncorrelated

  // Parameters
  PARAMETER_VECTOR(beta);      // Fixed effects coefficients
  PARAMETER_VECTOR(u);         // Random effects (all groups concatenated)
  PARAMETER(log_sigma);        // Log residual standard deviation

  // Variance components
  PARAMETER_VECTOR(log_sigma_u); // Log random effects standard deviations
  PARAMETER_VECTOR(theta);       // Cholesky correlation parameters (if correlated)

  // Transform parameters
  Type sigma = exp(log_sigma);
  vector<Type> sigma_u = exp(log_sigma_u.array());

  // Build variance-covariance matrix for random effects
  matrix<Type> Sigma_u(n_random, n_random);

  if (correlated == 1 && n_random > 1) {
    // Construct correlation matrix from Cholesky factor
    int n_theta = n_random * (n_random - 1) / 2;
    vector<Type> theta_full(n_theta);
    for (int i = 0; i < n_theta; i++) {
      theta_full(i) = theta(i);
    }

    // Build Cholesky factor
    matrix<Type> L(n_random, n_random);
    L.setZero();

    int idx = 0;
    for (int i = 0; i < n_random; i++) {
      L(i, i) = Type(1.0);
      for (int j = 0; j < i; j++) {
        L(i, j) = theta_full(idx);
        idx++;
      }
    }

    // Correlation matrix: R = L * L'
    matrix<Type> R = L * L.transpose();

    // Scale to get covariance matrix
    for (int i = 0; i < n_random; i++) {
      for (int j = 0; j < n_random; j++) {
        Sigma_u(i, j) = sigma_u(i) * sigma_u(j) * R(i, j);
      }
    }
  } else {
    // Uncorrelated random effects (diagonal)
    Sigma_u.setZero();
    for (int i = 0; i < n_random; i++) {
      Sigma_u(i, i) = sigma_u(i) * sigma_u(i);
    }
  }

  // Compute inverse and determinant for multivariate normal
  matrix<Type> Sigma_u_inv = Sigma_u.inverse();
  Type log_det_Sigma_u = atomic::logdet(Sigma_u);

  // Initialize negative log-likelihood
  Type nll = 0.0;

  // Prior for random effects: u ~ MVN(0, Sigma_u)
  for (int j = 0; j < n_groups; j++) {
    vector<Type> u_j(n_random);
    for (int k = 0; k < n_random; k++) {
      int idx = j * n_random + k;
      u_j(k) = u(idx);
    }

    // Multivariate normal log-density
    Type quad_form = (u_j * (Sigma_u_inv * u_j)).sum();
    nll += 0.5 * (Type(n_random) * log(2.0 * M_PI) + log_det_Sigma_u + quad_form);
  }

  // Likelihood for observations
  for (int i = 0; i < n_obs; i++) {
    // Linear predictor: fixed effects + random effects
    Type eta = 0.0;

    // Add fixed effects
    for (int j = 0; j < n_fixed; j++) {
      eta += X(i, j) * beta(j);
    }

    // Add random effects
    int g = groups(i);  // Group for this observation
    for (int k = 0; k < n_random; k++) {
      int u_idx = g * n_random + k;
      eta += Z.coeff(i, k) * u(u_idx);
    }

    // Gaussian likelihood: y ~ N(eta, sigma^2)
    nll -= dnorm(y(i), eta, sigma, true);
  }

  // Report fitted values and variance components
  vector<Type> fitted(n_obs);
  for (int i = 0; i < n_obs; i++) {
    Type eta = 0.0;
    for (int j = 0; j < n_fixed; j++) {
      eta += X(i, j) * beta(j);
    }
    int g = groups(i);
    for (int k = 0; k < n_random; k++) {
      int u_idx = g * n_random + k;
      eta += Z.coeff(i, k) * u(u_idx);
    }
    fitted(i) = eta;
  }

  ADREPORT(fitted);
  ADREPORT(sigma);
  ADREPORT(sigma_u);
  ADREPORT(Sigma_u);

  return nll;
}
