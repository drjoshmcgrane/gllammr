// Mixed Response GLMM
// Multiple outcomes of different types sharing random effects

#ifndef GLLAMM_MIXED_RESPONSE_HPP
#define GLLAMM_MIXED_RESPONSE_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_mixed_response(objective_function<Type>* obj)
{
  // Data inputs
  DATA_VECTOR(y1);             // Continuous outcome
  DATA_IVECTOR(y2);            // Binary outcome
  DATA_VECTOR(y3);             // Count outcome
  DATA_MATRIX(X1);             // Design matrix for outcome 1
  DATA_MATRIX(X2);             // Design matrix for outcome 2
  DATA_MATRIX(X3);             // Design matrix for outcome 3
  DATA_SPARSE_MATRIX(Z);       // Shared random effects design
  DATA_IVECTOR(groups);        // Group indices
  DATA_INTEGER(n_groups);      // Number of groups
  DATA_INTEGER(n1);            // Number of obs for outcome 1
  DATA_INTEGER(n2);            // Number of obs for outcome 2
  DATA_INTEGER(n3);            // Number of obs for outcome 3
  DATA_INTEGER(n_fixed1);      // Number of fixed effects for outcome 1
  DATA_INTEGER(n_fixed2);      // Number of fixed effects for outcome 2
  DATA_INTEGER(n_fixed3);      // Number of fixed effects for outcome 3
  DATA_INTEGER(n_random);      // Number of random effects
  DATA_INTEGER(has_y1);        // Indicator: is y1 present?
  DATA_INTEGER(has_y2);        // Indicator: is y2 present?
  DATA_INTEGER(has_y3);        // Indicator: is y3 present?

  // Parameters
  PARAMETER_VECTOR(beta1);     // Fixed effects for outcome 1
  PARAMETER_VECTOR(beta2);     // Fixed effects for outcome 2
  PARAMETER_VECTOR(beta3);     // Fixed effects for outcome 3
  PARAMETER_VECTOR(u);         // Shared random effects
  PARAMETER(log_sigma1);       // Log residual SD for outcome 1
  PARAMETER_VECTOR(log_sigma_u); // Log random effects SDs
  PARAMETER_VECTOR(theta);     // Cholesky correlation parameters

  // Transform parameters
  Type sigma1 = exp(log_sigma1);
  vector<Type> sigma_u = exp(log_sigma_u.array());

  // Build variance-covariance matrix for random effects
  matrix<Type> Sigma_u(n_random, n_random);

  if (n_random > 1) {
    matrix<Type> L(n_random, n_random);
    L.setZero();

    int idx = 0;
    for (int i = 0; i < n_random; i++) {
      L(i, i) = Type(1.0);
      for (int j = 0; j < i; j++) {
        L(i, j) = theta(idx);
        idx++;
      }
    }

    // L L' rescaled to unit diagonal so sigma_u are genuine standard
    // deviations (otherwise the scale is unidentified)
    matrix<Type> R = L * L.transpose();

    for (int i = 0; i < n_random; i++) {
      for (int j = 0; j < n_random; j++) {
        Type rij = R(i, j) / sqrt(R(i, i) * R(j, j));
        Sigma_u(i, j) = sigma_u(i) * sigma_u(j) * rij;
      }
    }
  } else {
    Sigma_u(0, 0) = sigma_u(0) * sigma_u(0);
  }

  matrix<Type> Sigma_u_inv = Sigma_u.inverse();
  Type log_det_Sigma_u = atomic::logdet(Sigma_u);

  // Initialize negative log-likelihood
  Type nll = 0.0;

  // Prior for random effects
  for (int j = 0; j < n_groups; j++) {
    vector<Type> u_j(n_random);
    for (int k = 0; k < n_random; k++) {
      u_j(k) = u(j * n_random + k);
    }
    Type quad_form = (u_j * (Sigma_u_inv * u_j)).sum();
    nll += 0.5 * (Type(n_random) * log(2.0 * M_PI) + log_det_Sigma_u + quad_form);
  }

  // Likelihood for outcome 1 (Gaussian)
  if (has_y1 == 1) {
    for (int i = 0; i < n1; i++) {
      Type eta = 0.0;
      for (int j = 0; j < n_fixed1; j++) {
        eta += X1(i, j) * beta1(j);
      }
      int g = groups(i);
      for (int k = 0; k < n_random; k++) {
        eta += Z.coeff(i, k) * u(g * n_random + k);
      }
      nll -= dnorm(y1(i), eta, sigma1, true);
    }
  }

  // Likelihood for outcome 2 (Binomial/logit)
  if (has_y2 == 1) {
    for (int i = 0; i < n2; i++) {
      Type eta = 0.0;
      for (int j = 0; j < n_fixed2; j++) {
        eta += X2(i, j) * beta2(j);
      }
      int g = groups(i);
      for (int k = 0; k < n_random; k++) {
        eta += Z.coeff(i, k) * u(g * n_random + k);
      }
      Type p = invlogit(eta);
      nll -= y2(i) * log(p + Type(1e-10)) + (Type(1.0) - y2(i)) * log(Type(1.0) - p + Type(1e-10));
    }
  }

  // Likelihood for outcome 3 (Poisson/log)
  if (has_y3 == 1) {
    for (int i = 0; i < n3; i++) {
      Type eta = 0.0;
      for (int j = 0; j < n_fixed3; j++) {
        eta += X3(i, j) * beta3(j);
      }
      int g = groups(i);
      for (int k = 0; k < n_random; k++) {
        eta += Z.coeff(i, k) * u(g * n_random + k);
      }
      Type lambda = exp(eta);
      nll -= dpois(y3(i), lambda, true);
    }
  }

  // Report
  ADREPORT(sigma1);
  ADREPORT(sigma_u);
  ADREPORT(Sigma_u);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_MIXED_RESPONSE_HPP
