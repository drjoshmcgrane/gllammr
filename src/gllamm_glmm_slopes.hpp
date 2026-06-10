// GLMM with random intercepts and slopes (gaussian, binomial, poisson)
// Multiple random effects per group with optional correlation structure

#ifndef GLLAMM_GLMM_SLOPES_HPP
#define GLLAMM_GLMM_SLOPES_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_glmm_slopes(objective_function<Type>* obj)
{
  // Data inputs
  DATA_VECTOR(y);              // Response vector
  DATA_MATRIX(X);              // Fixed effects design matrix
  DATA_SPARSE_MATRIX(Z);       // Random effects design matrix (n_obs x n_random)
  DATA_IVECTOR(groups);        // Group indices (0-indexed)
  DATA_INTEGER(n_groups);      // Number of groups
  DATA_INTEGER(n_obs);         // Number of observations
  DATA_INTEGER(n_random);      // Number of random effects per group
  DATA_INTEGER(correlated);    // 1 if correlated, 0 if uncorrelated
  DATA_INTEGER(family);        // 0 = gaussian, 1 = binomial, 2 = poisson
  DATA_INTEGER(link);          // 1 = canonical (identity/logit/log),
                               // 2 = probit, 3 = cloglog (binomial only)
  DATA_VECTOR(weights);        // Case weights (fweights or pweights)

  // Parameters
  PARAMETER_VECTOR(beta);        // Fixed effects coefficients
  PARAMETER_VECTOR(u);           // Random effects (all groups concatenated)
  PARAMETER(log_sigma);          // Log residual SD (gaussian; mapped off otherwise)
  PARAMETER_VECTOR(log_sigma_u); // Log random effects standard deviations
  PARAMETER_VECTOR(theta);       // Cholesky correlation parameters (if correlated)

  Type sigma = exp(log_sigma);
  vector<Type> sigma_u = exp(log_sigma_u.array());

  // Build variance-covariance matrix for random effects
  matrix<Type> Sigma_u(n_random, n_random);

  if (correlated == 1 && n_random > 1) {
    // Unit-diagonal Cholesky factor of the correlation matrix
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
    matrix<Type> R = L * L.transpose();

    // Normalize to a proper correlation matrix, then scale by SDs
    for (int i = 0; i < n_random; i++) {
      for (int j = 0; j < n_random; j++) {
        Type rij = R(i, j) / sqrt(R(i, i) * R(j, j));
        Sigma_u(i, j) = sigma_u(i) * sigma_u(j) * rij;
      }
    }
  } else {
    Sigma_u.setZero();
    for (int i = 0; i < n_random; i++) {
      Sigma_u(i, i) = sigma_u(i) * sigma_u(i);
    }
  }

  matrix<Type> Sigma_u_inv = atomic::matinv(Sigma_u);
  Type log_det_Sigma_u = atomic::logdet(Sigma_u);

  // parallel_accumulator splits the likelihood across OpenMP threads
  // when available (no-op on single-threaded builds)
  parallel_accumulator<Type> nll(obj);

  // Prior for random effects: u_j ~ MVN(0, Sigma_u)
  for (int j = 0; j < n_groups; j++) {
    vector<Type> u_j = u.segment(j * n_random, n_random);
    Type quad_form = (u_j * (Sigma_u_inv * u_j)).sum();
    nll += 0.5 * (Type(n_random) * log(2.0 * M_PI) + log_det_Sigma_u + quad_form);
  }

  // Linear predictor: X beta computed once, RE contribution per observation
  vector<Type> eta = X * beta;
  for (int i = 0; i < n_obs; i++) {
    int g = groups(i);
    for (int k = 0; k < n_random; k++) {
      eta(i) += Z.coeff(i, k) * u(g * n_random + k);
    }
  }

  // Likelihood (single pass; fitted values reuse the same eta)
  vector<Type> fitted(n_obs);
  for (int i = 0; i < n_obs; i++) {
    Type w_i = weights(i);

    if (family == 0) {
      // Gaussian
      nll -= w_i * dnorm(y(i), eta(i), sigma, true);
      fitted(i) = eta(i);
    } else if (family == 1) {
      // Binomial
      Type p;
      if (link == 2) {
        p = pnorm(eta(i));
      } else if (link == 3) {
        p = Type(1.0) - exp(-exp(eta(i)));
      } else {
        p = invlogit(eta(i));
      }
      Type ll_i = y(i) * log(p + Type(1e-10)) +
                  (Type(1.0) - y(i)) * log(Type(1.0) - p + Type(1e-10));
      nll -= w_i * ll_i;
      fitted(i) = p;
    } else {
      // Poisson (log link)
      Type lambda = exp(eta(i));
      nll -= w_i * dpois(y(i), lambda, true);
      fitted(i) = lambda;
    }
  }

  REPORT(fitted);
  ADREPORT(sigma);
  ADREPORT(sigma_u);
  ADREPORT(Sigma_u);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_GLMM_SLOPES_HPP
