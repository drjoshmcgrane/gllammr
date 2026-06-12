// Binomial GLMM with random effects
// Supports logit, probit, and complementary log-log links

#ifndef GLLAMM_BINOMIAL_HPP
#define GLLAMM_BINOMIAL_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_binomial(objective_function<Type>* obj)
{
  // Data inputs
  DATA_VECTOR(y);              // Response vector (0/1 or counts)
  DATA_MATRIX(X);              // Fixed effects design matrix
  DATA_SPARSE_MATRIX(Z);       // Random effects design matrix
  DATA_IVECTOR(groups);        // Group indices (0-indexed)
  DATA_INTEGER(n_groups);      // Number of groups
  DATA_INTEGER(n_obs);         // Number of observations
  DATA_INTEGER(n_fixed);       // Number of fixed effects
  DATA_INTEGER(n_random);      // Number of random effects per group
  DATA_INTEGER(link);          // Link function: 1=logit, 2=probit, 3=cloglog
  DATA_INTEGER(correlated);    // 1 if correlated, 0 if uncorrelated
  DATA_VECTOR(weights);
  DATA_VECTOR(group_weights);  // Level-2 weights (one per group; 1 = unweighted)        // Observation weights (pweights or fweights)

  // Parameters
  PARAMETER_VECTOR(beta);      // Fixed effects coefficients
  PARAMETER_VECTOR(u);         // Random effects
  PARAMETER_VECTOR(log_sigma_u); // Log random effects standard deviations
  PARAMETER_VECTOR(theta);     // Cholesky correlation parameters

  // Transform parameters
  vector<Type> sigma_u = exp(log_sigma_u.array());

  // Build variance-covariance matrix for random effects
  matrix<Type> Sigma_u(n_random, n_random);

  if (correlated == 1 && n_random > 1) {
    int n_theta = n_random * (n_random - 1) / 2;
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
    Sigma_u.setZero();
    for (int i = 0; i < n_random; i++) {
      Sigma_u(i, i) = sigma_u(i) * sigma_u(i);
    }
  }

  matrix<Type> Sigma_u_inv = Sigma_u.inverse();
  Type log_det_Sigma_u = atomic::logdet(Sigma_u);

  // Initialize negative log-likelihood
  // parallel_accumulator splits the likelihood across OpenMP threads
  // when available (no-op on single-threaded builds)
  parallel_accumulator<Type> nll(obj);

  // Prior for random effects, scaled by level-2 weights
  for (int j = 0; j < n_groups; j++) {
    vector<Type> u_j(n_random);
    for (int k = 0; k < n_random; k++) {
      u_j(k) = u(j * n_random + k);
    }
    Type quad_form = (u_j * (Sigma_u_inv * u_j)).sum();
    nll += group_weights(j) * Type(0.5) *
      (Type(n_random) * log(2.0 * M_PI) + log_det_Sigma_u + quad_form);
  }

  // Likelihood: single pass over observations (fitted values reuse the
  // same eta - a second loop would double the AD tape)
  vector<Type> eta_vec = X * beta;
  vector<Type> fitted(n_obs);
  for (int i = 0; i < n_obs; i++) {
    int g = groups(i);
    Type eta = eta_vec(i);
    for (int k = 0; k < n_random; k++) {
      eta += Z.coeff(i, k) * u(g * n_random + k);
    }

    Type p;
    if (link == 1) {
      p = invlogit(eta);
    } else if (link == 2) {
      p = pnorm(eta);
    } else {
      p = Type(1.0) - exp(-exp(eta));
    }

    Type w_i = weights(i) * group_weights(g);
    Type ll_i = y(i) * log(p + Type(1e-10)) +
                (Type(1.0) - y(i)) * log(Type(1.0) - p + Type(1e-10));
    nll -= w_i * ll_i;
    fitted(i) = p;
  }

  REPORT(fitted);
  ADREPORT(sigma_u);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_BINOMIAL_HPP
