// Multinomial GLMM with baseline category logit
// Supports unordered categorical responses

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data inputs
  DATA_IVECTOR(y);             // Nominal response (0, 1, ..., K-1)
  DATA_MATRIX(X);              // Fixed effects design matrix
  DATA_SPARSE_MATRIX(Z);       // Random effects design matrix
  DATA_IVECTOR(groups);        // Group indices (0-indexed)
  DATA_INTEGER(n_groups);      // Number of groups
  DATA_INTEGER(n_obs);         // Number of observations
  DATA_INTEGER(n_fixed);       // Number of fixed effects PER category
  DATA_INTEGER(n_random);      // Number of random effects per group
  DATA_INTEGER(n_categories);  // Number of nominal categories
  DATA_INTEGER(correlated);    // 1 if correlated, 0 if uncorrelated

  // Parameters
  // beta is a matrix: (n_categories - 1) x n_fixed
  // Each row contains coefficients for one non-reference category
  PARAMETER_MATRIX(beta);      // Fixed effects coefficients
  PARAMETER_VECTOR(u);         // Random effects
  PARAMETER_VECTOR(log_sigma_u); // Log random effects standard deviations
  PARAMETER_VECTOR(theta);     // Cholesky correlation parameters

  // Transform parameters
  vector<Type> sigma_u = exp(log_sigma_u.array());

  // Build variance-covariance matrix for random effects
  matrix<Type> Sigma_u(n_random, n_random);

  if (correlated == 1 && n_random > 1) {
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

    for (int i = 0; i < n_random; i++) {
      for (int j = 0; j < n_random; j++) {
        Sigma_u(i, j) = sigma_u(i) * sigma_u(j) * R(i, j);
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

  // Likelihood for observations
  for (int i = 0; i < n_obs; i++) {
    int g = groups(i);
    int obs_cat = y(i);

    // Compute linear predictors for all categories
    vector<Type> eta(n_categories);
    eta.setZero();

    // Reference category (category 0) has eta = 0
    // For other categories, compute eta
    for (int cat = 1; cat < n_categories; cat++) {
      // Fixed effects for this category
      for (int j = 0; j < n_fixed; j++) {
        eta(cat) += X(i, j) * beta(cat - 1, j);
      }

      // Random effects for this category
      for (int k = 0; k < n_random; k++) {
        eta(cat) += Z.coeff(i, k) * u(g * n_random + k);
      }
    }

    // Compute probabilities using softmax
    Type sum_exp = Type(1.0); // For reference category
    for (int cat = 1; cat < n_categories; cat++) {
      sum_exp += exp(eta(cat));
    }

    Type prob_obs_cat;
    if (obs_cat == 0) {
      prob_obs_cat = Type(1.0) / sum_exp;
    } else {
      prob_obs_cat = exp(eta(obs_cat)) / sum_exp;
    }

    // Add to negative log-likelihood
    nll -= log(prob_obs_cat + Type(1e-10));
  }

  // Report
  ADREPORT(sigma_u);

  return nll;
}
