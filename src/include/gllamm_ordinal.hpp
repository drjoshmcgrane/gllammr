// Ordinal GLMM with proportional odds and cumulative probit
// Supports ordered categorical responses

#ifndef GLLAMM_ORDINAL_HPP
#define GLLAMM_ORDINAL_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_ordinal(objective_function<Type>* obj)
{
  // Data inputs
  DATA_IVECTOR(y);             // Ordinal response (1, 2, ..., K)
  DATA_MATRIX(X);              // Fixed effects design matrix
  DATA_SPARSE_MATRIX(Z);       // Random effects design matrix
  DATA_IVECTOR(groups);        // Group indices (0-indexed)
  DATA_INTEGER(n_groups);      // Number of groups
  DATA_INTEGER(n_obs);         // Number of observations
  DATA_INTEGER(n_fixed);       // Number of fixed effects
  DATA_INTEGER(n_random);      // Number of random effects per group
  DATA_INTEGER(n_categories);  // Number of ordinal categories
  DATA_INTEGER(link);          // Link: 1=logit, 2=probit, 3=acl, 4=crl_forward, 5=crl_backward, 6=ppo
  DATA_INTEGER(correlated);    // 1 if correlated, 0 if uncorrelated
  DATA_VECTOR(weights);        // Case weights (fweights or pweights)

  // Parameters
  PARAMETER_VECTOR(beta);      // Fixed effects coefficients (for links 1-5) or first threshold (for PPO)
  PARAMETER_VECTOR(u);         // Random effects
  PARAMETER_VECTOR(threshold); // Threshold parameters (length n_categories - 1)
  PARAMETER_VECTOR(log_sigma_u); // Log random effects standard deviations
  PARAMETER_VECTOR(theta);     // Cholesky correlation parameters
  PARAMETER_MATRIX(beta_ppo);  // PPO coefficients matrix (n_categories-1) x n_fixed (only for link=6)

  // Transform parameters
  vector<Type> sigma_u = exp(log_sigma_u.array());

  // Ensure threshold parameters are ordered (using cumulative sum of exp)
  vector<Type> ordered_threshold(n_categories - 1);
  ordered_threshold(0) = threshold(0);
  for (int k = 1; k < n_categories - 1; k++) {
    ordered_threshold(k) = ordered_threshold(k-1) + exp(threshold(k));
  }

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
    // Linear predictor
    Type eta = 0.0;

    // Fixed effects
    for (int j = 0; j < n_fixed; j++) {
      eta += X(i, j) * beta(j);
    }

    // Random effects
    int g = groups(i);
    for (int k = 0; k < n_random; k++) {
      eta += Z.coeff(i, k) * u(g * n_random + k);
    }

    // Observed category (1-indexed in R, convert to 0-indexed)
    int obs_cat = y(i) - 1;

    // Compute probability for this observation based on link function
    Type prob_cat;

    // Links 1-2: Cumulative (proportional odds/probit)
    if (link == 1 || link == 2) {
      if (obs_cat == 0) {
        // P(Y = 1) = P(Y <= 1)
        if (link == 1) {
          prob_cat = invlogit(ordered_threshold(0) - eta);
        } else {
          prob_cat = pnorm(ordered_threshold(0) - eta);
        }
      } else if (obs_cat == n_categories - 1) {
        // P(Y = K) = 1 - P(Y <= K-1)
        if (link == 1) {
          prob_cat = Type(1.0) - invlogit(ordered_threshold(n_categories - 2) - eta);
        } else {
          prob_cat = Type(1.0) - pnorm(ordered_threshold(n_categories - 2) - eta);
        }
      } else {
        // P(Y = k) = P(Y <= k) - P(Y <= k-1)
        Type prob_le_k, prob_le_k_minus_1;
        if (link == 1) {
          prob_le_k = invlogit(ordered_threshold(obs_cat) - eta);
          prob_le_k_minus_1 = invlogit(ordered_threshold(obs_cat - 1) - eta);
        } else {
          prob_le_k = pnorm(ordered_threshold(obs_cat) - eta);
          prob_le_k_minus_1 = pnorm(ordered_threshold(obs_cat - 1) - eta);
        }
        prob_cat = prob_le_k - prob_le_k_minus_1;
      }
    }

    // Link 3: Adjacent Category Logit (ACL)
    else if (link == 3) {
      // Compute log-odds for adjacent categories
      vector<Type> log_odds(n_categories - 1);
      for (int k = 0; k < n_categories - 1; k++) {
        log_odds(k) = ordered_threshold(k) + eta;
      }

      // Convert to probabilities via softmax
      // P(Y=0) = 1, P(Y=k) = exp(sum log_odds[0:k-1])
      vector<Type> log_prob(n_categories);
      log_prob(0) = Type(0.0);
      for (int k = 1; k < n_categories; k++) {
        log_prob(k) = log_prob(k-1) + log_odds(k-1);
      }

      // Normalize via log-sum-exp
      Type log_sum = log_prob(0);
      for (int k = 1; k < n_categories; k++) {
        log_sum = logspace_add(log_sum, log_prob(k));
      }

      prob_cat = exp(log_prob(obs_cat) - log_sum);
    }

    // Link 4: Continuation Ratio Logit (Forward)
    else if (link == 4) {
      // P(Y=k) = P(Y>=k) * P(Y=k | Y>=k)
      if (obs_cat == 0) {
        prob_cat = invlogit(ordered_threshold(0) - eta);
      } else if (obs_cat == n_categories - 1) {
        // Highest category: survived all previous thresholds
        Type prob_survive = Type(1.0);
        for (int j = 0; j < n_categories - 1; j++) {
          prob_survive *= (Type(1.0) - invlogit(ordered_threshold(j) - eta));
        }
        prob_cat = prob_survive;
      } else {
        Type prob_survive = Type(1.0);
        for (int j = 0; j < obs_cat; j++) {
          prob_survive *= (Type(1.0) - invlogit(ordered_threshold(j) - eta));
        }
        prob_cat = prob_survive * invlogit(ordered_threshold(obs_cat) - eta);
      }
    }

    // Link 5: Continuation Ratio Logit (Backward)
    else if (link == 5) {
      // P(Y=k) = P(Y<=k) * P(Y=k | Y<=k)
      if (obs_cat == n_categories - 1) {
        prob_cat = invlogit(ordered_threshold(n_categories - 2) - eta);
      } else {
        Type prob_at_or_below = invlogit(ordered_threshold(obs_cat) - eta);
        Type prob_below = Type(0.0);
        if (obs_cat > 0) {
          prob_below = invlogit(ordered_threshold(obs_cat - 1) - eta);
        }
        prob_cat = prob_at_or_below - prob_below;
      }
    }

    // Link 6: Partial Proportional Odds (PPO)
    else if (link == 6) {
      // Each threshold has its own covariate effects
      if (obs_cat == 0) {
        // Compute eta for first threshold
        Type eta_k = Type(0.0);
        for (int j = 0; j < n_fixed; j++) {
          eta_k += X(i, j) * beta_ppo(0, j);
        }
        // Add random effects
        int g = groups(i);
        for (int k = 0; k < n_random; k++) {
          eta_k += Z.coeff(i, k) * u(g * n_random + k);
        }
        prob_cat = invlogit(ordered_threshold(0) - eta_k);

      } else if (obs_cat == n_categories - 1) {
        // Compute eta for last threshold
        Type eta_k = Type(0.0);
        for (int j = 0; j < n_fixed; j++) {
          eta_k += X(i, j) * beta_ppo(n_categories - 2, j);
        }
        int g = groups(i);
        for (int k = 0; k < n_random; k++) {
          eta_k += Z.coeff(i, k) * u(g * n_random + k);
        }
        prob_cat = Type(1.0) - invlogit(ordered_threshold(n_categories - 2) - eta_k);

      } else {
        // Compute etas for both thresholds
        Type eta_k = Type(0.0);
        Type eta_k_minus_1 = Type(0.0);
        for (int j = 0; j < n_fixed; j++) {
          eta_k += X(i, j) * beta_ppo(obs_cat, j);
          eta_k_minus_1 += X(i, j) * beta_ppo(obs_cat - 1, j);
        }
        int g = groups(i);
        for (int k = 0; k < n_random; k++) {
          eta_k += Z.coeff(i, k) * u(g * n_random + k);
          eta_k_minus_1 += Z.coeff(i, k) * u(g * n_random + k);
        }
        Type prob_le_k = invlogit(ordered_threshold(obs_cat) - eta_k);
        Type prob_le_k_minus_1 = invlogit(ordered_threshold(obs_cat - 1) - eta_k_minus_1);
        prob_cat = prob_le_k - prob_le_k_minus_1;
      }
    }

    // Default: should not reach here
    else {
      prob_cat = Type(1e-10);  // Avoid uninitialized variable
      Rf_error("Invalid link function code");
    }

    // Add to negative log-likelihood (weighted)
    Type w_i = weights(i);
    nll -= w_i * log(prob_cat + Type(1e-10));
  }

  // Report
  ADREPORT(ordered_threshold);
  ADREPORT(sigma_u);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_ORDINAL_HPP
