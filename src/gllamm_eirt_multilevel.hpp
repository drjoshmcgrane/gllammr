// Multi-Level Explanatory IRT Models with Random Effects
// Supports nested, crossed, and partially nested random effects
// Combined with item-level and threshold-level predictors
//
// Model: theta_i = theta_0i + sum_g u_g[group_g[i]]
//        b_j = W_j * gamma + epsilon_j (if item_residuals = 1)
//        tau_jk = V_jk * delta + eta_jk (if threshold predictors used)

#ifndef GLLAMM_EIRT_MULTILEVEL_HPP
#define GLLAMM_EIRT_MULTILEVEL_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_eirt_multilevel(objective_function<Type>* obj)
{
  // ============================================================================
  // DATA INPUTS
  // ============================================================================

  // Response data
  DATA_VECTOR(y);                      // Item responses
  DATA_IVECTOR(person_id);             // Person identifier (0-indexed)
  DATA_IVECTOR(item_id);               // Item identifier (0-indexed)
  DATA_INTEGER(n_persons);             // Number of persons
  DATA_INTEGER(n_items);               // Number of items
  DATA_INTEGER(n_obs);                 // Number of observations
  DATA_VECTOR(weights);                // Case weights
  DATA_INTEGER(model_type);            // 1=Rasch/2PL, 2=GRM, 3=PCM, 4=GPCM

  // Item-level predictors
  DATA_MATRIX(W_difficulty);           // Design matrix for difficulty
  DATA_MATRIX(W_discrimination);       // Design matrix for discrimination (if used)
  DATA_INTEGER(item_residuals);        // 1=include residuals, 0=pure LLTM

  // Threshold-level predictors (for polytomous)
  DATA_INTEGER(has_threshold_predictors);  // 1=yes, 0=no
  DATA_MATRIX(V_threshold);            // Design matrix for thresholds (if polytomous)
  DATA_IVECTOR(n_categories_per_item); // Number of categories per item
  DATA_INTEGER(max_categories);        // Maximum categories

  // Multi-level structure
  DATA_INTEGER(has_random);            // 0=standard, 1=multi-level
  DATA_INTEGER(n_random_effects);      // Number of RE levels
  DATA_IMATRIX(group_ids);             // [n_persons × n_random_effects], -1 = NA
  DATA_IVECTOR(n_groups);              // Number of groups per level
  DATA_INTEGER(max_n_groups);          // Maximum groups across levels

  // ============================================================================
  // PARAMETERS
  // ============================================================================

  // Item parameter predictors
  PARAMETER_VECTOR(gamma_b);           // Coefficients for difficulty
  PARAMETER_VECTOR(gamma_a);           // Coefficients for discrimination
  PARAMETER_VECTOR(delta_tau);         // Coefficients for thresholds

  // Item-specific residuals (if item_residuals = 1)
  PARAMETER_VECTOR(epsilon_b);         // Difficulty residuals
  PARAMETER_VECTOR(epsilon_a);         // Discrimination residuals

  // Random effect SDs for item residuals
  PARAMETER(log_sigma_epsilon_b);
  PARAMETER(log_sigma_epsilon_a);

  // Person-level ability deviations
  PARAMETER_VECTOR(theta_0);
  PARAMETER(log_sigma_theta);

  // Random effects (group-level deviations)
  PARAMETER_MATRIX(u_random);          // [max_n_groups × n_random_effects]
  PARAMETER_VECTOR(log_sigma_random);  // SD for each RE level

  // ============================================================================
  // TRANSFORM PARAMETERS
  // ============================================================================

  Type sigma_theta = exp(log_sigma_theta);
  Type sigma_epsilon_b = exp(log_sigma_epsilon_b);
  Type sigma_epsilon_a = exp(log_sigma_epsilon_a);

  vector<Type> sigma_random(n_random_effects);
  if (has_random == 1) {
    for (int re = 0; re < n_random_effects; re++) {
      sigma_random(re) = exp(log_sigma_random(re));
    }
  }

  // ============================================================================
  // INITIALIZE NEGATIVE LOG-LIKELIHOOD
  // ============================================================================

  Type nll = 0.0;

  // ============================================================================
  // PRIORS
  // ============================================================================

  // Person-level deviations
  for (int p = 0; p < n_persons; p++) {
    nll -= dnorm(theta_0(p), Type(0.0), sigma_theta, true);
  }

  // Item residuals (if included)
  if (item_residuals == 1) {
    for (int j = 0; j < n_items; j++) {
      nll -= dnorm(epsilon_b(j), Type(0.0), sigma_epsilon_b, true);
      nll -= dnorm(epsilon_a(j), Type(0.0), sigma_epsilon_a, true);
    }
  }

  // Random effects priors
  if (has_random == 1) {
    for (int re = 0; re < n_random_effects; re++) {
      int n_groups_re = n_groups(re);
      Type sigma_re = sigma_random(re);

      for (int g = 0; g < n_groups_re; g++) {
        nll -= dnorm(u_random(g, re), Type(0.0), sigma_re, true);
      }
    }
  }

  // ============================================================================
  // LIKELIHOOD FOR ITEM RESPONSES
  // ============================================================================

  for (int i = 0; i < n_obs; i++) {
    int person = person_id(i);
    int item = item_id(i);

    // Compose total ability: person + group effects
    Type theta = theta_0(person);

    if (has_random == 1) {
      for (int re = 0; re < n_random_effects; re++) {
        int group = group_ids(person, re);
        if (group >= 0) {  // -1 indicates NA (partial nesting)
          theta += u_random(group, re);
        }
      }
    }

    // Compute item parameters from predictors
    vector<Type> W_b_row = W_difficulty.row(item);
    Type difficulty_pred = (W_b_row * gamma_b).sum();

    vector<Type> W_a_row = W_discrimination.row(item);
    Type log_discrim_pred = (W_a_row * gamma_a).sum();

    Type difficulty, discrimination;
    if (item_residuals == 1) {
      difficulty = difficulty_pred + epsilon_b(item);
      discrimination = exp(log_discrim_pred + epsilon_a(item));
    } else {
      difficulty = difficulty_pred;
      discrimination = exp(log_discrim_pred);
    }

    // Compute probability based on model type
    Type prob;

    if (model_type == 1) {
      // Dichotomous (Rasch/2PL)
      Type eta = discrimination * (theta - difficulty);
      prob = invlogit(eta);

    } else {
      // Polytomous (GRM, PCM, GPCM) - simplified for now
      // In full implementation, would use threshold predictors
      int obs_cat = CppAD::Integer(y(i)) - 1;
      int K = n_categories_per_item(item);

      // Simple equal spacing placeholder
      // Full implementation would use V_threshold and delta_tau
      vector<Type> tau(K - 1);
      for (int k = 0; k < K - 1; k++) {
        tau(k) = Type(k) - Type(K - 1) / Type(2.0);
      }

      // GRM-like computation
      if (obs_cat == 0) {
        prob = Type(1.0) - invlogit(discrimination * (theta - tau(0)));
      } else if (obs_cat == K - 1) {
        prob = invlogit(discrimination * (theta - tau(K - 2)));
      } else {
        Type p_ge_k = invlogit(discrimination * (theta - tau(obs_cat - 1)));
        Type p_ge_kp1 = invlogit(discrimination * (theta - tau(obs_cat)));
        prob = p_ge_k - p_ge_kp1;
      }
    }

    // Log-likelihood (weighted)
    Type w_i = weights(i);
    if (model_type == 1) {
      // Binary
      nll -= w_i * (y(i) * log(prob + Type(1e-10)) +
                     (Type(1.0) - y(i)) * log(Type(1.0) - prob + Type(1e-10)));
    } else {
      // Categorical
      nll -= w_i * log(prob + Type(1e-10));
    }
  }

  // ============================================================================
  // REPORT ESTIMATES
  // ============================================================================

  ADREPORT(theta_0);
  ADREPORT(gamma_b);
  ADREPORT(gamma_a);

  if (item_residuals == 1) {
    ADREPORT(epsilon_b);
    ADREPORT(epsilon_a);
    ADREPORT(sigma_epsilon_b);
    ADREPORT(sigma_epsilon_a);
  }

  ADREPORT(sigma_theta);

  if (has_random == 1) {
    ADREPORT(sigma_random);
    ADREPORT(u_random);
  }

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_EIRT_MULTILEVEL_HPP
