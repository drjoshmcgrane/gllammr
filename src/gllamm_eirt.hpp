// Explanatory IRT: Item parameters as functions of item covariates
// Supports dichotomous (Rasch, 2PL) and polytomous (GRM, PCM, GPCM)
//
// Polytomous model framework follows Kim & Wilson (2019) / De Boeck & Wilson EIRM:
//   GRM  (poly_model_type=1): Cumulative logit, ordered thresholds
//   PCM  (poly_model_type=2): Adjacent-categories logit, two-fold parameterization
//                              delta_im = b_i + s_im (MFRM approach)
//   GPCM (poly_model_type=3): As PCM with item-specific discrimination
//   PCM with threshold_formula (poly_model_type=4): Step-difficulty regression (LPCM+UISE)
//                              delta_im = b_i + sum_k xi_km * x_ik + e_im

#ifndef GLLAMM_EIRT_HPP
#define GLLAMM_EIRT_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_eirt(objective_function<Type>* obj)
{
  // ============================================================================
  // DATA INPUTS
  // ============================================================================

  DATA_VECTOR(y);                    // Item responses
  DATA_IVECTOR(person_id);           // Person identifier (0-indexed)
  DATA_IVECTOR(item_id);             // Item identifier (0-indexed)
  DATA_MATRIX(W_difficulty);         // Item covariate matrix for difficulty [n_items x p_diff]
  DATA_MATRIX(W_discrimination);     // Item covariate matrix for discrimination [n_items x p_disc]
  DATA_MATRIX(W_threshold);          // Item covariate matrix for threshold regression [n_items x p_thresh]
                                     // Used for LPCM; dummy (zero-column) for other models
  DATA_INTEGER(n_persons);           // Number of persons
  DATA_INTEGER(n_items);             // Number of items
  DATA_INTEGER(n_obs);               // Number of observations
  DATA_VECTOR(weights);              // Case weights (fweights or pweights)
  DATA_INTEGER(model_type);          // 1=Rasch/1PL, 2=2PL (dichotomous only)
  DATA_INTEGER(is_polytomous);       // 0=dichotomous, 1=polytomous
  DATA_INTEGER(poly_model_type);     // 1=GRM, 2=PCM, 3=GPCM, 4=LPCM (polytomous only)
  DATA_INTEGER(item_residuals);      // 1=include item residuals (LLTM+error), 0=pure LLTM
  DATA_INTEGER(uses_discrimination); // 1=model reads discrimination (2PL/GRM/GPCM), 0=otherwise

  // For polytomous models
  DATA_IVECTOR(n_categories_per_item);  // Number of categories per item
  DATA_INTEGER(max_categories);          // Maximum K

  // Multi-level structure (group-level ability random effects)
  DATA_INTEGER(has_random);            // 0 = standard, 1 = multi-level
  DATA_INTEGER(n_random_effects);      // Number of RE levels
  DATA_IMATRIX(group_ids);             // [n_persons x n_random_effects], -1 = NA
  DATA_IVECTOR(n_groups);              // Number of groups per level

  // ============================================================================
  // PARAMETERS
  // ============================================================================

  PARAMETER_VECTOR(theta);           // Person abilities
  PARAMETER_VECTOR(gamma);           // Difficulty regression coefficients
  PARAMETER_VECTOR(delta);           // Discrimination regression coefficients
  PARAMETER_VECTOR(epsilon_b);       // Item-specific difficulty residuals
  PARAMETER_VECTOR(epsilon_a);       // Item-specific discrimination residuals
  PARAMETER(log_sigma_epsilon_b);    // Log SD of difficulty residuals
  PARAMETER(log_sigma_epsilon_a);    // Log SD of discrimination residuals
  PARAMETER(log_sigma_theta);        // Log SD of ability distribution

  // For GRM polytomous: threshold spacing parameters [n_items x (max_categories-1)]
  // Column 0: first threshold offset from item location
  // Columns 1+: log-spacing for subsequent thresholds
  // For PCM/GPCM: step deviation parameters [n_items x (max_categories-2)]
  // (K-2 free deviations; last deviation = -sum of others, enforcing sum-to-zero)
  PARAMETER_MATRIX(step_param);

  // For LPCM: step-specific regression weights [p_thresh x (max_categories-1)]
  PARAMETER_MATRIX(xi);

  // For LPCM: item-step residuals [n_items x (max_categories-1)]
  PARAMETER_MATRIX(e_step);
  PARAMETER(log_sigma_e_step);       // Log SD of step residuals

  // Multi-level: group-level ability deviations
  PARAMETER_MATRIX(u_random);          // [max_n_groups x n_random_effects]
  PARAMETER_VECTOR(log_sigma_random);  // SD per RE level

  // ============================================================================
  // TRANSFORM PARAMETERS
  // ============================================================================

  Type sigma_epsilon_b = exp(log_sigma_epsilon_b);
  Type sigma_epsilon_a = exp(log_sigma_epsilon_a);
  Type sigma_theta = exp(log_sigma_theta);
  Type sigma_e_step = exp(log_sigma_e_step);

  // ============================================================================
  // INITIALIZE NEGATIVE LOG-LIKELIHOOD
  // ============================================================================

  // parallel_accumulator splits the likelihood across OpenMP threads
  // when available (no-op on single-threaded builds)
  parallel_accumulator<Type> nll(obj);

  // ============================================================================
  // PRIORS
  // ============================================================================

  // Prior for person abilities: theta ~ N(0, sigma_theta^2)
  for (int p = 0; p < n_persons; p++) {
    nll -= dnorm(theta(p), Type(0.0), sigma_theta, true);
  }

  // Prior for item residuals (only if item_residuals == 1). The epsilon_a
  // prior is skipped when the model never reads discrimination: evaluating
  // it on a fixed value would shift the log-likelihood by a constant.
  if (item_residuals == 1) {
    for (int j = 0; j < n_items; j++) {
      nll -= dnorm(epsilon_b(j), Type(0.0), sigma_epsilon_b, true);
      if (uses_discrimination == 1) {
        nll -= dnorm(epsilon_a(j), Type(0.0), sigma_epsilon_a, true);
      }
    }
  }

  // Prior for LPCM step residuals
  if (is_polytomous == 1 && poly_model_type == 4) {
    for (int j = 0; j < n_items; j++) {
      int K = n_categories_per_item(j);
      for (int m = 0; m < K - 1; m++) {
        nll -= dnorm(e_step(j, m), Type(0.0), sigma_e_step, true);
      }
    }
  }

  // Prior for group-level random effects
  if (has_random == 1) {
    for (int re = 0; re < n_random_effects; re++) {
      Type sigma_re = exp(log_sigma_random(re));
      for (int g = 0; g < n_groups(re); g++) {
        nll -= dnorm(u_random(g, re), Type(0.0), sigma_re, true);
      }
    }
  }

  // Effective ability: person deviation plus group-level effects
  vector<Type> theta_eff(n_persons);
  for (int p = 0; p < n_persons; p++) {
    theta_eff(p) = theta(p);
    if (has_random == 1) {
      for (int re = 0; re < n_random_effects; re++) {
        int g = group_ids(p, re);
        if (g >= 0) {  // -1 indicates NA (partial nesting)
          theta_eff(p) += u_random(g, re);
        }
      }
    }
  }

  // ============================================================================
  // COMPUTE ITEM PARAMETERS FROM COVARIATES
  // ============================================================================

  vector<Type> difficulty(n_items);
  vector<Type> discrimination(n_items);

  int p_diff = W_difficulty.cols();
  int p_disc = W_discrimination.cols();

  for (int j = 0; j < n_items; j++) {
    // Item location (overall difficulty) as linear function of item covariates
    Type difficulty_pred = 0.0;
    for (int p = 0; p < p_diff; p++) {
      difficulty_pred += gamma(p) * W_difficulty(j, p);
    }

    // Log-discrimination as linear function of item covariates
    Type log_discrim_pred = 0.0;
    for (int p = 0; p < p_disc; p++) {
      log_discrim_pred += delta(p) * W_discrimination(j, p);
    }

    // Add item residuals if requested (LLTM+error vs pure LLTM)
    if (item_residuals == 1) {
      difficulty(j) = difficulty_pred + epsilon_b(j);
      discrimination(j) = exp(log_discrim_pred + epsilon_a(j));
    } else {
      difficulty(j) = difficulty_pred;
      discrimination(j) = exp(log_discrim_pred);
    }
  }

  // ============================================================================
  // LIKELIHOOD FOR ITEM RESPONSES
  // ============================================================================

  if (is_polytomous == 0) {
    // ========================================================================
    // DICHOTOMOUS IRT MODELS
    // ========================================================================

    for (int i = 0; i < n_obs; i++) {
      int person = person_id(i);
      int item = item_id(i);

      Type prob;

      if (model_type == 1) {
        // Rasch model: P(Y=1) = logit^{-1}(theta - b)
        Type eta = theta_eff(person) - difficulty(item);
        prob = invlogit(eta);

      } else {
        // 2PL model: P(Y=1) = logit^{-1}(a*(theta - b))
        Type eta = discrimination(item) * (theta_eff(person) - difficulty(item));
        prob = invlogit(eta);
      }

      // Bernoulli log-likelihood (weighted)
      Type w_i = weights(i);
      nll -= w_i * (y(i) * log(prob + Type(1e-10)) + (Type(1.0) - y(i)) * log(Type(1.0) - prob + Type(1e-10)));
    }

  } else {
    // ========================================================================
    // POLYTOMOUS IRT MODELS
    // ========================================================================

    int p_thresh = W_threshold.cols();

    for (int i = 0; i < n_obs; i++) {
      int person = person_id(i);
      int item = item_id(i);
      int obs_cat = CppAD::Integer(y(i)) - 1;  // Convert to 0-indexed
      int K = n_categories_per_item(item);

      Type prob_cat = Type(0.0);

      // ====================================================================
      // GRM: Cumulative logit with ordered thresholds
      // b_i provides item location, step_param provides spacing
      // ====================================================================
      if (poly_model_type == 1) {
        // Build ordered thresholds:
        //   tau_1 = b_i + step_param(item, 0)
        //   tau_k = tau_{k-1} + exp(step_param(item, k-1))  for k >= 2
        vector<Type> ordered_threshold(K - 1);
        ordered_threshold(0) = difficulty(item) + step_param(item, 0);
        for (int k = 1; k < K - 1; k++) {
          ordered_threshold(k) = ordered_threshold(k-1) + exp(step_param(item, k));
        }

        Type a = discrimination(item);
        // Standard IRT GRM: P(Y >= k) = invlogit(a*(theta - tau_{k-1}))
        // P(Y = k) = P(Y >= k) - P(Y >= k+1), requires ordered thresholds (tau increasing)
        if (obs_cat == 0) {
          // P(Y = 0) = 1 - P(Y >= 1) = 1 - invlogit(a*(theta - tau_0))
          prob_cat = Type(1.0) - invlogit(a * (theta_eff(person) - ordered_threshold(0)));
        } else if (obs_cat == K - 1) {
          // P(Y = K-1) = P(Y >= K-1) = invlogit(a*(theta - tau_{K-2}))
          prob_cat = invlogit(a * (theta_eff(person) - ordered_threshold(K - 2)));
        } else {
          // P(Y = k) = invlogit(a*(theta - tau_{k-1})) - invlogit(a*(theta - tau_k))
          Type p_ge_k = invlogit(a * (theta_eff(person) - ordered_threshold(obs_cat - 1)));
          Type p_ge_k_plus_1 = invlogit(a * (theta_eff(person) - ordered_threshold(obs_cat)));
          prob_cat = p_ge_k - p_ge_k_plus_1;
        }
      }

      // ====================================================================
      // PCM: Adjacent-categories logit, two-fold decomposition (MFRM)
      // delta_im = b_i + s_im
      //   b_i = difficulty(item)  [item location from difficulty_formula]
      //   s_im = step deviations with sum-to-zero constraint
      //   s_{i,K-2} = -(s_{i,0} + ... + s_{i,K-3})
      // ====================================================================
      else if (poly_model_type == 2) {
        vector<Type> cumsum(K);
        cumsum(0) = Type(0.0);

        for (int m = 1; m < K; m++) {
          // Compute step deviation with sum-to-zero constraint
          Type s_im;
          if (m <= K - 2) {
            s_im = step_param(item, m - 1);
          } else {
            // Last step deviation = negative sum of all others
            Type sum_s = Type(0.0);
            for (int q = 0; q < K - 2; q++) {
              sum_s += step_param(item, q);
            }
            s_im = -sum_s;
          }
          Type delta_im = difficulty(item) + s_im;
          cumsum(m) = cumsum(m-1) + (theta_eff(person) - delta_im);
        }

        Type denom = Type(0.0);
        for (int m = 0; m < K; m++) {
          denom += exp(cumsum(m));
        }
        prob_cat = exp(cumsum(obs_cat)) / denom;
      }

      // ====================================================================
      // GPCM: Adjacent-categories logit with discrimination, two-fold
      // Same as PCM but with item-specific discrimination a_i
      // ====================================================================
      else if (poly_model_type == 3) {
        Type a = discrimination(item);

        vector<Type> cumsum(K);
        cumsum(0) = Type(0.0);

        for (int m = 1; m < K; m++) {
          Type s_im;
          if (m <= K - 2) {
            s_im = step_param(item, m - 1);
          } else {
            Type sum_s = Type(0.0);
            for (int q = 0; q < K - 2; q++) {
              sum_s += step_param(item, q);
            }
            s_im = -sum_s;
          }
          Type delta_im = difficulty(item) + s_im;
          cumsum(m) = cumsum(m-1) + a * (theta_eff(person) - delta_im);
        }

        Type denom = Type(0.0);
        for (int m = 0; m < K; m++) {
          denom += exp(cumsum(m));
        }
        prob_cat = exp(cumsum(obs_cat)) / denom;
      }

      // ====================================================================
      // LPCM: Threshold-difficulty regression (Kim & Wilson LPCM+UISE)
      // delta_im = b_i + sum_k xi_{k,m} * x_{i,k} + e_{i,m}
      // b_i = item location (from difficulty_formula)
      // W_threshold [n_items x p_thresh], xi [p_thresh x (max_categories-1)]
      // ====================================================================
      else if (poly_model_type == 4) {
        vector<Type> cumsum(K);
        cumsum(0) = Type(0.0);

        for (int m = 1; m < K; m++) {
          // Threshold difficulty = item location + threshold-specific regression
          Type delta_im = difficulty(item);  // b_i (from difficulty_formula)
          for (int p = 0; p < p_thresh; p++) {
            delta_im += W_threshold(item, p) * xi(p, m - 1);  // threshold covariates
          }
          delta_im += e_step(item, m - 1);  // threshold residual

          cumsum(m) = cumsum(m-1) + (theta_eff(person) - delta_im);
        }

        Type denom = Type(0.0);
        for (int m = 0; m < K; m++) {
          denom += exp(cumsum(m));
        }
        prob_cat = exp(cumsum(obs_cat)) / denom;
      }

      // Categorical log-likelihood (weighted)
      Type w_i = weights(i);
      nll -= w_i * log(prob_cat + Type(1e-10));
    }
  }

  // ============================================================================
  // REPORT ESTIMATES
  // ============================================================================

  REPORT(theta);
  ADREPORT(gamma);
  ADREPORT(delta);
  ADREPORT(difficulty);
  ADREPORT(discrimination);
  ADREPORT(sigma_epsilon_b);
  ADREPORT(sigma_epsilon_a);
  ADREPORT(sigma_theta);

  if (is_polytomous == 1) {
    if (poly_model_type == 4) {
      ADREPORT(xi);
      ADREPORT(sigma_e_step);
    }
  }

  if (has_random == 1) {
    vector<Type> sigma_random = exp(log_sigma_random.array());
    ADREPORT(sigma_random);
  }

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_EIRT_HPP
