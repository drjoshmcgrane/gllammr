// Multi-Level Polytomous IRT Models: GRM, PCM, GPCM, NRM with Random Effects
// Supports nested, crossed, and partially nested random effects
//
// Model: theta_i = theta_0i + sum_g u_g[group_g[i]]
//        theta_0i ~ N(0, sigma_theta²)
//        u_g ~ N(0, sigma_g²)

#ifndef GLLAMM_IRT_POLY_MULTILEVEL_HPP
#define GLLAMM_IRT_POLY_MULTILEVEL_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_irt_poly_multilevel(objective_function<Type>* obj)
{
  // ============================================================================
  // DATA INPUTS
  // ============================================================================

  // Standard polytomous IRT data
  DATA_VECTOR(y);                         // Item responses (1, 2, ..., K)
  DATA_IVECTOR(person_id);                // Person identifier (0-indexed)
  DATA_IVECTOR(item_id);                  // Item identifier (0-indexed)
  DATA_IVECTOR(n_categories_per_item);    // Number of categories for each item
  DATA_INTEGER(max_categories);           // Maximum number of categories
  DATA_INTEGER(n_persons);                // Number of persons
  DATA_INTEGER(n_items);                  // Number of items
  DATA_INTEGER(n_obs);                    // Number of observations
  DATA_VECTOR(weights);                   // Case weights
  DATA_INTEGER(model_type);               // 1=GRM, 2=PCM, 3=GPCM, 4=NRM

  // Multi-level structure
  DATA_INTEGER(has_random);               // 0=standard IRT, 1=multi-level
  DATA_INTEGER(n_random_effects);         // Number of RE levels
  DATA_IMATRIX(group_ids);                // [n_persons × n_random_effects], -1 = NA
  DATA_IVECTOR(n_groups);                 // Number of groups per level
  DATA_INTEGER(max_n_groups);             // Maximum groups across levels

  // ============================================================================
  // PARAMETERS
  // ============================================================================

  // Item parameters
  PARAMETER_MATRIX(threshold_raw);        // Item threshold parameters [n_items × (max_categories-1)]
  PARAMETER_VECTOR(discrimination);       // Item discriminations

  // Person-level ability deviations
  PARAMETER_VECTOR(theta_0);
  PARAMETER(log_sigma_theta);

  // Random effects (group-level deviations)
  PARAMETER_MATRIX(u_random);             // [max_n_groups × n_random_effects]
  PARAMETER_VECTOR(log_sigma_random);     // SD for each RE level

  // ============================================================================
  // TRANSFORM PARAMETERS
  // ============================================================================

  Type sigma_theta = exp(log_sigma_theta);

  vector<Type> sigma_random(n_random_effects);
  if (has_random == 1) {
    for (int re = 0; re < n_random_effects; re++) {
      sigma_random(re) = exp(log_sigma_random(re));
    }
  }

  // ============================================================================
  // INITIALIZE NEGATIVE LOG-LIKELIHOOD
  // ============================================================================

  // parallel_accumulator splits the likelihood across OpenMP threads
  // when available (no-op on single-threaded builds)
  parallel_accumulator<Type> nll(obj);

  // ============================================================================
  // PRIORS
  // ============================================================================

  // Person-level deviations: theta_0 ~ N(0, sigma_theta²)
  for (int p = 0; p < n_persons; p++) {
    nll -= dnorm(theta_0(p), Type(0.0), sigma_theta, true);
  }

  // Random effects priors: u_g ~ N(0, sigma_g²)
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
    int obs_cat = CppAD::Integer(y(i)) - 1;  // Convert to 0-indexed

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

    // Get item parameters
    int K = n_categories_per_item(item);
    Type a = discrimination(item);

    // Build ordered thresholds for this item
    vector<Type> tau(K - 1);
    if (model_type == 1) {
      // GRM: ordered parameterization (threshold_raw[item,0] + sum of exp(threshold_raw))
      tau(0) = threshold_raw(item, 0);
      for (int k = 1; k < K - 1; k++) {
        tau(k) = tau(k-1) + exp(threshold_raw(item, k));
      }
    } else {
      // PCM/GPCM: free step difficulties
      for (int k = 0; k < K - 1; k++) {
        tau(k) = threshold_raw(item, k);
      }
    }

    // Compute category probabilities
    Type prob_cat;

    if (model_type == 1) {
      // GRM: Graded Response Model
      // P(Y = k) = P(Y >= k) - P(Y >= k+1)
      // P(Y >= k) = logit^(-1)(a * (theta - tau[k-1]))

      if (obs_cat == 0) {
        // P(Y = 0) = 1 - P(Y >= 1)
        prob_cat = Type(1.0) - invlogit(a * (theta - tau(0)));
      } else if (obs_cat == K - 1) {
        // P(Y = K-1) = P(Y >= K-1)
        prob_cat = invlogit(a * (theta - tau(K - 2)));
      } else {
        // P(Y = k) = P(Y >= k) - P(Y >= k+1)
        Type p_ge_k = invlogit(a * (theta - tau(obs_cat - 1)));
        Type p_ge_kp1 = invlogit(a * (theta - tau(obs_cat)));
        prob_cat = p_ge_k - p_ge_kp1;
      }

    } else if (model_type == 2 || model_type == 3) {
      // PCM/GPCM: Partial Credit Model
      // For PCM: a = 1 (constrained in R code)
      // For GPCM: a is free

      // Compute numerators for each category
      vector<Type> numerators(K);
      numerators(0) = Type(0.0);  // log-numerator for category 0

      for (int k = 1; k < K; k++) {
        Type sum_tau = Type(0.0);
        for (int m = 0; m < k; m++) {
          sum_tau += tau(m);
        }
        numerators(k) = a * (Type(k) * theta - sum_tau);
      }

      // Log-sum-exp for denominator
      Type log_denom = numerators(0);
      for (int k = 1; k < K; k++) {
        log_denom = logspace_add(log_denom, numerators(k));
      }

      // Category probability
      prob_cat = exp(numerators(obs_cat) - log_denom);

    } else {
      // NRM: Nominal Response Model (future implementation)
      prob_cat = Type(1.0) / Type(K);  // Placeholder
    }

    // Bernoulli/categorical log-likelihood (weighted)
    Type w_i = weights(i);
    nll -= w_i * log(prob_cat + Type(1e-10));
  }

  // ============================================================================
  // REPORT ESTIMATES
  // ============================================================================

  REPORT(theta_0);
  ADREPORT(threshold_raw);
  ADREPORT(discrimination);
  ADREPORT(sigma_theta);

  if (has_random == 1) {
    ADREPORT(sigma_random);
    REPORT(u_random);
  }

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_IRT_POLY_MULTILEVEL_HPP
