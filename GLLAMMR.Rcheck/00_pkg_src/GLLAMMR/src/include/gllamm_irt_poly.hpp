// Polytomous IRT Models: GRM, PCM, GPCM, NRM
// Template for ordered and unordered categorical item responses

#ifndef GLLAMM_IRT_POLY_HPP
#define GLLAMM_IRT_POLY_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_irt_poly(objective_function<Type>* obj)
{
  // ============================================================================
  // DATA INPUTS
  // ============================================================================

  DATA_VECTOR(y);                     // Item responses (1, 2, ..., K)
  DATA_IVECTOR(person_id);            // Person identifier (0-indexed)
  DATA_IVECTOR(item_id);              // Item identifier (0-indexed)
  DATA_IVECTOR(n_categories_per_item); // Number of categories per item
  DATA_INTEGER(max_categories);       // Maximum K across all items
  DATA_INTEGER(n_persons);            // Number of persons
  DATA_INTEGER(n_items);              // Number of items
  DATA_INTEGER(n_obs);                // Number of observations
  DATA_VECTOR(weights);               // Case weights (fweights or pweights)
  DATA_INTEGER(model_type);           // 1=GRM, 2=PCM, 3=GPCM, 4=NRM

  // ============================================================================
  // PARAMETERS
  // ============================================================================

  PARAMETER_VECTOR(theta);            // Person abilities (latent trait)
  PARAMETER_MATRIX(threshold_raw);    // Raw threshold parameters [n_items x (max_categories-1)]
  PARAMETER_VECTOR(discrimination);   // Item discriminations (for GRM, GPCM, NRM)
  PARAMETER(log_sigma_theta);         // Log SD of ability distribution

  // ============================================================================
  // TRANSFORM PARAMETERS
  // ============================================================================

  Type sigma_theta = exp(log_sigma_theta);

  // ============================================================================
  // INITIALIZE NEGATIVE LOG-LIKELIHOOD
  // ============================================================================

  // parallel_accumulator splits the likelihood across OpenMP threads
  // when available (no-op on single-threaded builds)
  parallel_accumulator<Type> nll(obj);

  // ============================================================================
  // PRIOR FOR PERSON ABILITIES
  // ============================================================================

  // Prior: theta ~ N(0, sigma_theta^2)
  for (int p = 0; p < n_persons; p++) {
    nll -= dnorm(theta(p), Type(0.0), sigma_theta, true);
  }

  // ============================================================================
  // LIKELIHOOD FOR ITEM RESPONSES
  // ============================================================================

  for (int i = 0; i < n_obs; i++) {
    int person = person_id(i);
    int item = item_id(i);
    int obs_cat = CppAD::Integer(y(i)) - 1;  // Convert to 0-indexed
    int K = n_categories_per_item(item);      // Number of categories for this item

    Type prob_cat;

    // ========================================================================
    // MODEL 1: GRADED RESPONSE MODEL (GRM)
    // ========================================================================
    if (model_type == 1) {
      // Enforce ordered thresholds using cumulative exponential
      vector<Type> ordered_threshold(K - 1);
      ordered_threshold(0) = threshold_raw(item, 0);
      for (int k = 1; k < K - 1; k++) {
        ordered_threshold(k) = ordered_threshold(k-1) + exp(threshold_raw(item, k));
      }

      // Get discrimination for this item
      Type a = discrimination(item);

      // Compute probability for observed category
      // P(Y = k) = P(Y >= k) - P(Y >= k+1)
      // where P(Y >= k) = invlogit(a * (theta - tau_k))

      // Standard IRT GRM: P(Y >= k) = invlogit(a*(theta - tau_{k-1}))
      // P(Y = k) = P(Y >= k) - P(Y >= k+1), requires ordered thresholds
      if (obs_cat == 0) {
        // P(Y = 0) = 1 - P(Y >= 1) = 1 - invlogit(a*(theta - tau_0))
        prob_cat = Type(1.0) - invlogit(a * (theta(person) - ordered_threshold(0)));

      } else if (obs_cat == K - 1) {
        // P(Y = K-1) = P(Y >= K-1) = invlogit(a*(theta - tau_{K-2}))
        prob_cat = invlogit(a * (theta(person) - ordered_threshold(K - 2)));

      } else {
        // P(Y = k) = invlogit(a*(theta - tau_{k-1})) - invlogit(a*(theta - tau_k))
        Type p_ge_k = invlogit(a * (theta(person) - ordered_threshold(obs_cat - 1)));
        Type p_ge_k_plus_1 = invlogit(a * (theta(person) - ordered_threshold(obs_cat)));
        prob_cat = p_ge_k - p_ge_k_plus_1;
      }
    }

    // ========================================================================
    // MODEL 2: PARTIAL CREDIT MODEL (PCM) - Rasch for polytomous
    // ========================================================================
    else if (model_type == 2) {
      // PCM uses adjacent-categories logits: log[P(Y=m)/P(Y=m-1)] = theta - delta_m
      // Step difficulties delta_m are FREE (no ordering constraint).
      // P(Y = k) = exp(sum_{m=1}^k (theta - delta_m)) / sum_j exp(sum_{m=1}^j (theta - delta_m))

      vector<Type> cumsum(K);
      cumsum(0) = Type(0.0);

      for (int m = 1; m < K; m++) {
        // threshold_raw(item, m-1) = delta_{i,m} directly (free step difficulty)
        cumsum(m) = cumsum(m-1) + (theta(person) - threshold_raw(item, m - 1));
      }

      Type denom = Type(0.0);
      for (int m = 0; m < K; m++) {
        denom += exp(cumsum(m));
      }

      prob_cat = exp(cumsum(obs_cat)) / denom;
    }

    // ========================================================================
    // MODEL 3: GENERALIZED PARTIAL CREDIT MODEL (GPCM)
    // ========================================================================
    else if (model_type == 3) {
      // GPCM: adjacent-categories logit with item-specific discrimination
      // log[P(Y=m)/P(Y=m-1)] = a * (theta - delta_m)
      // Step difficulties delta_m are FREE (no ordering constraint).
      Type a = discrimination(item);

      vector<Type> cumsum(K);
      cumsum(0) = Type(0.0);

      for (int m = 1; m < K; m++) {
        // threshold_raw(item, m-1) = delta_{i,m} directly (free step difficulty)
        cumsum(m) = cumsum(m-1) + a * (theta(person) - threshold_raw(item, m - 1));
      }

      Type denom = Type(0.0);
      for (int m = 0; m < K; m++) {
        denom += exp(cumsum(m));
      }

      prob_cat = exp(cumsum(obs_cat)) / denom;
    }

    // ========================================================================
    // MODEL 4: NOMINAL RESPONSE MODEL (NRM)
    // ========================================================================
    else if (model_type == 4) {
      // NRM: No ordering constraints
      // Each category has its own slope and intercept
      // P(Y = k) = exp(a_k * theta + c_k) / sum_j exp(a_j * theta + c_j)

      // Use threshold_raw for intercepts, discrimination for slopes
      // For identification: set first category to 0
      vector<Type> eta(K);
      eta(0) = Type(0.0);  // Reference category

      for (int k = 1; k < K; k++) {
        Type a_k = discrimination(item);  // Item discrimination (shared or category-specific)
        Type c_k = threshold_raw(item, k - 1);  // Category intercept
        eta(k) = a_k * theta(person) + c_k;
      }

      // Compute denominator
      Type denom = Type(0.0);
      for (int k = 0; k < K; k++) {
        denom += exp(eta(k));
      }

      // Probability for observed category
      prob_cat = exp(eta(obs_cat)) / denom;
    }

    // ========================================================================
    // ADD TO NEGATIVE LOG-LIKELIHOOD
    // ========================================================================

    // Categorical log-likelihood (weighted)
    Type w_i = weights(i);
    nll -= w_i * log(prob_cat + Type(1e-10));
  }

  // ============================================================================
  // REPORT ESTIMATES
  // ============================================================================

  REPORT(theta);
  ADREPORT(discrimination);
  ADREPORT(sigma_theta);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_IRT_POLY_HPP
