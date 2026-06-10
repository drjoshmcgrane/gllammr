// Multi-Level IRT Models: Rasch, 2PL, 3PL with Random Effects
// Supports nested, crossed, and partially nested random effects
//
// Model: theta_i = theta_0i + sum_g u_g[group_g[i]]
//        theta_0i ~ N(0, sigma_theta²)
//        u_g ~ N(0, sigma_g²)

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // ============================================================================
  // DATA INPUTS
  // ============================================================================

  // Standard IRT data
  DATA_VECTOR(y);              // Item responses (0/1)
  DATA_IVECTOR(person_id);     // Person identifier (0-indexed)
  DATA_IVECTOR(item_id);       // Item identifier (0-indexed)
  DATA_INTEGER(n_persons);     // Number of persons
  DATA_INTEGER(n_items);       // Number of items
  DATA_INTEGER(n_obs);         // Number of observations
  DATA_VECTOR(weights);        // Case weights
  DATA_INTEGER(model_type);    // 1=Rasch, 2=2PL, 3=3PL
  DATA_IVECTOR(mc_items);      // 1 if item has guessing, 0 otherwise

  // Multi-level structure
  DATA_INTEGER(has_random);        // 0=standard IRT, 1=multi-level
  DATA_INTEGER(n_random_effects);  // Number of RE levels
  DATA_IMATRIX(group_ids);         // [n_persons × n_random_effects], -1 = NA
  DATA_IVECTOR(n_groups);          // Number of groups per level
  DATA_INTEGER(max_n_groups);      // Maximum groups across levels

  // ============================================================================
  // PARAMETERS
  // ============================================================================

  // Item parameters
  PARAMETER_VECTOR(difficulty);     // Item difficulties
  PARAMETER_VECTOR(discrimination); // Item discriminations (2PL/3PL)
  PARAMETER_VECTOR(guessing);       // Guessing parameters (3PL)

  // Person-level ability deviations
  PARAMETER_VECTOR(theta_0);
  PARAMETER(log_sigma_theta);

  // Random effects (group-level deviations)
  PARAMETER_MATRIX(u_random);       // [max_n_groups × n_random_effects]
  PARAMETER_VECTOR(log_sigma_random); // SD for each RE level

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

  Type nll = 0.0;

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

    // Compute probability based on model type
    Type prob;

    if (model_type == 1) {
      // Rasch model: P(Y=1) = logit^{-1}(theta - b)
      Type eta = theta - difficulty(item);
      prob = invlogit(eta);

    } else if (model_type == 2) {
      // 2PL model: P(Y=1) = logit^{-1}(a*(theta - b))
      Type eta = discrimination(item) * (theta - difficulty(item));
      prob = invlogit(eta);

    } else {
      // 3PL model: P(Y=1) = c + (1-c) * logit^{-1}(a*(theta - b))
      Type eta = discrimination(item) * (theta - difficulty(item));

      if (mc_items(item) == 1) {
        // MC item: apply guessing parameter
        Type c = guessing(item);
        prob = c + (Type(1.0) - c) * invlogit(eta);
      } else {
        // Non-MC item: no guessing (same as 2PL)
        prob = invlogit(eta);
      }
    }

    // Bernoulli log-likelihood (weighted)
    Type w_i = weights(i);
    nll -= w_i * (y(i) * log(prob + Type(1e-10)) +
                   (Type(1.0) - y(i)) * log(Type(1.0) - prob + Type(1e-10)));
  }

  // ============================================================================
  // REPORT ESTIMATES
  // ============================================================================

  ADREPORT(theta_0);
  ADREPORT(difficulty);

  if (model_type >= 2) {
    ADREPORT(discrimination);
  }

  if (model_type == 3) {
    ADREPORT(guessing);
  }

  ADREPORT(sigma_theta);

  if (has_random == 1) {
    ADREPORT(sigma_random);
    ADREPORT(u_random);
  }

  return nll;
}
