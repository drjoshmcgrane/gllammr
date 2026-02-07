// Explanatory IRT: Item parameters as functions of item covariates
// Supports dichotomous and polytomous models

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // ============================================================================
  // DATA INPUTS
  // ============================================================================

  DATA_VECTOR(y);                    // Item responses
  DATA_IVECTOR(person_id);           // Person identifier (0-indexed)
  DATA_IVECTOR(item_id);             // Item identifier (0-indexed)
  DATA_MATRIX(W_difficulty);         // Item covariate matrix for difficulty [n_items x p_diff]
  DATA_MATRIX(W_discrimination);     // Item covariate matrix for discrimination [n_items x p_disc]
  DATA_INTEGER(n_persons);           // Number of persons
  DATA_INTEGER(n_items);             // Number of items
  DATA_INTEGER(n_obs);               // Number of observations
  DATA_INTEGER(model_type);          // 1=Rasch/1PL, 2=2PL, 3=3PL
  DATA_INTEGER(is_polytomous);       // 0=dichotomous, 1=polytomous

  // For polytomous models
  DATA_IVECTOR(n_categories_per_item);  // Number of categories per item
  DATA_INTEGER(max_categories);          // Maximum K

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

  // For polytomous: threshold parameters
  PARAMETER_MATRIX(threshold_resid); // Threshold residuals [n_items x (max_categories-1)]

  // ============================================================================
  // TRANSFORM PARAMETERS
  // ============================================================================

  Type sigma_epsilon_b = exp(log_sigma_epsilon_b);
  Type sigma_epsilon_a = exp(log_sigma_epsilon_a);
  Type sigma_theta = exp(log_sigma_theta);

  // ============================================================================
  // INITIALIZE NEGATIVE LOG-LIKELIHOOD
  // ============================================================================

  Type nll = 0.0;

  // ============================================================================
  // PRIORS
  // ============================================================================

  // Prior for person abilities: theta ~ N(0, sigma_theta^2)
  for (int p = 0; p < n_persons; p++) {
    nll -= dnorm(theta(p), Type(0.0), sigma_theta, true);
  }

  // Prior for item residuals
  for (int j = 0; j < n_items; j++) {
    nll -= dnorm(epsilon_b(j), Type(0.0), sigma_epsilon_b, true);
    nll -= dnorm(epsilon_a(j), Type(0.0), sigma_epsilon_a, true);
  }

  // ============================================================================
  // COMPUTE ITEM PARAMETERS FROM COVARIATES
  // ============================================================================

  vector<Type> difficulty(n_items);
  vector<Type> discrimination(n_items);

  int p_diff = W_difficulty.cols();
  int p_disc = W_discrimination.cols();

  for (int j = 0; j < n_items; j++) {
    // Difficulty as linear function of item covariates
    Type difficulty_pred = 0.0;
    for (int p = 0; p < p_diff; p++) {
      difficulty_pred += gamma(p) * W_difficulty(j, p);
    }
    difficulty(j) = difficulty_pred + epsilon_b(j);

    // Log-discrimination as linear function of item covariates
    Type log_discrim_pred = 0.0;
    for (int p = 0; p < p_disc; p++) {
      log_discrim_pred += delta(p) * W_discrimination(j, p);
    }
    discrimination(j) = exp(log_discrim_pred + epsilon_a(j));
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
        Type eta = theta(person) - difficulty(item);
        prob = invlogit(eta);

      } else if (model_type == 2) {
        // 2PL model: P(Y=1) = logit^{-1}(a*(theta - b))
        Type eta = discrimination(item) * (theta(person) - difficulty(item));
        prob = invlogit(eta);

      } else {
        // 3PL model: P(Y=1) = c + (1-c) * logit^{-1}(a*(theta - b))
        // For simplicity, guessing parameter not included in EIRT
        Type eta = discrimination(item) * (theta(person) - difficulty(item));
        prob = invlogit(eta);
      }

      // Bernoulli log-likelihood
      nll -= y(i) * log(prob + Type(1e-10)) + (Type(1.0) - y(i)) * log(Type(1.0) - prob + Type(1e-10));
    }

  } else {
    // ========================================================================
    // POLYTOMOUS IRT MODELS (GRM-like with covariates)
    // ========================================================================

    for (int i = 0; i < n_obs; i++) {
      int person = person_id(i);
      int item = item_id(i);
      int obs_cat = CppAD::Integer(y(i)) - 1;  // Convert to 0-indexed
      int K = n_categories_per_item(item);

      // Build ordered thresholds with covariate adjustments
      vector<Type> ordered_threshold(K - 1);

      // Base threshold from difficulty + residuals
      ordered_threshold(0) = difficulty(item) + threshold_resid(item, 0);

      for (int k = 1; k < K - 1; k++) {
        ordered_threshold(k) = ordered_threshold(k-1) + exp(threshold_resid(item, k));
      }

      // Compute probability using GRM-like formulation
      Type prob_cat;

      if (obs_cat == 0) {
        prob_cat = invlogit(discrimination(item) * (theta(person) - ordered_threshold(0)));
      } else if (obs_cat == K - 1) {
        prob_cat = Type(1.0) - invlogit(discrimination(item) * (theta(person) - ordered_threshold(K - 2)));
      } else {
        Type p_le_k = invlogit(discrimination(item) * (theta(person) - ordered_threshold(obs_cat)));
        Type p_le_k_minus_1 = invlogit(discrimination(item) * (theta(person) - ordered_threshold(obs_cat - 1)));
        prob_cat = p_le_k - p_le_k_minus_1;
      }

      nll -= log(prob_cat + Type(1e-10));
    }
  }

  // ============================================================================
  // REPORT ESTIMATES
  // ============================================================================

  ADREPORT(theta);
  ADREPORT(gamma);
  ADREPORT(delta);
  ADREPORT(difficulty);
  ADREPORT(discrimination);
  ADREPORT(sigma_epsilon_b);
  ADREPORT(sigma_epsilon_a);
  ADREPORT(sigma_theta);

  return nll;
}
