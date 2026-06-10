// Latent Class Analysis
// Finite mixture model with discrete latent classes

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data inputs
  DATA_MATRIX(Y);              // Response matrix (n x p) - items in columns
  DATA_INTEGER(n_obs);         // Number of observations
  DATA_INTEGER(n_items);       // Number of items/variables
  DATA_INTEGER(n_classes);     // Number of latent classes
  DATA_VECTOR(weights);        // Case weights (fweights or pweights)

  // Parameters
  PARAMETER_MATRIX(item_probs); // Item response probabilities (n_items x n_classes)
  PARAMETER_VECTOR(class_logits); // Class membership log-odds (length n_classes - 1)

  // Transform class probabilities from logits (softmax)
  vector<Type> class_probs(n_classes);
  Type sum_exp = Type(1.0); // For reference class
  for (int k = 0; k < n_classes - 1; k++) {
    sum_exp += exp(class_logits(k));
  }
  class_probs(n_classes - 1) = Type(1.0) / sum_exp; // Reference class
  for (int k = 0; k < n_classes - 1; k++) {
    class_probs(k) = exp(class_logits(k)) / sum_exp;
  }

  // Initialize negative log-likelihood
  Type nll = 0.0;

  // For each observation, compute likelihood as mixture
  for (int i = 0; i < n_obs; i++) {
    Type obs_likelihood = 0.0;

    // Sum over latent classes
    for (int k = 0; k < n_classes; k++) {
      Type class_likelihood = Type(1.0);

      // Product over items (conditional independence given class)
      for (int j = 0; j < n_items; j++) {
        Type p = item_probs(j, k);
        Type y_ij = Y(i, j);

        // Bernoulli likelihood for this item
        class_likelihood *= pow(p, y_ij) * pow(Type(1.0) - p, Type(1.0) - y_ij);
      }

      // Weight by class probability
      obs_likelihood += class_probs(k) * class_likelihood;
    }

    // Add to total log-likelihood (weighted)
    Type w_i = weights(i);
    nll -= w_i * log(obs_likelihood + Type(1e-10));
  }

  // Constraints: probabilities must be in (0,1)
  // Apply penalty for probabilities near boundaries
  for (int j = 0; j < n_items; j++) {
    for (int k = 0; k < n_classes; k++) {
      Type p = item_probs(j, k);
      if (p < Type(0.01) || p > Type(0.99)) {
        nll += Type(100.0); // Penalty
      }
    }
  }

  // Report
  ADREPORT(class_probs);
  ADREPORT(item_probs);

  return nll;
}
