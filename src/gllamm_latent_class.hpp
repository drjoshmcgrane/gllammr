// Latent Class Analysis
// Finite mixture model with discrete latent classes.
// Item response probabilities are logit-parameterized (unconstrained
// optimization, no boundary penalties) and the mixture is accumulated in
// log space for numerical stability.

#ifndef GLLAMM_LATENT_CLASS_HPP
#define GLLAMM_LATENT_CLASS_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_latent_class(objective_function<Type>* obj)
{
  // Data inputs
  DATA_MATRIX(Y);              // Response matrix (n x p) - items in columns
  DATA_INTEGER(n_obs);         // Number of observations
  DATA_INTEGER(n_items);       // Number of items/variables
  DATA_INTEGER(n_classes);     // Number of latent classes
  DATA_VECTOR(weights);        // Case weights (fweights or pweights)

  // Parameters
  PARAMETER_MATRIX(item_logits);  // logit item response probs (n_items x n_classes)
  PARAMETER_VECTOR(class_logits); // Class membership log-odds (length n_classes - 1)

  // Class probabilities from logits (softmax, last class is reference)
  vector<Type> log_class_probs(n_classes);
  {
    Type sum_exp = Type(1.0);
    for (int k = 0; k < n_classes - 1; k++) {
      sum_exp += exp(class_logits(k));
    }
    Type log_denom = log(sum_exp);
    log_class_probs(n_classes - 1) = -log_denom;
    for (int k = 0; k < n_classes - 1; k++) {
      log_class_probs(k) = class_logits(k) - log_denom;
    }
  }

  // parallel_accumulator splits the likelihood across OpenMP threads
  // when available (no-op on single-threaded builds)
  parallel_accumulator<Type> nll(obj);

  for (int i = 0; i < n_obs; i++) {
    // Log-likelihood of observation i within each class:
    // y*log(p) + (1-y)*log(1-p) with p = invlogit(x) reduces to
    // y*x - log(1 + exp(x)), which is stable for any x
    vector<Type> log_joint(n_classes);
    for (int k = 0; k < n_classes; k++) {
      Type ll_k = Type(0.0);
      for (int j = 0; j < n_items; j++) {
        Type x = item_logits(j, k);
        ll_k += Y(i, j) * x - logspace_add(Type(0.0), x);
      }
      log_joint(k) = log_class_probs(k) + ll_k;
    }

    // log-sum-exp over classes
    Type m = log_joint(0);
    for (int k = 1; k < n_classes; k++) {
      m = logspace_add(m, log_joint(k));
    }

    nll -= weights(i) * m;
  }

  // Report on the probability scale
  matrix<Type> item_probs(n_items, n_classes);
  for (int j = 0; j < n_items; j++) {
    for (int k = 0; k < n_classes; k++) {
      item_probs(j, k) = invlogit(item_logits(j, k));
    }
  }
  vector<Type> class_probs = exp(log_class_probs.array());

  ADREPORT(class_probs);
  ADREPORT(item_probs);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_LATENT_CLASS_HPP
