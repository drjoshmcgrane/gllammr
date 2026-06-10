// Latent Class Analysis with mixed indicator types
// item_type: 0 = binary (Bernoulli), 1 = categorical (multinomial, coded
// 1..K), 2 = continuous (gaussian). All parameters are unconstrained
// (logits / log-SDs) and the mixture is accumulated in log space.

#ifndef GLLAMM_LATENT_CLASS_HPP
#define GLLAMM_LATENT_CLASS_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_latent_class(objective_function<Type>* obj)
{
  // Data inputs
  DATA_MATRIX(Y);              // Response matrix (n x p); categorical coded 1..K
  DATA_INTEGER(n_obs);
  DATA_INTEGER(n_items);
  DATA_INTEGER(n_classes);
  DATA_IVECTOR(item_type);     // 0 = binary, 1 = categorical, 2 = gaussian
  DATA_IVECTOR(n_cats);        // categories per categorical item (0 otherwise)
  DATA_INTEGER(max_cats);      // max categories across categorical items (>= 2)
  DATA_VECTOR(weights);        // Case weights

  // Parameters
  PARAMETER_MATRIX(item_logits);   // binary: logit P(Y=1) (n_items x n_classes)
  PARAMETER_ARRAY(cat_logits);     // categorical: logits vs reference category
                                   // (n_items x n_classes x (max_cats - 1))
  PARAMETER_MATRIX(item_means);    // gaussian means (n_items x n_classes)
  PARAMETER_MATRIX(item_log_sds);  // gaussian log-SDs (n_items x n_classes)
  PARAMETER_VECTOR(class_logits);  // class log-odds (n_classes - 1)

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
    vector<Type> log_joint(n_classes);
    for (int k = 0; k < n_classes; k++) {
      Type ll_k = Type(0.0);
      for (int j = 0; j < n_items; j++) {
        if (item_type(j) == 0) {
          // Binary: y*x - log(1 + exp(x)) with x = logit p
          Type x = item_logits(j, k);
          ll_k += Y(i, j) * x - logspace_add(Type(0.0), x);
        } else if (item_type(j) == 1) {
          // Categorical: softmax with category 1 as reference
          int K = n_cats(j);
          int y_cat = CppAD::Integer(Y(i, j)) - 1;   // 0-based
          Type log_denom = Type(0.0);                // ref category logit 0
          for (int m = 0; m < K - 1; m++) {
            log_denom = logspace_add(log_denom, cat_logits(j, k, m));
          }
          Type eta_y = (y_cat == 0) ? Type(0.0) : cat_logits(j, k, y_cat - 1);
          ll_k += eta_y - log_denom;
        } else {
          // Gaussian
          Type sd_jk = exp(item_log_sds(j, k));
          ll_k += dnorm(Y(i, j), item_means(j, k), sd_jk, true);
        }
      }
      log_joint(k) = log_class_probs(k) + ll_k;
    }

    Type m = log_joint(0);
    for (int k = 1; k < n_classes; k++) {
      m = logspace_add(m, log_joint(k));
    }
    nll -= weights(i) * m;
  }

  // Report on natural scales
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
