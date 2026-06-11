// Multinomial (baseline-category logit) model with multiple random-effects
// terms (crossed and/or nested).
//
// Random-effects layout mirrors gllamm_glmm_multi.hpp: Z is one sparse
// [n_obs x q_total] matrix, u term-major and group-major within term, each
// term with its own (possibly correlated) covariance. As in the
// single-term template (gllamm_multinomial.hpp), the random effects act as
// a common shifter added to every non-reference category's linear
// predictor.

#ifndef GLLAMM_MULTINOMIAL_MULTI_HPP
#define GLLAMM_MULTINOMIAL_MULTI_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_multinomial_multi(objective_function<Type>* obj)
{
  DATA_IVECTOR(y);               // Nominal response (0, 1, ..., K-1)
  DATA_MATRIX(X);                // Fixed effects design matrix
  DATA_SPARSE_MATRIX(Z);         // Combined RE design [n_obs x q_total]
  DATA_INTEGER(n_obs);
  DATA_INTEGER(n_terms);
  DATA_IVECTOR(term_n_random);
  DATA_IVECTOR(term_n_groups);
  DATA_IVECTOR(term_correlated);
  DATA_INTEGER(n_fixed);         // Fixed effects per category
  DATA_INTEGER(n_categories);
  DATA_VECTOR(weights);

  PARAMETER_MATRIX(beta);        // (n_categories-1) x n_fixed
  PARAMETER_VECTOR(u);           // Full RE vector, term-major
  PARAMETER_VECTOR(log_sigma_u); // SDs, concatenated across terms
  PARAMETER_VECTOR(theta);       // Cholesky parameters, concatenated

  Type nll = 0.0;

  // ---- Random-effects priors, term by term (as in glmm_multi) ----
  int u_offset = 0, sd_offset = 0, th_offset = 0;
  for (int t = 0; t < n_terms; t++) {
    int nr = term_n_random(t);
    int ng = term_n_groups(t);

    if (nr == 1) {
      Type sd_t = exp(log_sigma_u(sd_offset));
      for (int g = 0; g < ng; g++) {
        nll -= dnorm(u(u_offset + g), Type(0.0), sd_t, true);
      }
    } else {
      vector<Type> sigma_t(nr);
      for (int k = 0; k < nr; k++) {
        sigma_t(k) = exp(log_sigma_u(sd_offset + k));
      }
      matrix<Type> Sigma_t(nr, nr);
      if (term_correlated(t) == 1) {
        matrix<Type> L(nr, nr);
        L.setZero();
        int idx = th_offset;
        for (int i = 0; i < nr; i++) {
          L(i, i) = Type(1.0);
          for (int j = 0; j < i; j++) {
            L(i, j) = theta(idx);
            idx++;
          }
        }
        matrix<Type> R = L * L.transpose();
        for (int i = 0; i < nr; i++) {
          for (int j = 0; j < nr; j++) {
            Type rij = R(i, j) / sqrt(R(i, i) * R(j, j));
            Sigma_t(i, j) = sigma_t(i) * sigma_t(j) * rij;
          }
        }
      } else {
        Sigma_t.setZero();
        for (int k = 0; k < nr; k++) {
          Sigma_t(k, k) = sigma_t(k) * sigma_t(k);
        }
      }
      matrix<Type> Sigma_inv = atomic::matinv(Sigma_t);
      Type log_det = atomic::logdet(Sigma_t);
      for (int g = 0; g < ng; g++) {
        vector<Type> u_g = u.segment(u_offset + g * nr, nr);
        Type quad = (u_g * (Sigma_inv * u_g)).sum();
        nll += 0.5 * (Type(nr) * log(2.0 * M_PI) + log_det + quad);
      }
    }

    u_offset += ng * nr;
    sd_offset += nr;
    if (nr > 1 && term_correlated(t) == 1) {
      th_offset += nr * (nr - 1) / 2;
    }
  }

  // ---- Common random shifter: one sparse product covers all terms ----
  vector<Type> Zu = Z * u;

  // ---- Baseline-category logit likelihood ----
  for (int i = 0; i < n_obs; i++) {
    int obs_cat = y(i);

    vector<Type> eta(n_categories);
    eta.setZero();
    for (int cat = 1; cat < n_categories; cat++) {
      for (int j = 0; j < n_fixed; j++) {
        eta(cat) += X(i, j) * beta(cat - 1, j);
      }
      eta(cat) += Zu(i);
    }

    Type sum_exp = Type(1.0);
    for (int cat = 1; cat < n_categories; cat++) {
      sum_exp += exp(eta(cat));
    }

    Type prob_obs_cat = (obs_cat == 0)
      ? Type(1.0) / sum_exp
      : exp(eta(obs_cat)) / sum_exp;

    nll -= weights(i) * log(prob_obs_cat + Type(1e-10));
  }

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif