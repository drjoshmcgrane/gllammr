// Ordinal (cumulative / adjacent-category / continuation-ratio) model
// with multiple random-effects terms (crossed and/or nested).
//
// Random-effects layout mirrors gllamm_glmm_multi.hpp: Z is one sparse
// [n_obs x q_total] matrix, u is term-major and group-major within term,
// each term with its own (possibly correlated) covariance. The ordinal
// likelihood (links 1-5) mirrors gllamm_ordinal.hpp; the proportional-
// odds-relaxing PPO link remains single-term only.

#ifndef GLLAMM_ORDINAL_MULTI_HPP
#define GLLAMM_ORDINAL_MULTI_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_ordinal_multi(objective_function<Type>* obj)
{
  DATA_IVECTOR(y);               // Ordinal response (1, 2, ..., K)
  DATA_MATRIX(X);                // Fixed effects (intercept dropped)
  DATA_SPARSE_MATRIX(Z);         // Combined RE design [n_obs x q_total]
  DATA_INTEGER(n_obs);
  DATA_INTEGER(n_terms);
  DATA_IVECTOR(term_n_random);   // Coefficients per group, per term
  DATA_IVECTOR(term_n_groups);   // Number of groups, per term
  DATA_IVECTOR(term_correlated); // 1 = correlated coefficients within group
  DATA_INTEGER(n_categories);
  DATA_INTEGER(link);            // 1=logit 2=probit 3=acl 4=crl_fwd 5=crl_bwd
  DATA_VECTOR(weights);

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(u);           // Full RE vector, term-major
  PARAMETER_VECTOR(threshold);   // K-1 thresholds (log-spacing form)
  PARAMETER_VECTOR(log_sigma_u); // SDs, concatenated across terms
  PARAMETER_VECTOR(theta);       // Cholesky parameters, concatenated

  int n_fixed = X.cols();

  // Ordered thresholds via cumulative sum of exponentials
  vector<Type> ordered_threshold(n_categories - 1);
  ordered_threshold(0) = threshold(0);
  for (int k = 1; k < n_categories - 1; k++) {
    ordered_threshold(k) = ordered_threshold(k-1) + exp(threshold(k));
  }

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

  // ---- Linear predictor: one sparse product covers all terms ----
  vector<Type> Zu = Z * u;

  // ---- Ordinal likelihood (links 1-5, as in gllamm_ordinal.hpp) ----
  for (int i = 0; i < n_obs; i++) {
    Type eta = Zu(i);
    for (int j = 0; j < n_fixed; j++) {
      eta += X(i, j) * beta(j);
    }

    int obs_cat = y(i) - 1;
    Type prob_cat;

    if (link == 1 || link == 2) {
      // Cumulative logit / probit
      if (obs_cat == 0) {
        prob_cat = (link == 1)
          ? invlogit(ordered_threshold(0) - eta)
          : pnorm(ordered_threshold(0) - eta);
      } else if (obs_cat == n_categories - 1) {
        prob_cat = (link == 1)
          ? Type(1.0) - invlogit(ordered_threshold(n_categories - 2) - eta)
          : Type(1.0) - pnorm(ordered_threshold(n_categories - 2) - eta);
      } else {
        Type p_le_k, p_le_km1;
        if (link == 1) {
          p_le_k = invlogit(ordered_threshold(obs_cat) - eta);
          p_le_km1 = invlogit(ordered_threshold(obs_cat - 1) - eta);
        } else {
          p_le_k = pnorm(ordered_threshold(obs_cat) - eta);
          p_le_km1 = pnorm(ordered_threshold(obs_cat - 1) - eta);
        }
        prob_cat = p_le_k - p_le_km1;
      }
    } else if (link == 3) {
      // Adjacent-category logit
      vector<Type> log_prob(n_categories);
      log_prob(0) = Type(0.0);
      for (int k = 1; k < n_categories; k++) {
        log_prob(k) = log_prob(k-1) + ordered_threshold(k-1) + eta;
      }
      Type log_sum = log_prob(0);
      for (int k = 1; k < n_categories; k++) {
        log_sum = logspace_add(log_sum, log_prob(k));
      }
      prob_cat = exp(log_prob(obs_cat) - log_sum);
    } else if (link == 4) {
      // Continuation-ratio (forward)
      if (obs_cat == 0) {
        prob_cat = invlogit(ordered_threshold(0) - eta);
      } else if (obs_cat == n_categories - 1) {
        Type surv = Type(1.0);
        for (int j = 0; j < n_categories - 1; j++) {
          surv *= (Type(1.0) - invlogit(ordered_threshold(j) - eta));
        }
        prob_cat = surv;
      } else {
        Type surv = Type(1.0);
        for (int j = 0; j < obs_cat; j++) {
          surv *= (Type(1.0) - invlogit(ordered_threshold(j) - eta));
        }
        prob_cat = surv * invlogit(ordered_threshold(obs_cat) - eta);
      }
    } else {
      // Continuation-ratio (backward), as in gllamm_ordinal.hpp:
      // P(c) = b_c * prod_{j>c} (1 - b_j) with b_c = invlogit(tau_{c-1}-eta)
      Type surv = Type(1.0);
      for (int m = obs_cat; m <= n_categories - 2; m++) {
        surv *= (Type(1.0) - invlogit(ordered_threshold(m) - eta));
      }
      prob_cat = (obs_cat == 0)
        ? surv
        : invlogit(ordered_threshold(obs_cat - 1) - eta) * surv;
    }

    nll -= weights(i) * log(prob_cat + Type(1e-10));
  }

  ADREPORT(ordered_threshold);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif