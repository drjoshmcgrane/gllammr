// Two-level GLMM estimated by adaptive Gauss-Hermite quadrature.
//
// The scalar random intercept is integrated by Gauss-Hermite quadrature
// with group-specific centers m_j and scales s_j (supplied as data and
// updated between optimization rounds by the R driver - the classic
// adaptive scheme). With u = m_j + sqrt(2) s_j x:
//   log L_j = log(sqrt(2) s_j)
//           + logsumexp_q [ log w_q + x_q^2
//                           + log phi(u_q; 0, sigma_u)
//                           + sum_{i in j} loglik_i(eta_i + u_q) ]
// No Laplace approximation is involved; group-level weights multiply the
// whole-group log marginal (exact pseudo-likelihood weighting).

#ifndef GLLAMM_GLMM_AGHQ_HPP
#define GLLAMM_GLMM_AGHQ_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_glmm_aghq(objective_function<Type>* obj)
{
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_IVECTOR(groups);        // 0-indexed
  DATA_INTEGER(n_obs);
  DATA_INTEGER(n_groups);
  DATA_INTEGER(family);        // 0 = gaussian, 1 = binomial, 2 = poisson
  DATA_INTEGER(link);          // 1 = canonical, 2 = probit, 3 = cloglog
  DATA_VECTOR(weights);        // Level-1 weights
  DATA_VECTOR(group_weights);  // Level-2 weights (applied OUTSIDE the integral)
  DATA_VECTOR(gh_x);           // Gauss-Hermite nodes
  DATA_VECTOR(gh_logw);        // log(w_q) + x_q^2 (precomputed)
  DATA_VECTOR(center);         // m_j: adaptation centers per group
  DATA_VECTOR(scale);          // s_j: adaptation scales per group

  PARAMETER_VECTOR(beta);
  PARAMETER(log_sigma);        // Residual SD (gaussian; mapped otherwise)
  PARAMETER(log_sigma_u);

  Type sigma = exp(log_sigma);
  Type sigma_u = exp(log_sigma_u);
  int Q = gh_x.size();

  vector<Type> xb = X * beta;

  // Per-group, per-node accumulated conditional log-likelihoods
  matrix<Type> node_ll(n_groups, Q);
  node_ll.setZero();

  Type sqrt2 = sqrt(Type(2.0));

  for (int i = 0; i < n_obs; i++) {
    int g = groups(i);
    Type w_i = weights(i);
    for (int q = 0; q < Q; q++) {
      Type u_q = center(g) + sqrt2 * scale(g) * gh_x(q);
      Type eta = xb(i) + u_q;
      Type ll;
      if (family == 0) {
        ll = dnorm(y(i), eta, sigma, true);
      } else if (family == 1) {
        if (link == 2) {
          Type p = pnorm(eta);
          ll = y(i) * log(p + Type(1e-12)) +
               (Type(1.0) - y(i)) * log(Type(1.0) - p + Type(1e-12));
        } else if (link == 3) {
          Type p = Type(1.0) - exp(-exp(eta));
          ll = y(i) * log(p + Type(1e-12)) +
               (Type(1.0) - y(i)) * log(Type(1.0) - p + Type(1e-12));
        } else {
          ll = y(i) * eta - logspace_add(Type(0.0), eta);
        }
      } else {
        ll = dpois(y(i), exp(eta), true);
      }
      node_ll(g, q) += w_i * ll;
    }
  }

  Type nll = 0.0;
  for (int g = 0; g < n_groups; g++) {
    Type m = Type(-1e30);
    for (int q = 0; q < Q; q++) {
      Type u_q = center(g) + sqrt2 * scale(g) * gh_x(q);
      Type term = gh_logw(q) + dnorm(u_q, Type(0.0), sigma_u, true) +
                  node_ll(g, q);
      m = logspace_add(m, term);
    }
    Type log_Lg = log(sqrt2 * scale(g)) + m;
    nll -= group_weights(g) * log_Lg;
  }

  ADREPORT(sigma_u);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_GLMM_AGHQ_HPP
