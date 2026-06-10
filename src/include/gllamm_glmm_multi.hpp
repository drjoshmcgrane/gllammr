// GLMM with multiple random-effects terms (crossed and/or nested),
// gaussian/binomial/poisson families.
//
// lme4-style layout: Z is one sparse [n_obs x q_total] matrix mapping the
// full random-effects vector u to observations; u is laid out term-major,
// and group-major within term:
//   u[offset_t + g * n_random_t + k] = coefficient k of group g in term t.
// Each term has its own (possibly correlated) covariance.

#ifndef GLLAMM_GLMM_MULTI_HPP
#define GLLAMM_GLMM_MULTI_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_glmm_multi(objective_function<Type>* obj)
{
  DATA_VECTOR(y);                // Response vector
  DATA_MATRIX(X);                // Fixed effects design matrix
  DATA_SPARSE_MATRIX(Z);         // Combined RE design matrix [n_obs x q_total]
  DATA_INTEGER(n_obs);
  DATA_INTEGER(n_terms);         // Number of RE terms
  DATA_IVECTOR(term_n_random);   // Coefficients per group, per term
  DATA_IVECTOR(term_n_groups);   // Number of groups, per term
  DATA_IVECTOR(term_correlated); // 1 = correlated coefficients within group
  DATA_INTEGER(family);          // 0 = gaussian, 1 = binomial, 2 = poisson
  DATA_INTEGER(link);            // 1 = canonical, 2 = probit, 3 = cloglog
  DATA_VECTOR(weights);          // Case weights

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(u);           // Full RE vector, term-major
  PARAMETER(log_sigma);          // Residual SD (gaussian; mapped otherwise)
  PARAMETER_VECTOR(log_sigma_u); // SDs, concatenated across terms
  PARAMETER_VECTOR(theta);       // Cholesky parameters, concatenated across terms

  Type sigma = exp(log_sigma);

  Type nll = 0.0;

  // ---- Random-effects priors, term by term ----
  int u_offset = 0;      // position in u
  int sd_offset = 0;     // position in log_sigma_u
  int th_offset = 0;     // position in theta

  for (int t = 0; t < n_terms; t++) {
    int nr = term_n_random(t);
    int ng = term_n_groups(t);

    if (nr == 1) {
      // Scalar random effect: simple normal prior
      Type sd_t = exp(log_sigma_u(sd_offset));
      for (int g = 0; g < ng; g++) {
        nll -= dnorm(u(u_offset + g), Type(0.0), sd_t, true);
      }
    } else {
      // Vector random effect: MVN with (optionally) correlated components
      vector<Type> sigma_t(nr);
      for (int k = 0; k < nr; k++) {
        sigma_t(k) = exp(log_sigma_u(sd_offset + k));
      }

      matrix<Type> Sigma_t(nr, nr);
      if (term_correlated(t) == 1) {
        // Unit-diagonal Cholesky, normalized to a correlation matrix
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
  vector<Type> eta = X * beta + Z * u;

  // ---- Likelihood ----
  vector<Type> fitted(n_obs);
  for (int i = 0; i < n_obs; i++) {
    Type w_i = weights(i);

    if (family == 0) {
      nll -= w_i * dnorm(y(i), eta(i), sigma, true);
      fitted(i) = eta(i);
    } else if (family == 1) {
      Type p;
      if (link == 2) {
        p = pnorm(eta(i));
      } else if (link == 3) {
        p = Type(1.0) - exp(-exp(eta(i)));
      } else {
        p = invlogit(eta(i));
      }
      Type ll_i = y(i) * log(p + Type(1e-10)) +
                  (Type(1.0) - y(i)) * log(Type(1.0) - p + Type(1e-10));
      nll -= w_i * ll_i;
      fitted(i) = p;
    } else if (family == 2) {
      Type lambda = exp(eta(i));
      nll -= w_i * dpois(y(i), lambda, true);
      fitted(i) = lambda;
    } else {
      // Gamma: mean mu, dispersion phi = exp(log_sigma);
      // shape = 1/phi, scale = mu*phi (variance = phi * mu^2)
      Type mu;
      if (link == 2) {
        mu = Type(1.0) / eta(i);
      } else if (link == 3) {
        mu = eta(i);
      } else {
        mu = exp(eta(i));
      }
      Type phi = exp(log_sigma);
      Type shape = Type(1.0) / phi;
      nll -= w_i * dgamma(y(i), shape, mu * phi, true);
      fitted(i) = mu;
    }
  }

  REPORT(fitted);
  ADREPORT(sigma);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_GLMM_MULTI_HPP
