// Nonparametric maximum likelihood (NPML) two-level GLMM:
// the random-intercept distribution is a discrete distribution on K
// estimated mass points. The likelihood is fully marginal (a finite
// mixture over groups) - no Laplace approximation involved.
//
// Identification: X must NOT contain an intercept; the K locations play
// the role of class-specific intercepts.

#ifndef GLLAMM_NPML_HPP
#define GLLAMM_NPML_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_npml(objective_function<Type>* obj)
{
  DATA_VECTOR(y);
  DATA_MATRIX(X);              // No intercept column
  DATA_IVECTOR(groups);        // 0-indexed group per observation
  DATA_INTEGER(n_obs);
  DATA_INTEGER(n_groups);
  DATA_INTEGER(K);             // Number of mass points
  DATA_INTEGER(family);        // 0 = gaussian, 1 = binomial, 2 = poisson
  DATA_VECTOR(weights);        // Observation weights

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(locations);   // K mass-point locations
  PARAMETER_VECTOR(mass_logits); // K-1 free logits (last mass is reference)
  PARAMETER(log_sigma);          // Residual SD (gaussian only)

  Type sigma = exp(log_sigma);

  // Log masses (softmax)
  vector<Type> log_mass(K);
  {
    Type denom = Type(0.0);   // reference logit 0
    for (int k = 0; k < K - 1; k++) {
      denom = logspace_add(denom, mass_logits(k));
    }
    log_mass(K - 1) = -denom;
    for (int k = 0; k < K - 1; k++) {
      log_mass(k) = mass_logits(k) - denom;
    }
  }

  vector<Type> xb = X * beta;

  // Per-group, per-mass-point conditional log-likelihoods
  matrix<Type> group_ll(n_groups, K);
  group_ll.setZero();

  for (int i = 0; i < n_obs; i++) {
    int g = groups(i);
    Type w_i = weights(i);
    for (int k = 0; k < K; k++) {
      Type eta = xb(i) + locations(k);
      Type ll;
      if (family == 0) {
        ll = dnorm(y(i), eta, sigma, true);
      } else if (family == 1) {
        ll = y(i) * eta - logspace_add(Type(0.0), eta);
      } else {
        ll = dpois(y(i), exp(eta), true);
      }
      group_ll(g, k) += w_i * ll;
    }
  }

  // Marginal likelihood: mixture over mass points per group
  Type nll = 0.0;
  for (int g = 0; g < n_groups; g++) {
    Type m = log_mass(0) + group_ll(g, 0);
    for (int k = 1; k < K; k++) {
      m = logspace_add(m, log_mass(k) + group_ll(g, k));
    }
    nll -= m;
  }

  vector<Type> masses = exp(log_mass.array());
  ADREPORT(locations);
  ADREPORT(masses);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_NPML_HPP
