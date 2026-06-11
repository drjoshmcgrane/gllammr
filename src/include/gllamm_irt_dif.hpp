// Confirmatory IRT DIF model (IRT-LR DIF; Thissen, Steinberg & Wainer).
//
// Dichotomous Rasch/2PL MML model with
//   - latent regression (impact): theta_p ~ N(z_p' gamma, sigma^2),
//     z_p the person-covariate design (no intercept; the reference
//     profile has latent mean 0)
//   - uniform DIF for studied items: the logit gains z_p' delta_i
//   - nonuniform DIF (2PL): discrimination scaled by exp(z_p' kappa_i)
//
// Anchors (items without DIF parameters) identify DIF against impact.
//
//   eta_pi = a_i exp(z_p' kappa_i) (theta_p - b_i) + z_p' delta_i

#ifndef GLLAMM_IRT_DIF_HPP
#define GLLAMM_IRT_DIF_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_irt_dif(objective_function<Type>* obj)
{
  DATA_VECTOR(y);              // Item responses (0/1), long format
  DATA_IVECTOR(person_id);     // 0-indexed
  DATA_IVECTOR(item_id);       // 0-indexed
  DATA_INTEGER(n_persons);
  DATA_INTEGER(n_items);
  DATA_INTEGER(n_obs);
  DATA_MATRIX(Zp);             // Person covariates [n_persons x q], no intercept
  DATA_IVECTOR(dif_item);      // For each item: row of delta/kappa (0-based),
                               //   or -1 if the item is an anchor
  DATA_INTEGER(model_type);    // 1 = Rasch, 2 = 2PL
  DATA_INTEGER(nonuniform);    // 1 = kappa active (2PL only)

  PARAMETER_VECTOR(theta);       // Person abilities (random)
  PARAMETER_VECTOR(difficulty);  // b_i
  PARAMETER_VECTOR(discrimination); // a_i (mapped off for Rasch)
  PARAMETER(log_sigma_theta);    // free for Rasch; mapped to 0 for 2PL
  PARAMETER_VECTOR(gamma_impact); // latent regression (impact) [q]
  PARAMETER_MATRIX(delta);       // uniform DIF [n_dif x q]
  PARAMETER_MATRIX(kappa);       // nonuniform DIF [n_dif x q]

  int q = Zp.cols();
  Type sigma_theta = exp(log_sigma_theta);

  Type nll = 0.0;

  // Impact: theta_p ~ N(z_p' gamma, sigma^2)
  for (int p = 0; p < n_persons; p++) {
    Type mu_p = 0.0;
    for (int k = 0; k < q; k++) mu_p += Zp(p, k) * gamma_impact(k);
    nll -= dnorm(theta(p), mu_p, sigma_theta, true);
  }

  for (int i = 0; i < n_obs; i++) {
    int p = person_id(i);
    int j = item_id(i);
    int d = dif_item(j);

    Type a = (model_type == 1) ? Type(1.0) : discrimination(j);
    Type shift = 0.0;
    if (d >= 0) {
      if (nonuniform == 1) {
        Type ka = 0.0;
        for (int k = 0; k < q; k++) ka += Zp(p, k) * kappa(d, k);
        a = a * exp(ka);
      }
      for (int k = 0; k < q; k++) shift += Zp(p, k) * delta(d, k);
    }

    Type eta = a * (theta(p) - difficulty(j)) + shift;
    Type prob = invlogit(eta);
    nll -= y(i) * log(prob + Type(1e-10)) +
           (Type(1.0) - y(i)) * log(Type(1.0) - prob + Type(1e-10));
  }

  ADREPORT(gamma_impact);
  if (delta.rows() > 0) {
    ADREPORT(delta);
    if (nonuniform == 1) ADREPORT(kappa);
  }

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif