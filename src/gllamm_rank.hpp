// Rank-ordered (exploded) logit with group-level taste heterogeneity.
//
// Rows are alternatives, sorted by (case, rank) with rank 1 = most
// preferred; ties are not supported. The ranking likelihood explodes into
// sequential conditional logits: at stage s the alternative ranked s is
// chosen from all alternatives ranked >= s. Unranked alternatives appear
// in every choice set but contribute no stage of their own.
//
// A constant-within-case random intercept cancels in conditional logits,
// so the random effect enters as a random coefficient: utility gains
// u_g * Zu(row), where Zu is an alternative-varying attribute.

#ifndef GLLAMM_RANK_HPP
#define GLLAMM_RANK_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_rank(objective_function<Type>* obj)
{
  DATA_MATRIX(X);              // Alternative-level covariates (n_rows x p)
  DATA_VECTOR(Zu);             // Attribute carrying the random coefficient
  DATA_IVECTOR(case_start);    // First row of each case (0-indexed)
  DATA_IVECTOR(case_n_alts);   // Alternatives per case
  DATA_IVECTOR(case_n_ranked); // Ranked alternatives per case
  DATA_IVECTOR(case_group);    // Group index per case (0-indexed)
  DATA_INTEGER(n_cases);
  DATA_INTEGER(n_groups);
  DATA_VECTOR(case_weights);   // One weight per case

  PARAMETER_VECTOR(beta);      // Preference coefficients
  PARAMETER_VECTOR(u);         // Group random coefficients (scalar per group)
  PARAMETER(log_sigma_u);

  Type sigma_u = exp(log_sigma_u);

  Type nll = 0.0;

  // Random-effects prior
  for (int g = 0; g < n_groups; g++) {
    nll -= dnorm(u(g), Type(0.0), sigma_u, true);
  }

  // Utilities
  vector<Type> eta = X * beta;

  for (int c = 0; c < n_cases; c++) {
    int start = case_start(c);
    int n_alt = case_n_alts(c);
    int n_ranked = case_n_ranked(c);
    int g = case_group(c);
    Type w_c = case_weights(c);

    for (int a = 0; a < n_alt; a++) {
      eta(start + a) += Zu(start + a) * u(g);
    }

    // Exploded logit: stage s chooses row (start + s) among rows
    // start+s .. start+n_alt-1 (rows sorted by rank within case).
    // A full ranking's last stage is uninformative (one alternative left).
    int n_stages = (n_ranked < n_alt) ? n_ranked : (n_ranked - 1);
    for (int s = 0; s < n_stages; s++) {
      Type denom = eta(start + s);
      for (int a = s + 1; a < n_alt; a++) {
        denom = logspace_add(denom, eta(start + a));
      }
      nll -= w_c * (eta(start + s) - denom);
    }
  }

  ADREPORT(sigma_u);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_RANK_HPP
