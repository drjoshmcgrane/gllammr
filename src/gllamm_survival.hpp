// Survival/Time-to-Event Models with Random Effects
// Supports right censoring, Weibull and exponential distributions

#ifndef GLLAMM_SURVIVAL_HPP
#define GLLAMM_SURVIVAL_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type gllamm_survival(objective_function<Type>* obj)
{
  // Data inputs
  DATA_VECTOR(time);           // Observed times
  DATA_IVECTOR(event);         // Event indicator: 1=event, 0=censored
  DATA_MATRIX(X);              // Covariates
  DATA_SPARSE_MATRIX(Z);       // Random effects design
  DATA_IVECTOR(groups);        // Group indices
  DATA_INTEGER(n_groups);      // Number of groups
  DATA_INTEGER(n_obs);         // Number of observations
  DATA_INTEGER(n_fixed);       // Number of fixed effects
  DATA_INTEGER(n_random);      // Number of random effects
  DATA_INTEGER(distribution);  // 1=exponential, 2=Weibull
  DATA_INTEGER(correlated);    // Random effects correlation
  DATA_VECTOR(weights);        // Case weights (fweights or pweights)

  // Parameters
  PARAMETER_VECTOR(beta);      // Fixed effects (log hazard ratios)
  PARAMETER_VECTOR(u);         // Random effects
  PARAMETER(log_shape);        // Log shape parameter (for Weibull)
  PARAMETER_VECTOR(log_sigma_u); // Log random effects SDs
  PARAMETER_VECTOR(theta);     // Cholesky correlation parameters

  // Transform parameters
  Type shape = exp(log_shape);
  vector<Type> sigma_u = exp(log_sigma_u.array());

  // Build variance-covariance matrix
  matrix<Type> Sigma_u(n_random, n_random);

  if (correlated == 1 && n_random > 1) {
    matrix<Type> L(n_random, n_random);
    L.setZero();

    int idx = 0;
    for (int i = 0; i < n_random; i++) {
      L(i, i) = Type(1.0);
      for (int j = 0; j < i; j++) {
        L(i, j) = theta(idx);
        idx++;
      }
    }

    matrix<Type> R = L * L.transpose();

    for (int i = 0; i < n_random; i++) {
      for (int j = 0; j < n_random; j++) {
        Sigma_u(i, j) = sigma_u(i) * sigma_u(j) * R(i, j);
      }
    }
  } else {
    Sigma_u.setZero();
    for (int i = 0; i < n_random; i++) {
      Sigma_u(i, i) = sigma_u(i) * sigma_u(i);
    }
  }

  matrix<Type> Sigma_u_inv = Sigma_u.inverse();
  Type log_det_Sigma_u = atomic::logdet(Sigma_u);

  // Initialize negative log-likelihood
  Type nll = 0.0;

  // Prior for random effects
  for (int j = 0; j < n_groups; j++) {
    vector<Type> u_j(n_random);
    for (int k = 0; k < n_random; k++) {
      u_j(k) = u(j * n_random + k);
    }
    Type quad_form = (u_j * (Sigma_u_inv * u_j)).sum();
    nll += 0.5 * (Type(n_random) * log(2.0 * M_PI) + log_det_Sigma_u + quad_form);
  }

  // Likelihood for survival times
  for (int i = 0; i < n_obs; i++) {
    // Linear predictor
    Type eta = 0.0;
    for (int j = 0; j < n_fixed; j++) {
      eta += X(i, j) * beta(j);
    }

    int g = groups(i);
    for (int k = 0; k < n_random; k++) {
      eta += Z.coeff(i, k) * u(g * n_random + k);
    }

    Type t = time(i);
    int d = event(i);
    Type w_i = weights(i);

    if (distribution == 1) {
      // Exponential: lambda = exp(eta)
      Type lambda = exp(eta);

      if (d == 1) {
        // Event occurred: f(t) = lambda * exp(-lambda * t) (weighted)
        nll -= w_i * (log(lambda) - lambda * t);
      } else {
        // Censored: S(t) = exp(-lambda * t) (weighted)
        nll -= w_i * (-lambda * t);
      }

    } else {
      // Weibull: lambda = exp(eta), shape parameter
      Type lambda = exp(eta);

      if (d == 1) {
        // Event: f(t) = shape * lambda * (lambda*t)^(shape-1) * exp(-(lambda*t)^shape) (weighted)
        nll -= w_i * (log(shape) + log(lambda) + (shape - Type(1.0)) * log(lambda * t) -
               pow(lambda * t, shape));
      } else {
        // Censored: S(t) = exp(-(lambda*t)^shape) (weighted)
        nll -= w_i * (-pow(lambda * t, shape));
      }
    }
  }

  // Report
  ADREPORT(sigma_u);
  if (distribution == 2) {
    ADREPORT(shape);
  }

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif // GLLAMM_SURVIVAL_HPP
