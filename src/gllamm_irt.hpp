// IRT Models: Rasch, 2PL, 3PL
// Item Response Theory with latent ability as random effect

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data inputs
  DATA_VECTOR(y);              // Item responses (0/1)
  DATA_IVECTOR(person_id);     // Person identifier (0-indexed)
  DATA_IVECTOR(item_id);       // Item identifier (0-indexed)
  DATA_INTEGER(n_persons);     // Number of persons
  DATA_INTEGER(n_items);       // Number of items
  DATA_INTEGER(n_obs);         // Number of observations
  DATA_INTEGER(model_type);    // 1=Rasch, 2=2PL, 3=3PL

  // Parameters
  PARAMETER_VECTOR(theta);     // Person abilities (latent trait)
  PARAMETER_VECTOR(difficulty); // Item difficulties (b parameters)
  PARAMETER_VECTOR(discrimination); // Item discriminations (a parameters, for 2PL/3PL)
  PARAMETER_VECTOR(guessing);  // Guessing parameters (c parameters, for 3PL)
  PARAMETER(log_sigma_theta);  // Log SD of ability distribution

  // Transform parameters
  Type sigma_theta = exp(log_sigma_theta);

  // Initialize negative log-likelihood
  Type nll = 0.0;

  // Prior for person abilities: theta ~ N(0, sigma_theta^2)
  for (int p = 0; p < n_persons; p++) {
    nll -= dnorm(theta(p), Type(0.0), sigma_theta, true);
  }

  // Likelihood for item responses
  for (int i = 0; i < n_obs; i++) {
    int person = person_id(i);
    int item = item_id(i);

    Type prob;

    if (model_type == 1) {
      // Rasch model: P(Y=1) = logit^{-1}(theta - b)
      Type eta = theta(person) - difficulty(item);
      prob = invlogit(eta);

    } else if (model_type == 2) {
      // 2PL model: P(Y=1) = logit^{-1}(a*(theta - b))
      Type eta = discrimination(item) * (theta(person) - difficulty(item));
      prob = invlogit(eta);

    } else {
      // 3PL model: P(Y=1) = c + (1-c) * logit^{-1}(a*(theta - b))
      Type eta = discrimination(item) * (theta(person) - difficulty(item));
      Type c = guessing(item);
      prob = c + (Type(1.0) - c) * invlogit(eta);
    }

    // Bernoulli log-likelihood
    nll -= y(i) * log(prob + Type(1e-10)) + (Type(1.0) - y(i)) * log(Type(1.0) - prob + Type(1e-10));
  }

  // Report estimates
  ADREPORT(theta);
  ADREPORT(difficulty);
  if (model_type >= 2) {
    ADREPORT(discrimination);
  }
  if (model_type == 3) {
    ADREPORT(guessing);
  }
  ADREPORT(sigma_theta);

  return nll;
}
