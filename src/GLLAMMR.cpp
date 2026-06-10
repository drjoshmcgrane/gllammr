// GLLAMMR: single TMB objective function dispatching to all model templates.
// Each header defines one model function gllamm_<name>(objective_function<Type>*)
// using the TMB_OBJECTIVE_PTR convention, so DATA_/PARAMETER_ macros resolve
// against the passed objective. Exactly one operator() may exist per DLL.
#include <TMB.hpp>

#include "include/gllamm_gaussian.hpp"
#include "include/gllamm_glmm_slopes.hpp"
#include "include/gllamm_glmm_multi.hpp"
#include "include/gllamm_binomial.hpp"
#include "include/gllamm_poisson.hpp"
#include "include/gllamm_ordinal.hpp"
#include "include/gllamm_multinomial.hpp"
#include "include/gllamm_irt.hpp"
#include "include/gllamm_irt_multilevel.hpp"
#include "include/gllamm_irt_poly.hpp"
#include "include/gllamm_irt_poly_multilevel.hpp"
#include "include/gllamm_eirt.hpp"
#include "include/gllamm_latent_class.hpp"
#include "include/gllamm_mixed_response.hpp"
#include "include/gllamm_sem.hpp"
#include "include/gllamm_survival.hpp"
#include "include/gllamm_rank.hpp"
#include "include/gllamm_npml.hpp"
#include "include/gllamm_glmm_aghq.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(model_name);

  if (model_name == "gaussian")                 return gllamm_gaussian(this);
  else if (model_name == "glmm_slopes")         return gllamm_glmm_slopes(this);
  else if (model_name == "glmm_multi")          return gllamm_glmm_multi(this);
  else if (model_name == "binomial")            return gllamm_binomial(this);
  else if (model_name == "poisson")             return gllamm_poisson(this);
  else if (model_name == "ordinal")             return gllamm_ordinal(this);
  else if (model_name == "multinomial")         return gllamm_multinomial(this);
  else if (model_name == "irt")                 return gllamm_irt(this);
  else if (model_name == "irt_multilevel")      return gllamm_irt_multilevel(this);
  else if (model_name == "irt_poly")            return gllamm_irt_poly(this);
  else if (model_name == "irt_poly_multilevel") return gllamm_irt_poly_multilevel(this);
  else if (model_name == "eirt")                return gllamm_eirt(this);
  else if (model_name == "latent_class")        return gllamm_latent_class(this);
  else if (model_name == "mixed_response")      return gllamm_mixed_response(this);
  else if (model_name == "sem")                 return gllamm_sem(this);
  else if (model_name == "survival")            return gllamm_survival(this);
  else if (model_name == "rank")                return gllamm_rank(this);
  else if (model_name == "npml")                return gllamm_npml(this);
  else if (model_name == "glmm_aghq")           return gllamm_glmm_aghq(this);
  else Rf_error("Unknown model_name '%s'", model_name.c_str());

  return Type(0);
}
