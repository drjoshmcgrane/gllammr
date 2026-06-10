#' Fit Generalized Linear Latent and Mixed Models
#'
#' Main function for fitting GLLAMM models. Supports multilevel generalized
#' linear models, factor models, IRT models, latent class models, and more.
#'
#' @param formula A two-sided formula with syntax: \code{y ~ x + (terms | group)}.
#'   The right-hand side contains fixed effects and random effects specifications.
#'   Random effects are specified using \code{(term | group)} for correlated
#'   random effects or \code{(term || group)} for uncorrelated random effects.
#'   Nested random effects can be specified as \code{(1 | level1/level2)}.
#'
#' @param data A data frame containing the variables in the formula.
#'
#' @param family A GLM family object specifying the error distribution and link
#'   function. Supports:
#'   \itemize{
#'     \item \code{gaussian()} - Normal distribution (default)
#'     \item \code{binomial()} - Binary/binomial outcomes
#'     \item \code{poisson()} - Count data
#'     \item \code{ordinal()} - Ordered categorical responses (see \code{?ordinal})
#'   }
#'
#' @param start Optional named list of starting values for parameters.
#'
#' @param control A list of control parameters for the optimization algorithm:
#'   \describe{
#'     \item{eval.max}{Maximum number of function evaluations (default: 1000)}
#'     \item{iter.max}{Maximum number of iterations (default: 500)}
#'     \item{trace}{Integer controlling printed output (default: 0)}
#'   }
#'
#' @param ... Additional arguments (reserved for future use).
#'
#' @return An object of class \code{gllamm} with components:
#'   \item{coefficients}{List with \code{fixed} (fixed effects coefficients)
#'     and \code{random_var} (random effects variance components)}
#'   \item{vcov}{List with \code{fixed} (variance-covariance matrix of fixed
#'     effects) and \code{all} (full variance-covariance matrix)}
#'   \item{random_effects}{List of random effects predictions by group}
#'   \item{fitted.values}{Fitted values on the response scale}
#'   \item{residuals}{Response residuals}
#'   \item{y}{Response vector}
#'   \item{X}{Fixed effects design matrix}
#'   \item{logLik}{Log-likelihood at convergence}
#'   \item{AIC}{Akaike Information Criterion}
#'   \item{BIC}{Bayesian Information Criterion}
#'   \item{n_obs}{Number of observations}
#'   \item{n_params}{Number of parameters}
#'   \item{n_groups}{Number of groups for each random effects term}
#'   \item{convergence}{List with convergence information}
#'   \item{call}{The matched call}
#'   \item{formula}{The model formula}
#'   \item{family}{The GLM family}
#'   \item{data}{The data frame (if requested)}
#'   \item{random_terms}{List of random effects specifications}
#'
#' @details
#' GLLAMMR uses Template Model Builder (TMB) for efficient computation via
#' automatic differentiation. The random effects are integrated out using
#' Laplace approximation, providing fast and accurate inference.
#'
#' The formula syntax follows lme4 conventions:
#' \itemize{
#'   \item \code{(1 | group)} - Random intercept for group
#'   \item \code{(x | group)} - Random intercept and slope for x
#'   \item \code{(x || group)} - Uncorrelated random intercept and slope
#'   \item \code{(1 | level1/level2)} - Nested random effects
#'   \item \code{(1 | group1) + (1 | group2)} - Crossed random effects
#' }
#'
#' @examples
#' \dontrun{
#' # Basic random intercept model
#' data(sleepstudy, package = "lme4")
#' fit1 <- gllamm(Reaction ~ Days + (1 | Subject),
#'                data = sleepstudy)
#' summary(fit1)
#'
#' # Random intercept and slope
#' fit2 <- gllamm(Reaction ~ Days + (Days | Subject),
#'                data = sleepstudy)
#' summary(fit2)
#'
#' # Extract components
#' fixef(fit2)        # Fixed effects
#' ranef(fit2)        # Random effects
#' VarCorr(fit2)      # Variance components
#' fitted(fit2)       # Fitted values
#' residuals(fit2)    # Residuals
#'
#' # Ordinal regression (proportional odds)
#' data$satisfaction <- ordered(sample(1:5, nrow(data), replace = TRUE))
#' fit3 <- gllamm(satisfaction ~ x + (1 | group),
#'                data = data,
#'                family = ordinal(link = "logit"))
#'
#' # Ordinal with adjacent category logit
#' fit4 <- gllamm(satisfaction ~ x + (1 | group),
#'                data = data,
#'                family = ordinal(link = "acl"))
#' }
#'
#' @references
#' Rabe-Hesketh, S., Skrondal, A., & Pickles, A. (2004). GLLAMM Manual.
#' U.C. Berkeley Division of Biostatistics Working Paper Series.
#'
#' Skrondal, A., & Rabe-Hesketh, S. (2004). Generalized Latent Variable
#' Modeling: Multilevel, Longitudinal, and Structural Equation Models.
#' Chapman & Hall/CRC.
#'
#' @export
gllamm <- function(formula,
                   data,
                   family = gaussian(),
                   weights = NULL,
                   random = NULL,
                   start = NULL,
                   control = list(),
                   ...) {

  # Capture call
  call <- match.call()

  # Validate inputs
  if (missing(formula)) {
    stop("Argument 'formula' is required")
  }

  # Matrix-response families dispatch before formula validation: the first
  # argument is the response matrix, `data` (optional) carries person-level
  # variables, `random` the person-level random-effects formula
  if (inherits(family, "irt_family")) {
    if (!is.matrix(formula) && !is.data.frame(formula)) {
      stop("For family = irt(), the first argument must be the persons x ",
           "items response matrix")
    }
    return(fit_irt(as.matrix(formula),
                   model = family$model,
                   person_data = if (missing(data)) NULL else data,
                   random = random,
                   weights = weights,
                   mc_items = family$mc_items,
                   start = start,
                   control = control))
  }

  if (inherits(family, "lca_family")) {
    if (!is.matrix(formula) && !is.data.frame(formula)) {
      stop("For family = lca(), the first argument must be the matrix of ",
           "binary manifest variables")
    }
    return(fit_lca(as.matrix(formula),
                   nclass = family$nclass,
                   weights = weights,
                   start = start,
                   control = control))
  }

  if (missing(data)) {
    stop("Argument 'data' is required")
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (inherits(family, "multinomial_family")) {
    if (!is.null(weights)) {
      warning("Weights are not yet supported for multinomial models; ignoring.")
    }
    return(fit_multinomial(formula = formula,
                           data = data,
                           reference = family$reference,
                           start = start,
                           control = control))
  }

  # Validate weights if provided
  if (!is.null(weights)) {
    if (length(weights) != nrow(data)) {
      stop("Length of weights (", length(weights), ") must match number of observations (", nrow(data), ")")
    }
    if (any(weights < 0, na.rm = TRUE)) {
      stop("All weights must be non-negative")
    }
    if (any(is.na(weights))) {
      stop("weights cannot contain missing values")
    }
  }

  # Validate formula
  validate_formula(formula, data)

  # Dispatch to specialized fitting functions based on family type
  # This provides a unified interface: gllamm(formula, data, family = ordinal(...))

  if (inherits(family, "ordinal_family")) {
    # Ordinal regression models
    return(fit_ordinal(formula = formula,
                      data = data,
                      link = family$link,
                      weights = weights,
                      start = start,
                      control = control))
  }

  if (inherits(family, "binomial_family")) {
    # Binomial regression models with custom link (logit, probit, cloglog)
    return(fit_binomial(formula = formula,
                       data = data,
                       link = family$link,
                       weights = weights,
                       start = start,
                       control = control))
  }

  # Parse formula
  parsed <- parse_formula(formula, data)

  if (length(parsed$random_terms) == 0) {
    stop("No random effects specified. Use glm() for fixed effects only models.")
  }

  # Create model matrices
  model_data <- make_model_matrices(parsed, data)

  # Fit model using TMB (use v2 interface if available)
  if (exists("fit_tmb_gllamm_v2")) {
    fit_result <- fit_tmb_gllamm_v2(
      model_data = model_data,
      family = family,
      random_terms = parsed$random_terms,
      start_params = start,
      control = control,
      weights = weights
    )
  } else {
    fit_result <- fit_tmb_gllamm(
      model_data = model_data,
      family = family,
      start_params = start,
      control = control,
      weights = weights
    )
  }

  # Calculate residuals
  resids <- model_data$y - fit_result$fitted.values

  # Calculate AIC and BIC
  n_obs <- model_data$n_obs
  n_params <- fit_result$n_params
  loglik <- fit_result$logLik

  aic <- -2 * loglik + 2 * n_params
  bic <- -2 * loglik + log(n_obs) * n_params

  # Construct gllamm object
  result <- list(
    coefficients = fit_result$coefficients,
    vcov = fit_result$vcov,
    random_effects = fit_result$random_effects,
    fitted.values = fit_result$fitted.values,
    residuals = resids,
    y = model_data$y,
    X = model_data$X,
    logLik = loglik,
    AIC = aic,
    BIC = bic,
    n_obs = n_obs,
    n_params = n_params,
    n_groups = model_data$n_groups,
    convergence = fit_result$convergence,
    call = call,
    formula = formula,
    family = family,
    data = data,
    random_terms = parsed$random_terms,
    tmb_obj = fit_result$tmb_obj,
    tmb_opt = fit_result$tmb_opt,
    tmb_sdr = fit_result$tmb_sdr
  )

  class(result) <- "gllamm"

  return(result)
}
