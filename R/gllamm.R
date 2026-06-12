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
#' @param family A family object selecting the model class. Every model in
#'   the package is reachable here:
#'   \itemize{
#'     \item \code{gaussian()} (default), \code{binomial()},
#'       \code{poisson()}, \code{Gamma()} - GLMMs with any number of
#'       crossed/nested random-effects terms and random slopes
#'     \item \code{ordinal(link)} - cumulative, adjacent-category,
#'       continuation-ratio, and partial-proportional-odds models
#'     \item \code{multinomial(reference)} - baseline-category logit
#'     \item \code{irt(model)} / \code{eirt(item_data, ...)} - (explanatory)
#'       item response models; first argument is the response matrix
#'     \item \code{lca(nclass, ordering)} / \code{cdm(Q, model)} - latent
#'       class and cognitive diagnosis models; first argument is the
#'       response matrix
#'     \item \code{sem(measurement, structural)} - structural equation
#'       models; first argument is the data frame
#'     \item \code{mixed_response(...)} - joint mixed-type outcomes; first
#'       argument is the shared random-effects formula
#'     \item \code{survival_family(distribution)} - parametric frailty
#'       survival with \code{Surv(time, event)} on the left-hand side
#'     \item \code{ranking(case)} - rank-ordered (exploded) logit
#'   }
#'   The latent distribution is normal with Laplace integration by
#'   default; \code{integration = aghq(k)} requests adaptive quadrature
#'   and \code{integration = npml(k)} a nonparametric (mass-point)
#'   latent distribution.
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
#' @param weights Optional case weights: a numeric vector of observation
#'   (level-1) weights, or a list with elements \code{level1} and/or
#'   \code{level2} for survey designs with weights at both levels.
#'   Under the default Laplace fit, level-2 weights must be integer
#'   frequency weights; they are implemented by exact replication of
#'   whole groups, so results are identical to fitting the duplicated
#'   data. Non-integer level-2 weights require
#'   \code{integration = aghq(k)}, which weights each group's log
#'   marginal likelihood directly and supports arbitrary weights.
#' @param random For matrix-response families (\code{irt()}, \code{lca()}):
#'   an optional person-level random-effects formula such as
#'   \code{~ (1 | class)}.
#' @param integration Optional integration specification; \code{aghq(k)}
#'   selects adaptive Gauss-Hermite quadrature with \code{k} nodes for
#'   two-level random-intercept models (default is the Laplace
#'   approximation).
#' @param ... Additional arguments (reserved for future use).
#'
#' @section Missing data:
#' Formula-based models (GLMMs, ordinal, multinomial, survival, NPML,
#' mixed responses) listwise-delete rows with missing values in any model
#' variable, with a warning; supplied weights are aligned automatically.
#' Matrix-response latent variable models (\code{irt()}, \code{lca()},
#' \code{cdm()}, and the corresponding \code{fit_*} functions) use all
#' observed responses - item-level missingness is handled under the MAR
#' assumption by the marginal likelihood itself. \code{fit_sem} offers
#' full-information maximum likelihood via \code{missing = "fiml"}.
#' \code{fit_rank} treats missing ranks as deliberately unranked
#' alternatives (partial rankings), not as missing data.
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
#' gllammr uses Template Model Builder (TMB) for efficient computation via
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
                   integration = NULL,
                   start = NULL,
                   control = list(),
                   ...) {

  # Capture call
  call <- match.call()

  # Validate inputs
  if (missing(formula)) {
    stop("Argument 'formula' is required")
  }

  # SEM: the first argument is the data frame of indicators/covariates;
  # the model lives in the family object
  if (inherits(family, "sem_family")) {
    d_sem <- if (is.data.frame(formula)) formula
             else if (!missing(data)) data
             else stop("For family = sem(), pass the data frame as the ",
                       "first argument")
    return(fit_sem(measurement = family$measurement,
                   structural = family$structural,
                   data = d_sem,
                   missing = family$missing,
                   se = family$se,
                   start = start,
                   control = control))
  }

  # Mixed responses: the first argument is the shared random-effects
  # formula (e.g. ~ 1 | group); the outcome formulas live in the family
  if (inherits(family, "mixed_family")) {
    return(fit_mixed(formulas = family$formulas,
                     random = formula,
                     data = data,
                     start = start,
                     control = control))
  }

  # Parametric frailty survival: Surv(time, event) on the left-hand side
  if (inherits(family, "survival_family")) {
    return(fit_survival(formula = formula,
                        data = data,
                        distribution = family$distribution,
                        weights = weights,
                        start = start,
                        control = control))
  }

  # Rank-ordered (exploded) logit
  if (inherits(family, "rank_family")) {
    return(fit_rank(formula = formula,
                    case = family$case,
                    data = data,
                    random = random,
                    weights = weights,
                    start = start,
                    control = control))
  }

  # Matrix-response families dispatch before formula validation: the first
  # argument is the response matrix, `data` (optional) carries person-level
  # variables, `random` the person-level random-effects formula
  if (inherits(family, "eirt_family")) {
    if (!is.matrix(formula) && !is.data.frame(formula)) {
      stop("For family = eirt(), the first argument must be the persons x ",
           "items response matrix")
    }
    return(fit_eirt(as.matrix(formula),
                    item_data = family$item_data,
                    difficulty_formula = family$difficulty_formula,
                    discrimination_formula = family$discrimination_formula,
                    threshold_formula = family$threshold_formula,
                    person_data = if (missing(data)) NULL else data,
                    random = random,
                    weights = weights,
                    model = family$model,
                    item_residuals = family$item_residuals,
                    start = start,
                    control = control))
  }

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
                   ordering = family$ordering %||% "none",
                   weights = weights,
                   start = start,
                   control = control))
  }

  if (inherits(family, "cdm_family")) {
    if (!is.matrix(formula) && !is.data.frame(formula)) {
      stop("For family = cdm(), the first argument must be the persons x ",
           "items binary response matrix")
    }
    return(fit_cdm(as.matrix(formula),
                   Q = family$Q,
                   model = family$model,
                   hierarchy = family$hierarchy,
                   monotone = family$monotone,
                   weights = weights,
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

  # Validate weights if provided. Lists carry level-specific survey weights
  # (validated downstream against the grouping structure).
  if (!is.null(weights) && !is.list(weights)) {
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
    if (is.list(weights)) {
      stop("Level-specific weights are not supported for ordinal models; ",
           "supply a vector of observation weights")
    }
    # Ordinal regression models
    return(fit_ordinal(formula = formula,
                      data = data,
                      link = family$link,
                      weights = weights,
                      start = start,
                      control = control))
  }

  if (inherits(family, "binomial_family")) {
    # Binomial regression models with custom link (logit, probit, cloglog).
    # fit_binomial handles the single-random-term Laplace case; crossed/
    # multiple random-effects terms and non-default integration (aghq,
    # npml) route through the general engines, which support all three
    # links.
    n_re <- length(parse_formula(formula, data)$random_terms)
    # List (level-specific) weights need the general engine: fit_binomial
    # only understands observation-weight vectors
    if (n_re <= 1 && is.null(integration) && !is.list(weights)) {
      return(fit_binomial(formula = formula,
                          data = data,
                          link = family$link,
                          weights = weights,
                          start = start,
                          control = control))
    }
    family <- stats::binomial(link = family$link)
  }

  # Parse formula
  parsed <- parse_formula(formula, data)

  if (length(parsed$random_terms) == 0) {
    stop("No random effects specified. Use glm() for fixed effects only models.")
  }

  # Level-2 frequency weights under Laplace are implemented by exact
  # replication of whole groups (each copy gets its own random effect),
  # which reproduces duplicated-data fits exactly. Scaling each group's
  # likelihood-plus-prior contribution instead is only approximate under
  # the Laplace normalization and its objective is unbounded in principle
  # (every weighted group contributes -(w - 1) * log(sigma_u)). Adaptive
  # quadrature weights each group's log marginal likelihood outside the
  # integral, which is exact for arbitrary weights, so aghq fits pass
  # level-2 weights straight through.
  is_aghq <- inherits(integration, "gllamm_integration") &&
    identical(integration$method, "aghq")
  if (is.list(weights) && !is.null(weights$level2) && !is_aghq) {
    if (length(parsed$random_terms) != 1) {
      stop("Level-specific weights are currently supported for two-level ",
           "models (a single random-effects term)")
    }
    gvar <- parsed$random_terms[[1]]$grouping
    if (is.null(gvar) || !gvar %in% names(data)) {
      stop("Level-2 weights require a simple grouping variable in the data")
    }
    gfac <- factor(data[[gvar]])
    w2 <- as.numeric(weights$level2)
    if (length(w2) == nlevels(gfac)) {
      w2_obs <- w2[as.integer(gfac)]
    } else if (length(w2) == nrow(data)) {
      if (any(tapply(w2, gfac, function(v) length(unique(v))) > 1)) {
        stop("weights$level2 must be constant within each group")
      }
      w2_obs <- w2
    } else {
      stop("weights$level2 must have length ", nlevels(gfac),
           " (one per group) or ", nrow(data), " (one per observation)")
    }
    if (any(is.na(w2_obs)) || any(w2_obs < 0)) {
      stop("weights$level2 must be non-negative and complete")
    }
    if (any(abs(w2_obs - round(w2_obs)) > 1e-8)) {
      stop("Non-integer level-2 weights are not supported under the ",
           "Laplace approximation. Use integration = aghq(k), which ",
           "weights each group's log marginal likelihood exactly.")
    }
    w2_grp <- round(tapply(w2_obs, gfac, function(v) v[1]))
    idx <- unlist(lapply(seq_len(nlevels(gfac)), function(j) {
      rows <- which(as.integer(gfac) == j)
      rep(list(rows), w2_grp[j])
    }), recursive = FALSE)
    copy_id <- rep(seq_along(idx), lengths(idx))
    idx <- unlist(idx)
    data <- data[idx, , drop = FALSE]
    data[[gvar]] <- factor(paste0(as.character(gfac)[idx], ".rep", copy_id))
    w1 <- weights$level1
    weights <- if (!is.null(w1)) as.numeric(w1)[idx] else NULL
    parsed <- parse_formula(formula, data)
  }

  # Create model matrices
  model_data <- make_model_matrices(parsed, data)

  # Nonparametric maximum likelihood: discrete mass points instead of a
  # normal latent distribution
  if (inherits(integration, "gllamm_integration") &&
      identical(integration$method, "npml")) {
    return(fit_npml(formula = formula,
                    data = data,
                    k = integration$k,
                    family = family,
                    weights = if (is.list(weights)) {
                      stop("Level-specific weights are not supported ",
                           "under NPML")
                    } else weights,
                    start = start,
                    control = control))
  }

  # Adaptive quadrature integration (Laplace is the default)
  if (inherits(integration, "gllamm_integration") &&
      identical(integration$method, "aghq")) {
    fit_result <- fit_tmb_gllamm_aghq(
      model_data = model_data,
      family = family,
      random_terms = parsed$random_terms,
      k = integration$k,
      start_params = start,
      control = control,
      weights = weights
    )
  } else if (exists("fit_tmb_gllamm_v2")) {
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
