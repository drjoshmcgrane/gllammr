#' Family Objects for GLLAMM Models
#'
#' Family constructors for generalized linear latent and mixed models
#'
#' @name families
NULL

# Register the S4 inheritance so objects from GLLAMMR's binomial() remain
# usable by S4-based consumers (e.g. lme4 slots typed "family")
methods::setOldClass(c("binomial_family", "family"))
methods::setOldClass(c("ordinal_family", "family"))
methods::setOldClass(c("irt_family", "family"))
methods::setOldClass(c("lca_family", "family"))
methods::setOldClass(c("multinomial_family", "family"))


#' Ordinal Family for Proportional and Non-Proportional Odds Models
#'
#' Create a family object for ordinal regression models with various link functions
#'
#' @param link Link function to use:
#'   \describe{
#'     \item{"logit"}{Proportional odds (cumulative logit) - default}
#'     \item{"probit"}{Cumulative probit}
#'     \item{"acl"}{Adjacent category logit}
#'     \item{"crl_forward"}{Forward continuation ratio logit}
#'     \item{"crl_backward"}{Backward continuation ratio logit}
#'     \item{"ppo"}{Partial proportional odds (non-proportional)}
#'   }
#'
#' @return A family object of class \code{ordinal_family} with components:
#'   \item{family}{Character: "ordinal"}
#'   \item{link}{Character: name of link function}
#'   \item{link_code}{Integer: numeric code for TMB (1-6)}
#'
#' @details
#' The ordinal family supports several models for ordered categorical responses:
#'
#' \strong{Proportional Odds (logit):}
#' \deqn{P(Y \le k | x) = \frac{1}{1 + \exp(-(\tau_k - x'\beta))}}
#'
#' \strong{Cumulative Probit:}
#' \deqn{P(Y \le k | x) = \Phi(\tau_k - x'\beta)}
#'
#' \strong{Adjacent Category Logit (ACL):}
#' Models the log-odds of adjacent categories:
#' \deqn{\log\frac{P(Y=k)}{P(Y=k-1)} = \alpha_k + x'\beta}
#'
#' \strong{Continuation Ratio Logit (CRL):}
#' Forward version models sequential decisions:
#' \deqn{\log\frac{P(Y=k | Y \ge k)}{P(Y>k | Y \ge k)} = \tau_k - x'\beta}
#'
#' Backward version reverses the conditioning.
#'
#' \strong{Partial Proportional Odds (PPO):}
#' Relaxes the proportional odds assumption by allowing different
#' covariate effects per threshold:
#' \deqn{P(Y \le k | x) = F(\tau_k - x'\beta_k)}
#'
#' @examples
#' \dontrun{
#' # Proportional odds model (default)
#' family1 <- ordinal()
#' family2 <- ordinal(link = "logit")
#'
#' # Adjacent category logit
#' family3 <- ordinal(link = "acl")
#'
#' # Partial proportional odds
#' family4 <- ordinal(link = "ppo")
#'
#' # Use with gllamm() - recommended interface
#' fit <- gllamm(rating ~ temp + (1 | judge),
#'               data = wine,
#'               family = ordinal(link = "logit"))
#'
#' # Or use fit_ordinal() directly
#' fit2 <- fit_ordinal(rating ~ temp + (1 | judge),
#'                     data = wine,
#'                     link = "acl")
#' }
#'
#' @export
ordinal <- function(link = c("logit", "probit", "acl", "crl_forward",
                             "crl_backward", "ppo")) {
  link <- match.arg(link)

  # Map link function to numeric code for TMB
  link_code <- switch(link,
    logit = 1L,
    probit = 2L,
    acl = 3L,
    crl_forward = 4L,
    crl_backward = 5L,
    ppo = 6L
  )

  structure(
    list(
      family = "ordinal",
      link = link,
      link_code = link_code
    ),
    class = c("ordinal_family", "family")
  )
}


#' Print method for ordinal family
#' @keywords internal
#' @export
print.ordinal_family <- function(x, ...) {
  cat("Family: ordinal\n")
  cat("Link function:", x$link, "\n")
  invisible(x)
}


#' Binomial Family for Binary and Binomial Outcomes
#'
#' Create a family object for binomial regression models with various link functions
#'
#' @param link Link function to use:
#'   \describe{
#'     \item{"logit"}{Logistic regression (default) - symmetric S-curve}
#'     \item{"probit"}{Probit regression - based on normal distribution}
#'     \item{"cloglog"}{Complementary log-log - asymmetric, suitable for rare events and survival data}
#'   }
#'
#' @return A family object of class \code{binomial_family} with components:
#'   \item{family}{Character: "binomial"}
#'   \item{link}{Character: name of link function}
#'   \item{link_code}{Integer: numeric code for TMB (1-3)}
#'
#' @details
#' The binomial family is used for binary (0/1) or binomial (count/total) responses.
#'
#' \strong{Logit Link (default):}
#' \deqn{P(Y=1|x) = \frac{1}{1 + \exp(-x'\beta)}}
#' This is the standard logistic regression. The logit link is symmetric and
#' appropriate when the probability of success and failure are equally likely
#' to change with covariates.
#'
#' \strong{Probit Link:}
#' \deqn{P(Y=1|x) = \Phi(x'\beta)}
#' where \eqn{\Phi} is the standard normal CDF. This link is also symmetric
#' and yields similar results to logit but with slightly different tail behavior.
#'
#' \strong{Complementary Log-Log (cloglog) Link:}
#' \deqn{P(Y=1|x) = 1 - \exp(-\exp(x'\beta))}
#' This link is \emph{asymmetric} and is particularly useful for:
#' \itemize{
#'   \item Rare events (when P(Y=1) is small)
#'   \item Survival analysis with discrete time intervals
#'   \item Gompertz or extreme value distributions
#'   \item When the hazard is proportional (as in survival models)
#' }
#'
#' The cloglog link arises naturally when modeling grouped survival data
#' or when events follow a Poisson process over time intervals.
#'
#' @examples
#' \dontrun{
#' # Logistic regression (default)
#' family1 <- binomial()
#' family2 <- binomial(link = "logit")
#'
#' # Probit regression
#' family3 <- binomial(link = "probit")
#'
#' # Complementary log-log for rare events
#' family4 <- binomial(link = "cloglog")
#'
#' # Use with gllamm() - recommended interface
#' fit <- gllamm(outcome ~ age + treatment + (1 | clinic),
#'               data = mydata,
#'               family = binomial(link = "logit"))
#'
#' # Rare event with cloglog
#' fit_rare <- gllamm(rare_disease ~ exposure + (1 | region),
#'                    data = epi_data,
#'                    family = binomial(link = "cloglog"))
#' }
#'
#' @seealso \code{\link{gllamm}}, \code{\link{ordinal}}
#'
#' @export
binomial <- function(link = c("logit", "probit", "cloglog")) {
  link <- match.arg(link)

  # Build on the full stats family object so the result remains usable by
  # other packages (lme4, glm, ...) when GLLAMMR masks stats::binomial
  f <- stats::binomial(link = link)
  f$link_code <- switch(link,
    logit = 1L,
    probit = 2L,
    cloglog = 3L
  )
  class(f) <- c("binomial_family", class(f))
  f
}


#' IRT Family for Item Response Theory Models
#'
#' Create a family object for fitting IRT models through the unified
#' \code{gllamm()} interface. The response is a persons x items matrix
#' passed as the first argument of \code{gllamm()}.
#'
#' @param model IRT model type: "Rasch", "2PL", "3PL" (dichotomous) or
#'   "GRM", "PCM", "GPCM", "NRM" (polytomous)
#' @param mc_items For 3PL only: which items have guessing parameters
#'   (NULL = all; logical or integer index vector)
#'
#' @return A family object of class \code{irt_family}
#'
#' @examples
#' \dontrun{
#' fit <- gllamm(response_matrix, family = irt("2PL"))
#' # Multi-level IRT: persons nested in classes
#' fit_ml <- gllamm(response_matrix, data = person_data,
#'                  family = irt("Rasch"), random = ~ (1 | class))
#' }
#'
#' @export
irt <- function(model = c("Rasch", "2PL", "3PL", "GRM", "PCM", "GPCM", "NRM"),
                mc_items = NULL) {
  model <- match.arg(model)
  structure(
    list(family = "irt", model = model, mc_items = mc_items),
    class = c("irt_family", "family")
  )
}


#' Latent Class Family for Finite Mixture Models
#'
#' Create a family object for fitting latent class models through the unified
#' \code{gllamm()} interface. The response is a matrix of binary manifest
#' variables passed as the first argument of \code{gllamm()}.
#'
#' @param nclass Number of latent classes (default 2)
#' @param ordering Class order restriction passed to \code{\link{fit_lca}}:
#'   "none" (default), "increasing", or a list/matrix of class pairs
#'   defining a partial order
#'
#' @return A family object of class \code{lca_family}
#'
#' @examples
#' \dontrun{
#' fit <- gllamm(indicator_matrix, family = lca(nclass = 3))
#' fit_ord <- gllamm(indicator_matrix,
#'                   family = lca(nclass = 3, ordering = "increasing"))
#' }
#'
#' @export
lca <- function(nclass = 2, ordering = "none") {
  if (!is.numeric(nclass) || length(nclass) != 1 || nclass < 2) {
    stop("nclass must be a single integer >= 2")
  }
  structure(
    list(family = "lca", nclass = as.integer(nclass), ordering = ordering),
    class = c("lca_family", "family")
  )
}


#' Cognitive Diagnosis Family for Q-Matrix Models
#'
#' Create a family object for fitting cognitive diagnosis models through
#' the unified \code{gllamm()} interface. The response is a persons x items
#' binary matrix passed as the first argument of \code{gllamm()}; the
#' Q-matrix and model options are carried by the family object. See
#' \code{\link{fit_cdm}} for details of the models and arguments.
#'
#' @param Q Binary Q-matrix (items x attributes)
#' @param model "gdina" (default), "dina", or "dino"
#' @param hierarchy Optional attribute hierarchy (list of prerequisite
#'   pairs)
#' @param monotone Enforce monotonicity in the attributes (default TRUE)
#'
#' @return An object of class \code{cdm_family}
#'
#' @examples
#' \dontrun{
#' fit <- gllamm(Y, family = cdm(Q, model = "dina"))
#' }
#'
#' @export
cdm <- function(Q, model = c("gdina", "dina", "dino"),
                hierarchy = NULL, monotone = TRUE) {
  model <- match.arg(model)
  structure(
    list(family = "cdm", Q = Q, model = model,
         hierarchy = hierarchy, monotone = monotone),
    class = c("cdm_family", "family")
  )
}


#' Explanatory IRT Family
#'
#' Create a family object for explanatory item response models through the
#' unified \code{gllamm()} interface. The response is the persons x items
#' matrix passed as the first argument of \code{gllamm()}; \code{data}
#' (optional) carries person-level variables and \code{random} person-level
#' random effects. See \code{\link{fit_eirt}}.
#'
#' @param item_data Data frame of item covariates (one row per item)
#' @param difficulty_formula Item-covariate formula for difficulty
#' @param discrimination_formula Item-covariate formula for (log)
#'   discrimination (2PL/GRM/GPCM)
#' @param threshold_formula Optional threshold regression (LPCM framework)
#' @param model "Rasch", "2PL", "GRM", "PCM", or "GPCM"
#' @param item_residuals Random item residuals around the regression
#'   (LLTM-plus-error; default TRUE)
#'
#' @return An object of class \code{eirt_family}
#'
#' @examples
#' \dontrun{
#' fit <- gllamm(resp, family = eirt(item_data,
#'                                   difficulty_formula = ~ btype + mode))
#' }
#'
#' @export
eirt <- function(item_data, difficulty_formula = ~ 1,
                 discrimination_formula = ~ 1, threshold_formula = NULL,
                 model = c("Rasch", "2PL", "GRM", "PCM", "GPCM"),
                 item_residuals = TRUE) {
  model <- match.arg(model)
  structure(
    list(family = "eirt", item_data = item_data,
         difficulty_formula = difficulty_formula,
         discrimination_formula = discrimination_formula,
         threshold_formula = threshold_formula,
         model = model, item_residuals = item_residuals),
    class = c("eirt_family", "family")
  )
}


#' SEM Family for Structural Equation Models
#'
#' Create a family object for structural equation models through the
#' unified \code{gllamm()} interface. The data frame is passed as the
#' first argument of \code{gllamm()}; the measurement and structural
#' models live in the family object. See \code{\link{fit_sem}}.
#'
#' @param measurement Named list of one-sided indicator formulas
#' @param structural Optional list of structural formulas (latent and/or
#'   observed predictors)
#' @param missing "listwise" (default) or "fiml"
#' @param se Compute standard errors (default TRUE)
#'
#' @return An object of class \code{sem_family}
#'
#' @examples
#' \dontrun{
#' fit <- gllamm(d, family = sem(
#'   measurement = list(f1 = ~ x1 + x2 + x3, f2 = ~ y1 + y2 + y3),
#'   structural = list(f2 ~ f1)))
#' }
#'
#' @export
sem <- function(measurement, structural = NULL,
                missing = c("listwise", "fiml"), se = TRUE) {
  structure(
    list(family = "sem", measurement = measurement,
         structural = structural, missing = match.arg(missing), se = se),
    class = c("sem_family", "family")
  )
}


#' Mixed-Response Family for Joint Outcome Models
#'
#' Create a family object for joint models of mixed-type outcomes sharing
#' a random effect, through the unified \code{gllamm()} interface. The
#' first argument of \code{gllamm()} is the shared random-effects formula
#' (e.g. \code{~ 1 | group}). See \code{\link{fit_mixed}}.
#'
#' @param ... Named outcome formulas: \code{gaussian = y1 ~ x},
#'   \code{binomial = y2 ~ x}, \code{poisson = y3 ~ x} (any subset)
#'
#' @return An object of class \code{mixed_family}
#'
#' @examples
#' \dontrun{
#' fit <- gllamm(~ 1 | clinic, data = d,
#'               family = mixed_response(gaussian = severity ~ age,
#'                                       binomial = dropout ~ age))
#' }
#'
#' @export
mixed_response <- function(...) {
  formulas <- list(...)
  if (length(formulas) == 0 || is.null(names(formulas)) ||
      any(names(formulas) == "")) {
    stop("mixed_response() needs named outcome formulas, e.g. ",
         "mixed_response(gaussian = y1 ~ x, binomial = y2 ~ x)")
  }
  structure(
    list(family = "mixed_response", formulas = formulas),
    class = c("mixed_family", "family")
  )
}


#' Rank-Ordered Logit Family
#'
#' Create a family object for rank-ordered (exploded) logit models through
#' the unified \code{gllamm()} interface. See \code{\link{fit_rank}}.
#'
#' @param case Case (chooser) identifier: a one-sided formula
#'   (\code{~ chooser}) or a variable name
#'
#' @return An object of class \code{rank_family}
#'
#' @examples
#' \dontrun{
#' fit <- gllamm(rank ~ price + quality, data = d,
#'               family = ranking(case = ~ chooser),
#'               random = ~ (1 | chooser))
#' }
#'
#' @export
ranking <- function(case) {
  if (is.character(case)) {
    case <- stats::as.formula(paste("~", case))
  }
  if (!inherits(case, "formula")) {
    stop("case must be a one-sided formula or a variable name")
  }
  structure(
    list(family = "ranking", case = case),
    class = c("rank_family", "family")
  )
}


#' Parametric Frailty Survival Family
#'
#' Create a family object for parametric survival models with shared
#' (log-normal) frailties through the unified \code{gllamm()} interface;
#' the formula uses \code{Surv(time, event)} on the left-hand side. See
#' \code{\link{fit_survival}}.
#'
#' @param distribution "exponential" (default) or "weibull"
#'
#' @return An object of class \code{survival_family}
#'
#' @examples
#' \dontrun{
#' fit <- gllamm(Surv(time, status) ~ x + (1 | clinic), data = d,
#'               family = survival_family("weibull"))
#' }
#'
#' @export
survival_family <- function(distribution = c("exponential", "weibull")) {
  structure(
    list(family = "survival", distribution = match.arg(distribution)),
    class = c("survival_family", "family")
  )
}


#' Multinomial Family for Unordered Categorical Outcomes
#'
#' Create a family object for baseline-category multinomial logit models
#' through the unified \code{gllamm()} interface.
#'
#' @param reference Reference category (default: first level)
#'
#' @return A family object of class \code{multinomial_family}
#'
#' @examples
#' \dontrun{
#' fit <- gllamm(choice ~ x + (1 | region), data = d, family = multinomial())
#' }
#'
#' @export
multinomial <- function(reference = NULL) {
  structure(
    list(family = "multinomial", reference = reference),
    class = c("multinomial_family", "family")
  )
}


#' Print method for binomial family
#' @keywords internal
#' @export
print.binomial_family <- function(x, ...) {
  cat("Family: binomial\n")
  cat("Link function:", x$link, "\n")
  invisible(x)
}
