#' Family Objects for GLLAMM Models
#'
#' Family constructors for generalized linear latent and mixed models
#'
#' @name families
NULL


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

  # Map link function to numeric code for TMB
  link_code <- switch(link,
    logit = 1L,
    probit = 2L,
    cloglog = 3L
  )

  structure(
    list(
      family = "binomial",
      link = link,
      link_code = link_code
    ),
    class = c("binomial_family", "family")
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
