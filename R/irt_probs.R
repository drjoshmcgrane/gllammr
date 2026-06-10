#' Category response probabilities for polytomous IRT models
#'
#' Single source of truth for the category probability math, mirroring the
#' likelihood in the TMB template (gllamm_irt_poly.hpp / gllamm_eirt.hpp).
#' Used by plotting, DIF displays, and marginal predictions so the formulas
#' cannot drift apart.
#'
#' @param model One of "GRM", "PCM", "GPCM", "NRM"
#' @param theta Numeric vector of ability values
#' @param thresholds Numeric vector of item parameters: ordered thresholds
#'   (GRM), free step difficulties (PCM/GPCM), or category intercepts (NRM,
#'   reference category omitted)
#' @param discrimination Item discrimination (ignored for PCM)
#'
#' @return Matrix [length(theta) x K] of category probabilities; rows sum to 1
#' @keywords internal
irt_category_probs <- function(model, theta, thresholds, discrimination = 1) {
  n_theta <- length(theta)
  K <- length(thresholds) + 1
  a <- discrimination

  if (model == "GRM") {
    # Cumulative logits: P(Y >= k) = plogis(a * (theta - tau_k))
    cum <- vapply(seq_len(K - 1),
                  function(k) plogis(a * (theta - thresholds[k])),
                  numeric(n_theta))
    cum <- matrix(cum, n_theta, K - 1)
    probs <- matrix(0, n_theta, K)
    probs[, 1] <- 1 - cum[, 1]
    if (K > 2) {
      for (k in 2:(K - 1)) {
        probs[, k] <- cum[, k - 1] - cum[, k]
      }
    }
    probs[, K] <- cum[, K - 1]
    return(probs)
  }

  if (model %in% c("PCM", "GPCM")) {
    # Adjacent-categories logits: numerator_m = sum_{l<=m} a*(theta - delta_l)
    slope <- if (model == "PCM") 1 else a
    eta <- matrix(0, n_theta, K)
    for (m in 2:K) {
      eta[, m] <- eta[, m - 1] + slope * (theta - thresholds[m - 1])
    }
    eta <- eta - apply(eta, 1, max)   # guard against overflow
    expeta <- exp(eta)
    return(expeta / rowSums(expeta))
  }

  if (model == "NRM") {
    # Nominal: eta_1 = 0 (reference), eta_k = a*theta + c_k
    eta <- matrix(0, n_theta, K)
    for (k in 2:K) {
      eta[, k] <- a * theta + thresholds[k - 1]
    }
    eta <- eta - apply(eta, 1, max)
    expeta <- exp(eta)
    return(expeta / rowSums(expeta))
  }

  stop("Unknown polytomous model: ", model)
}
