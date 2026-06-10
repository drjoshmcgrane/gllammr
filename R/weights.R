#' Parse observation-level and group-level survey weights
#'
#' GLLAMM-style level-specific weights: \code{weights} may be a numeric
#' vector (level-1 / observation weights, as before) or a list with elements
#' \code{level1} (length n_obs) and/or \code{level2} (one weight per group,
#' or per observation but constant within group). Level-2 weights scale each
#' group's full likelihood contribution including its random-effects prior
#' (pseudo-likelihood for two-stage sampling designs).
#'
#' @param weights NULL, numeric vector, or list(level1=, level2=)
#' @param n_obs Number of observations
#' @param groups 0-indexed group index per observation
#' @param n_groups Number of groups
#' @return list(level1 = numeric(n_obs), level2 = numeric(n_groups))
#' @keywords internal
parse_level_weights <- function(weights, n_obs, groups, n_groups) {
  w1 <- rep(1.0, n_obs)
  w2 <- rep(1.0, n_groups)

  if (is.null(weights)) {
    return(list(level1 = w1, level2 = w2))
  }

  if (is.list(weights)) {
    bad <- setdiff(names(weights), c("level1", "level2"))
    if (length(bad) > 0 || is.null(names(weights))) {
      stop("weights list elements must be named 'level1' and/or 'level2'")
    }
    if (!is.null(weights$level1)) {
      if (length(weights$level1) != n_obs) {
        stop("weights$level1 must have length ", n_obs)
      }
      w1 <- as.numeric(weights$level1)
    }
    if (!is.null(weights$level2)) {
      wl2 <- as.numeric(weights$level2)
      if (length(wl2) == n_groups) {
        w2 <- wl2
      } else if (length(wl2) == n_obs) {
        # Per-observation: must be constant within group
        agg <- tapply(wl2, groups, function(v) {
          if (length(unique(v)) > 1) NA_real_ else v[1]
        })
        if (anyNA(agg)) {
          stop("weights$level2 must be constant within each group")
        }
        w2 <- as.numeric(agg[order(as.integer(names(agg)))])
      } else {
        stop("weights$level2 must have length ", n_groups,
             " (per group) or ", n_obs, " (per observation)")
      }
    }
  } else {
    if (length(weights) != n_obs) {
      stop("weights must have length ", n_obs)
    }
    w1 <- as.numeric(weights)
  }

  if (any(c(w1, w2) < 0) || anyNA(c(w1, w2))) {
    stop("weights must be non-negative and non-missing")
  }
  list(level1 = w1, level2 = w2)
}
