#' Fit Rank-Ordered Logit Models with Taste Heterogeneity
#'
#' Fits an exploded (rank-ordered) logit model for ranking data, with an
#' optional group-level random coefficient on one alternative attribute
#' (a "taste shifter"). A constant-within-case random intercept cancels in
#' ranking likelihoods, so heterogeneity must attach to an attribute that
#' varies across alternatives.
#'
#' @param formula Model formula \code{rank ~ x1 + x2} where the response is
#'   the rank of each alternative within its case (1 = most preferred).
#'   Unranked alternatives may be coded \code{NA}; they remain in every
#'   choice set without contributing a stage.
#' @param case One-sided formula naming the case identifier, e.g. \code{~ id}
#' @param data Data frame in long format (one row per alternative per case)
#' @param random Optional one-sided formula \code{~ (0 + attribute | group)}
#'   giving the attribute carrying the random coefficient and the grouping
#'   variable. NULL fits a fixed-effects rank-ordered logit (the random SD
#'   is fixed near zero).
#' @param weights Optional one weight per case (matched via the case id)
#' @param start Optional starting values
#' @param control Optimization control list
#'
#' @return An object of class \code{gllamm_rank}
#'
#' @examples
#' \dontrun{
#' fit <- fit_rank(rank ~ price + quality, case = ~ subject,
#'                 random = ~ (0 + price | region), data = d)
#' }
#'
#' @export
fit_rank <- function(formula, case, data, random = NULL,
                     weights = NULL, start = NULL, control = list()) {

  rank_var <- as.character(formula[[2]])
  case_var <- all.vars(case)
  if (length(case_var) != 1 || !case_var %in% names(data)) {
    stop("case must name one variable present in data")
  }

  # ---- Random-coefficient specification ----
  if (!is.null(random)) {
    re_term <- attr(terms(random), "term.labels")
    if (length(re_term) != 1 || !grepl("\\|", re_term)) {
      stop("random must be of the form ~ (0 + attribute | group)")
    }
    parts <- strsplit(gsub("[()]", "", re_term), "\\|")[[1]]
    attr_expr <- gsub("^\\s*0\\s*\\+\\s*", "", trimws(parts[1]))
    group_var <- trimws(parts[2])
    if (!group_var %in% names(data)) {
      stop("Grouping variable '", group_var, "' not in data")
    }
    Zu <- eval(parse(text = attr_expr), envir = data)
    if (length(unique(Zu)) <= 1) {
      stop("The random-coefficient attribute must vary across alternatives ",
           "(a constant cancels in ranking likelihoods)")
    }
    group_factor <- factor(data[[group_var]])
  } else {
    Zu <- rep(1, nrow(data))
    group_factor <- factor(rep(1, nrow(data)))
  }

  # ---- Sort rows by (case, rank), unranked (NA) last ----
  rank_vals <- data[[rank_var]]
  ord <- order(data[[case_var]], ifelse(is.na(rank_vals), Inf, rank_vals))
  d_sorted <- data[ord, , drop = FALSE]
  Zu <- as.numeric(Zu)[ord]
  rank_sorted <- d_sorted[[rank_var]]

  case_f <- factor(d_sorted[[case_var]], levels = unique(d_sorted[[case_var]]))
  n_cases <- nlevels(case_f)
  case_rows <- split(seq_len(nrow(d_sorted)), case_f)
  case_start <- vapply(case_rows, function(r) r[1] - 1L, integer(1))
  case_n_alts <- vapply(case_rows, length, integer(1))
  case_n_ranked <- vapply(case_rows, function(r) {
    sum(!is.na(rank_sorted[r]))
  }, integer(1))

  if (any(case_n_ranked < 2 & case_n_alts == case_n_ranked)) {
    stop("Each case must rank at least two alternatives")
  }
  # Validate rank coding 1..n_ranked within case
  for (cr in case_rows) {
    rk <- rank_sorted[cr]
    rk <- rk[!is.na(rk)]
    if (!identical(sort(as.integer(rk)), seq_along(rk))) {
      stop("Ranks must be 1, 2, ... within each case (no ties)")
    }
  }

  case_group <- vapply(case_rows, function(r) {
    as.integer(group_factor[ord][r[1]]) - 1L
  }, integer(1))
  n_groups <- nlevels(group_factor)

  if (is.null(weights)) {
    case_weights <- rep(1.0, n_cases)
  } else {
    if (length(weights) != n_cases) {
      stop("weights must have one entry per case (", n_cases, ")")
    }
    case_weights <- as.numeric(weights)
  }

  # Fixed-effects design: no intercept (cancels within choice sets)
  X <- model.matrix(update(formula, NULL ~ . - 1), data = d_sorted)

  tmb_data <- list(
    X = X,
    Zu = Zu,
    case_start = as.integer(case_start),
    case_n_alts = as.integer(case_n_alts),
    case_n_ranked = as.integer(case_n_ranked),
    case_group = as.integer(case_group),
    n_cases = as.integer(n_cases),
    n_groups = as.integer(n_groups),
    case_weights = case_weights,
    model_name = "rank"
  )

  if (is.null(start)) {
    tmb_params <- list(
      beta = rep(0, ncol(X)),
      u = rep(0, n_groups),
      log_sigma_u = log(0.5)
    )
  } else {
    tmb_params <- start
  }

  tmb_map <- list()
  if (is.null(random)) {
    # Fixed-effects model: freeze the (single, dummy) random effect at zero
    tmb_map$u <- factor(rep(NA, n_groups))
    tmb_map$log_sigma_u <- factor(NA)
    tmb_params$u <- rep(0, n_groups)
  }

  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    random = if (is.null(random)) NULL else "u",
    map = tmb_map,
    DLL = "GLLAMMR",
    silent = TRUE
  )

  control_defaults <- list(eval.max = 2000, iter.max = 1000, trace = 0)
  control <- modifyList(control_defaults, control)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)

  sdr <- try(TMB::sdreport(obj), silent = TRUE)
  par_full <- obj$env$last.par.best

  beta_hat <- par_full[names(par_full) == "beta"]
  names(beta_hat) <- colnames(X)
  sigma_u_hat <- if (is.null(random)) {
    NA_real_
  } else {
    exp(unname(par_full[names(par_full) == "log_sigma_u"]))
  }

  n_params <- length(beta_hat) + as.integer(!is.null(random))

  result <- list(
    coefficients = list(fixed = beta_hat, random_sd = sigma_u_hat),
    logLik = -opt$objective,
    AIC = 2 * opt$objective + 2 * n_params,
    BIC = 2 * opt$objective + log(n_cases) * n_params,
    convergence = list(converged = (opt$convergence == 0),
                       message = opt$message),
    n_cases = n_cases,
    n_groups = if (is.null(random)) NA_integer_ else n_groups,
    formula = formula,
    random = random,
    tmb_obj = obj,
    tmb_opt = opt,
    tmb_sdr = sdr
  )
  class(result) <- c("gllamm_rank", "gllamm")
  result
}


#' @export
print.gllamm_rank <- function(x, ...) {
  cat("Rank-Ordered (Exploded) Logit\n\n")
  cat("Cases:", x$n_cases, "\n")
  if (!is.na(x$n_groups)) {
    cat("Groups:", x$n_groups,
        " Random-coefficient SD:", round(x$coefficients$random_sd, 4), "\n")
  }
  cat("\nPreference coefficients:\n")
  print(round(x$coefficients$fixed, 4))
  cat("\nLog-likelihood:", round(x$logLik, 2), "\n")
  invisible(x)
}
