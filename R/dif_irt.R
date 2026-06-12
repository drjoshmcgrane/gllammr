#' Confirmatory model-based DIF (IRT likelihood-ratio tests)
#'
#' Tests differential item functioning inside the joint marginal-ML IRT
#' model (IRT-LR DIF; Thissen, Steinberg & Wainer 1993), the De Boeck &
#' Wilson item-by-covariate-interaction formulation. Group differences in
#' ability (impact) are modeled by a latent regression
#' \eqn{\theta_p \sim N(z_p'\gamma, \sigma^2)}, so DIF parameters measure
#' item bias beyond true ability differences - the model-based companion
#' to the observed-criterion screening tests in \code{\link{dif_test}}.
#'
#' @param response_matrix Binary item response matrix (persons x items;
#'   NA allowed)
#' @param dif DIF specification: a grouping vector or a one-sided formula
#'   over \code{person_data} (multiple variables and interactions
#'   supported, as in \code{\link{dif_test}})
#' @param person_data Data frame with the DIF variables (required when
#'   \code{dif} is a formula)
#' @param model "Rasch" (default) or "2PL"
#' @param items Item indices to test (default: all non-anchor items)
#' @param anchors Optional indices of items constrained DIF-free
#' @param type "uniform" (default; covariate shifts of the item logit,
#'   \eqn{z'\delta_i}) or "both" (additionally nonuniform DIF: covariate
#'   scaling of the discrimination, \eqn{a_i e^{z'\kappa_i}}; 2PL only)
#' @param method "lr" (default; per-item likelihood-ratio tests, refitting
#'   the joint model) or "wald" (one joint fit with all studied items
#'   free - requires explicit \code{anchors} for identification - and
#'   per-item block Wald tests)
#' @param purify Purified IRT-LR (default FALSE): after each round, items
#'   flagged so far keep free DIF parameters in both compared models, so
#'   their misfit cannot contaminate the impact estimate; iterate until
#'   the flag set stabilizes
#' @param alpha Significance level
#' @param p_adjust Multiple-testing correction (\code{p.adjust} method)
#' @param max_iter Maximum purification rounds
#' @param control Optimization control list
#'
#' @return An object of class \code{dif_irt}: \code{dif_results}
#'   (per item: LR or Wald chi-square, df, p, adjusted p, flag, and the
#'   estimated uniform DIF effects \eqn{\delta} on the logit metric with
#'   standard errors), \code{impact} (latent regression coefficients
#'   \eqn{\gamma} with SEs - the ability difference attributable to the
#'   covariates), \code{flagged_items}, and purification details.
#'
#' @details
#' For a studied item the compared models differ only in that item's DIF
#' parameters; all item difficulties (and discriminations) are estimated
#' jointly with the latent regression, so the test is a genuine
#' marginal-likelihood ratio with q (uniform) or 2q ("both") degrees of
#' freedom, q the number of DIF design columns. For the Rasch model with
#' uniform DIF this is likelihood-equivalent to the long-format GLMM
#' \code{y ~ 0 + item + z + item_j:z + (1 | person)} under the same
#' Laplace approximation (verified in the test suite).
#'
#' @examples
#' \dontrun{
#' # Screen, then confirm
#' screen <- dif_test(resp, dif = ~ gender * language, person_data = pd)
#' confirm <- dif_irt(resp, dif = ~ gender * language, person_data = pd,
#'                    items = screen$flagged_items)
#' summary(confirm)
#' }
#'
#' @export
dif_irt <- function(response_matrix, dif, person_data = NULL,
                    model = c("Rasch", "2PL"),
                    items = NULL, anchors = NULL,
                    type = c("uniform", "both"),
                    method = c("lr", "wald"),
                    purify = FALSE,
                    alpha = 0.05, p_adjust = "none",
                    max_iter = 10, control = list()) {
  model <- match.arg(model)
  type <- match.arg(type)
  method <- match.arg(method)
  response_matrix <- as.matrix(response_matrix)
  n_persons <- nrow(response_matrix)
  n_items <- ncol(response_matrix)
  item_names <- colnames(response_matrix) %||% paste0("Item", 1:n_items)

  vals <- response_matrix[!is.na(response_matrix)]
  if (!all(vals %in% c(0, 1))) {
    stop("dif_irt requires binary (0/1) responses")
  }
  if (type == "both" && model == "Rasch") {
    stop("Nonuniform DIF requires model = \"2PL\"")
  }

  # ---- DIF design (no intercept; reference profile has latent mean 0) ----
  if (inherits(dif, "formula")) {
    if (is.null(person_data)) {
      stop("person_data is required when dif is a formula")
    }
    if (nrow(person_data) != n_persons) {
      stop("person_data must have one row per person")
    }
    mm <- model.matrix(dif, data = person_data)
    dif_formula <- dif
  } else {
    if (length(dif) != n_persons) {
      stop("dif length must match the number of persons")
    }
    person_data <- data.frame(group = factor(dif))
    mm <- model.matrix(~ group, data = person_data)
    dif_formula <- ~ group
  }
  Zp <- mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
  q <- ncol(Zp)
  if (q == 0) stop("The DIF specification has no terms to test")
  dif_terms <- colnames(Zp)

  if (!is.null(anchors) && (any(anchors < 1) || any(anchors > n_items))) {
    stop("anchors must be valid item indices")
  }
  if (is.null(items)) {
    items <- setdiff(seq_len(n_items), anchors)
  }
  if (any(items < 1 | items > n_items)) {
    stop("items must be between 1 and ", n_items)
  }
  items <- setdiff(items, anchors)

  nonuniform <- as.integer(type == "both")
  df_test <- q * (1L + nonuniform)

  # ---- Long format ----
  y_long <- as.vector(response_matrix)
  pid <- rep(seq_len(n_persons), times = n_items) - 1L
  iid <- rep(seq_len(n_items), each = n_persons) - 1L
  keep <- !is.na(y_long)
  y_long <- y_long[keep]; pid <- pid[keep]; iid <- iid[keep]

  # ---- One joint MML fit with a given studied set ----
  fit_one <- function(studied, want_sdr = FALSE) {
    studied <- sort(unique(studied))
    n_dif <- length(studied)
    dif_item <- rep(-1L, n_items)
    if (n_dif) dif_item[studied] <- seq_len(n_dif) - 1L

    tmb_data <- list(
      y = as.numeric(y_long), person_id = pid, item_id = iid,
      n_persons = as.integer(n_persons), n_items = as.integer(n_items),
      n_obs = length(y_long), Zp = Zp, dif_item = dif_item,
      model_type = if (model == "Rasch") 1L else 2L,
      nonuniform = nonuniform, model_name = "irt_dif")

    nd <- max(n_dif, 1L)
    tmb_params <- list(
      theta = rep(0, n_persons),
      difficulty = rep(0, n_items),
      discrimination = rep(1, n_items),
      log_sigma_theta = 0,
      gamma_impact = rep(0, q),
      delta = matrix(0, nd, q),
      kappa = matrix(0, nd, q))

    tmb_map <- list()
    if (model == "Rasch") {
      tmb_map$discrimination <- factor(rep(NA, n_items))
    } else {
      tmb_map$log_sigma_theta <- factor(NA)   # scale on discriminations
    }
    if (n_dif == 0) tmb_map$delta <- factor(rep(NA, nd * q))
    if (n_dif == 0 || nonuniform == 0) {
      tmb_map$kappa <- factor(rep(NA, nd * q))
    }

    obj <- TMB::MakeADFun(data = tmb_data, parameters = tmb_params,
                          random = "theta", map = tmb_map,
                          DLL = "gllammr", silent = TRUE)
    ctl <- modifyList(list(eval.max = 2000, iter.max = 1000, trace = 0),
                      control)
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = ctl)
    out <- list(logLik = -opt$objective,
                converged = (opt$convergence == 0), obj = obj)
    pf <- obj$env$last.par.best
    out$gamma <- unname(pf[names(pf) == "gamma_impact"])
    if (n_dif) {
      out$delta <- matrix(pf[names(pf) == "delta"], n_dif, q)
      if (nonuniform == 1) {
        out$kappa <- matrix(pf[names(pf) == "kappa"], n_dif, q)
      }
    }
    if (want_sdr) {
      out$sdr <- try(TMB::sdreport(obj), silent = TRUE)
    }
    out
  }

  history <- list()
  converged_pur <- TRUE

  if (method == "wald") {
    if (is.null(anchors) || length(anchors) < 1) {
      stop("method = \"wald\" needs explicit anchors (the joint model ",
           "with every item studied is not identified against impact)")
    }
    items <- sort(items)
    joint <- fit_one(items, want_sdr = TRUE)
    if (inherits(joint$sdr, "try-error")) {
      stop("sdreport failed; use method = \"lr\"")
    }
    V <- joint$sdr$cov.fixed
    fixed_par <- joint$obj$env$last.par.best[-joint$obj$env$random]
    didx <- which(names(joint$obj$par) == "delta")
    kidx <- which(names(joint$obj$par) == "kappa")
    n_dif <- length(items)
    res_rows <- lapply(seq_along(items), function(r) {
      sel <- didx[r + n_dif * (seq_len(q) - 1)]      # column-major rows
      if (nonuniform == 1) {
        sel <- c(sel, kidx[r + n_dif * (seq_len(q) - 1)])
      }
      est <- fixed_par[sel]
      Vb <- V[sel, sel, drop = FALSE]
      chi <- tryCatch(as.numeric(t(est) %*% solve(Vb, est)),
                      error = function(e) NA_real_)
      data.frame(chisq = chi, df = df_test,
                 p_value = stats::pchisq(chi, df_test, lower.tail = FALSE))
    })
    res <- do.call(rbind, res_rows)
    res[paste0("delta_", dif_terms)] <- joint$delta
    res[paste0("se_delta_", dif_terms)] <-
      matrix(sqrt(pmax(diag(V)[didx], 0)), n_dif, q)
    gamma_fit <- joint
    final_flag_basis <- res
  } else {
    # ---- (Purified) IRT-LR loop ----
    flagged <- integer(0)
    iter <- 1
    repeat {
      F_set <- intersect(flagged, items)
      m_base <- fit_one(F_set)               # DIF free for flagged set
      base_cache <- list()
      rows <- lapply(items, function(j) {
        if (j %in% F_set) {
          m0 <- fit_one(setdiff(F_set, j))
          m1 <- m_base
        } else {
          m0 <- m_base
          m1 <- fit_one(union(F_set, j))
        }
        chi <- max(2 * (m1$logLik - m0$logLik), 0)
        dj <- match(j, sort(unique(union(F_set, j))))
        data.frame(chisq = chi, df = df_test,
                   p_value = stats::pchisq(chi, df_test,
                                           lower.tail = FALSE),
                   t(setNames(m1$delta[dj, ], paste0("delta_",
                                                     dif_terms))))
      })
      res <- do.call(rbind, rows)
      p_adj <- stats::p.adjust(res$p_value, method = p_adjust)
      new_flagged <- items[!is.na(p_adj) & p_adj < alpha]
      history[[iter]] <- new_flagged
      if (!purify || setequal(new_flagged, flagged)) {
        flagged <- new_flagged
        break
      }
      flagged <- new_flagged
      if (iter >= max_iter) {
        warning("Purification did not stabilize in ", max_iter,
                " iterations")
        converged_pur <- FALSE
        break
      }
      iter <- iter + 1
    }
    # Final joint fit with the flagged items studied, for impact + SEs
    gamma_fit <- fit_one(intersect(flagged, items), want_sdr = TRUE)
    delta_hat <- NULL
    delta_se <- NULL
    final_flag_basis <- res
  }

  res <- final_flag_basis
  res$item <- items
  res$name <- item_names[items]
  res$p_adj <- stats::p.adjust(res$p_value, method = p_adjust)
  res$flagged <- !is.na(res$p_adj) & res$p_adj < alpha
  flagged_items <- res$item[res$flagged]

  # Impact coefficients with SEs from the final fit
  gamma_se <- rep(NA_real_, q)
  if (!is.null(gamma_fit$sdr) && !inherits(gamma_fit$sdr, "try-error")) {
    ss <- summary(gamma_fit$sdr, "report")
    gi <- rownames(ss) == "gamma_impact"
    if (any(gi)) gamma_se <- ss[gi, "Std. Error"]
  }
  impact <- data.frame(term = dif_terms, gamma = gamma_fit$gamma,
                       se = gamma_se,
                       z = gamma_fit$gamma / gamma_se,
                       stringsAsFactors = FALSE)

  front <- c("item", "name", "chisq", "df", "p_value", "p_adj", "flagged")
  res <- res[, c(front, setdiff(names(res), front))]

  result <- list(
    dif_results = res,
    flagged_items = flagged_items,
    impact = impact,
    anchors = anchors,
    dif_formula = dif_formula,
    dif_terms = dif_terms,
    model = model,
    type = type,
    method = method,
    purification = list(purify = purify,
                        iterations = length(history),
                        history = history,
                        converged = converged_pur),
    alpha = alpha,
    p_adjust = p_adjust,
    logLik_final = gamma_fit$logLik
  )
  class(result) <- "dif_irt"
  result
}


#' @export
print.dif_irt <- function(x, ...) {
  cat("Confirmatory IRT DIF analysis (", x$model, ", ",
      if (x$type == "both") "uniform + nonuniform" else "uniform",
      " DIF, ", toupper(x$method), " tests)\n", sep = "")
  cat("DIF specification:", deparse(x$dif_formula), "\n")
  if (!is.null(x$anchors)) {
    cat("Anchors:", paste(x$anchors, collapse = ", "), "\n")
  }
  if (isTRUE(x$purification$purify)) {
    cat("Purification:", x$purification$iterations, "round(s),",
        if (x$purification$converged) "converged" else "NOT converged",
        "\n")
  }
  cat("\nImpact (latent regression of ability on the DIF variables):\n")
  imp <- x$impact
  imp[, -1] <- round(imp[, -1], 4)
  print(imp, row.names = FALSE)
  cat("\nItems flagged:", if (length(x$flagged_items))
    paste(x$flagged_items, collapse = ", ") else "none",
    "(alpha =", x$alpha, ")\n")
  fr <- x$dif_results
  num <- vapply(fr, is.numeric, TRUE)
  fr[num] <- lapply(fr[num], round, 4)
  print(fr, row.names = FALSE)
  invisible(x)
}


#' @export
summary.dif_irt <- function(object, ...) {
  print(object)
  cat("\nUniform DIF effects (delta) are on the logit metric: the shift\n")
  cat("in the item's log-odds for a unit change in the covariate, at\n")
  cat("equal ability. Impact (gamma) is the latent-mean difference.\n")
  invisible(object)
}