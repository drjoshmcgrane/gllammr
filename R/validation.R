#' Cross-package validation of gllammr estimates
#'
#' Fits canonical benchmark datasets with gllammr and with established
#' reference packages, and reports the agreement. Reference packages that use
#' the same Laplace approximation (lme4 with nAGQ = 1, ordinal::clmm) should
#' agree to numerical precision; packages using different integration schemes
#' (mirt, ltm EM quadrature) agree within small tolerances.
#'
#' All reference packages are Suggests; cases whose reference package is not
#' installed are skipped.
#'
#' @param cases Character vector of case names to run, or "all" (default).
#'   Available: "gaussian_sleepstudy", "binomial_toenail",
#'   "poisson_grouseticks", "ordinal_wine", "rasch_lsat", "twopl_simulated",
#'   "lca_carcinoma", "grm_science", "gamma_simulated",
#'   "survival_exponential", "sem_lavaan", "lca_polytomous",
#'   "npml_binomial", "aghq_binomial", "twopl_lsat_em", "eirt_verbagg",
#'   "eirt_verbagg_pcm", "cdm_fraction_dina", "ordinal_crossed",
#'   "dif_logistic", "dif_irt_glmm".
#' @param scale "standard" (default) runs the canonical-dataset cases;
#'   "large" runs the large-scale tier (n in the tens of thousands, long
#'   item batteries - sizes where quadrature grids and tolerances can fail
#'   silently); "all" runs both.
#' @param verbose Print progress messages (default TRUE)
#'
#' @return Data frame with one row per compared statistic: case, statistic,
#'   gllammr value, reference value, absolute and relative difference,
#'   tolerance, and pass/fail.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   gllammr_validate(cases = "gaussian_sleepstudy")
#' }
#' }
#'
#' @export
gllammr_validate <- function(cases = "all", scale = c("standard", "large", "all"),
                             verbose = TRUE) {
  scale <- match.arg(scale)
  standard_cases <- c("gaussian_sleepstudy", "binomial_toenail",
                 "poisson_grouseticks", "ordinal_wine", "rasch_lsat",
                 "twopl_simulated", "lca_carcinoma", "grm_science",
                 "gamma_simulated", "survival_exponential", "sem_lavaan",
                 "lca_polytomous", "npml_binomial", "aghq_binomial",
                 "twopl_lsat_em", "eirt_verbagg", "eirt_verbagg_pcm",
                 "cdm_fraction_dina", "ordinal_crossed", "dif_logistic",
                 "dif_irt_glmm")
  # Large-scale tier: numerical behavior at sizes where quadrature grids,
  # tolerances, and interpreted-loop costs can fail silently
  large_cases <- c("large_glmm_binomial", "large_grm_battery",
                   "large_lca", "large_sem")
  all_cases <- switch(scale,
                      standard = standard_cases,
                      large = large_cases,
                      all = c(standard_cases, large_cases))
  if (identical(cases, "all")) cases <- all_cases
  unknown <- setdiff(cases, all_cases)
  if (length(unknown) > 0) {
    stop("Unknown validation case(s): ", paste(unknown, collapse = ", "))
  }

  rows <- list()
  for (case in cases) {
    if (verbose) message("Validating: ", case)
    fn <- get(paste0(".validate_", case), mode = "function")
    res <- tryCatch(fn(),
      # A numerical breakdown *inside a reference package* (e.g. lme4's
      # "Downdated VtV is not positive definite" on some BLAS/Matrix builds)
      # is a platform artefact, not a gllammr defect: record it as a skip
      # (pass = NA), not a failure.
      gllammr_reference_skip = function(e) {
        if (verbose) message("  skipped: ", conditionMessage(e))
        data.frame(case = case, statistic = "SKIP",
                   gllammr = NA_real_, reference = NA_real_,
                   abs_diff = NA_real_, rel_diff = NA_real_,
                   tolerance = NA_real_, pass = NA,
                   note = conditionMessage(e))
      },
      error = function(e) {
        data.frame(case = case, statistic = "ERROR",
                   gllammr = NA_real_, reference = NA_real_,
                   abs_diff = NA_real_, rel_diff = NA_real_,
                   tolerance = NA_real_, pass = FALSE,
                   note = conditionMessage(e))
      })
    rows[[case]] <- res
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' Signal that a reference-package fit failed (a platform artefact, not a
#' gllammr defect). Caught by [gllammr_validate()] and turned into a skipped
#' row rather than a failure.
#' @keywords internal
.reference_skip <- function(msg) {
  stop(structure(
    class = c("gllammr_reference_skip", "error", "condition"),
    list(message = msg, call = NULL)))
}

#' Evaluate a reference-package model fit, converting a fitting error into a
#' platform skip. Numerical breakdowns inside a reference package (e.g. lme4's
#' "Downdated VtV is not positive definite" on some BLAS/Matrix builds) then
#' mark the case skipped, never failed. `expr` is evaluated lazily inside the
#' handler.
#' @keywords internal
.reference_fit <- function(expr) {
  tryCatch(expr, error = function(e)
    .reference_skip(paste0(
      "reference-package fit failed on this platform: ",
      conditionMessage(e))))
}


#' Build one validation result row
#' @keywords internal
.val_row <- function(case, statistic, gllammr, reference, tolerance,
                     relative = TRUE, note = "") {
  abs_diff <- abs(gllammr - reference)
  rel_diff <- abs_diff / max(abs(reference), .Machine$double.eps)
  pass <- if (relative) rel_diff <= tolerance else abs_diff <= tolerance
  data.frame(case = case, statistic = statistic,
             gllammr = gllammr, reference = reference,
             abs_diff = abs_diff, rel_diff = rel_diff,
             tolerance = tolerance, pass = pass, note = note)
}


#' @keywords internal
.validate_gaussian_sleepstudy <- function() {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(NULL)
  }
  data("sleepstudy", package = "lme4", envir = environment())

  fit <- gllamm(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  ref <- lme4::lmer(Reaction ~ Days + (Days | Subject),
                    data = sleepstudy, REML = FALSE)

  Sigma <- fit$coefficients$random_var[[1]]
  Sigma_ref <- as.matrix(Matrix::bdiag(lme4::VarCorr(ref)$Subject))

  rbind(
    .val_row("gaussian_sleepstudy", "beta_intercept",
             unname(coef(fit)$fixed[1]), unname(lme4::fixef(ref)[1]), 1e-4),
    .val_row("gaussian_sleepstudy", "beta_Days",
             unname(coef(fit)$fixed[2]), unname(lme4::fixef(ref)[2]), 1e-4),
    .val_row("gaussian_sleepstudy", "logLik",
             fit$logLik, as.numeric(logLik(ref)), 1e-3, relative = FALSE),
    .val_row("gaussian_sleepstudy", "var_intercept",
             Sigma[1, 1], Sigma_ref[1, 1], 1e-2),
    .val_row("gaussian_sleepstudy", "var_slope",
             Sigma[2, 2], Sigma_ref[2, 2], 1e-2),
    .val_row("gaussian_sleepstudy", "cov_int_slope",
             Sigma[1, 2], Sigma_ref[1, 2], 5e-2)
  )
}


#' @keywords internal
.validate_binomial_toenail <- function() {
  if (!requireNamespace("lme4", quietly = TRUE) ||
      !requireNamespace("HSAUR3", quietly = TRUE)) {
    return(NULL)
  }
  data("toenail", package = "HSAUR3", envir = environment())
  toenail$y <- as.integer(toenail$outcome == "moderate or severe")

  fit <- gllamm(y ~ treatment * time + (1 | patientID), data = toenail,
                family = stats::binomial())
  ref <- lme4::glmer(y ~ treatment * time + (1 | patientID), data = toenail,
                     family = stats::binomial(), nAGQ = 1)

  sigma_u <- sqrt(fit$coefficients$random_var[[1]][1, 1])
  sigma_ref <- attr(lme4::VarCorr(ref)$patientID, "stddev")

  rbind(
    .val_row("binomial_toenail", "beta_time",
             unname(coef(fit)$fixed["time"]),
             unname(lme4::fixef(ref)["time"]), 1e-2),
    .val_row("binomial_toenail", "beta_interaction",
             unname(coef(fit)$fixed[4]), unname(lme4::fixef(ref)[4]), 2e-2),
    .val_row("binomial_toenail", "logLik",
             fit$logLik, as.numeric(logLik(ref)), 0.5, relative = FALSE),
    .val_row("binomial_toenail", "sigma_u",
             unname(sigma_u), unname(sigma_ref), 2e-2)
  )
}


#' @keywords internal
.validate_poisson_grouseticks <- function() {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(NULL)
  }
  data("grouseticks", package = "lme4", envir = environment())

  fit <- gllamm(TICKS ~ YEAR + (1 | BROOD), data = grouseticks,
                family = stats::poisson())
  ref <- lme4::glmer(TICKS ~ YEAR + (1 | BROOD), data = grouseticks,
                     family = stats::poisson(), nAGQ = 1)

  sigma_u <- sqrt(fit$coefficients$random_var[[1]][1, 1])
  sigma_ref <- attr(lme4::VarCorr(ref)$BROOD, "stddev")

  rbind(
    .val_row("poisson_grouseticks", "beta_intercept",
             unname(coef(fit)$fixed[1]), unname(lme4::fixef(ref)[1]), 1e-2),
    .val_row("poisson_grouseticks", "logLik",
             fit$logLik, as.numeric(logLik(ref)), 0.5, relative = FALSE),
    .val_row("poisson_grouseticks", "sigma_u",
             unname(sigma_u), unname(sigma_ref), 1e-2)
  )
}


#' @keywords internal
.validate_ordinal_wine <- function() {
  if (!requireNamespace("ordinal", quietly = TRUE)) {
    return(NULL)
  }
  data("wine", package = "ordinal", envir = environment())
  wine$rating_num <- as.integer(wine$rating)

  fit <- fit_ordinal(rating_num ~ temp + contact + (1 | judge), data = wine,
                     link = "logit")
  ref <- ordinal::clmm(rating ~ temp + contact + (1 | judge), data = wine,
                       link = "logit")

  sigma_u <- sqrt(fit$coefficients$random_var)
  sigma_ref <- as.numeric(attr(ordinal::VarCorr(ref)$judge, "stddev"))

  rbind(
    .val_row("ordinal_wine", "beta_temp",
             unname(fit$coefficients$fixed["tempwarm"]),
             unname(coef(ref)["tempwarm"]), 1e-2),
    .val_row("ordinal_wine", "beta_contact",
             unname(fit$coefficients$fixed["contactyes"]),
             unname(coef(ref)["contactyes"]), 1e-2),
    .val_row("ordinal_wine", "threshold_1",
             unname(fit$coefficients$thresholds[1]),
             unname(ref$alpha[1]), 2e-2),
    .val_row("ordinal_wine", "logLik",
             fit$logLik, as.numeric(logLik(ref)), 0.1, relative = FALSE),
    .val_row("ordinal_wine", "sigma_judge",
             unname(sigma_u), sigma_ref, 2e-2)
  )
}


#' @keywords internal
.validate_rasch_lsat <- function() {
  if (!requireNamespace("lme4", quietly = TRUE) ||
      !requireNamespace("ltm", quietly = TRUE)) {
    return(NULL)
  }
  data("LSAT", package = "ltm", envir = environment())
  resp <- as.matrix(LSAT)

  # Pin the Laplace path: this case validates exact agreement with the
  # identical approximation in glmer (nAGQ = 1)
  fit <- fit_irt(resp, model = "Rasch", method = "laplace", se = FALSE)

  # De Boeck-style Rasch as a GLMM: same model, same Laplace approximation
  long <- data.frame(
    y = as.vector(resp),
    item = factor(rep(colnames(resp), each = nrow(resp))),
    id = factor(rep(seq_len(nrow(resp)), times = ncol(resp)))
  )
  ref <- .reference_fit(lme4::glmer(
    y ~ 0 + item + (1 | id), data = long,
    family = stats::binomial(), nAGQ = 1,
    control = lme4::glmerControl(optimizer = "bobyqa")))
  b_ref <- -lme4::fixef(ref)   # difficulty = -item easiness
  sigma_ref <- attr(lme4::VarCorr(ref)$id, "stddev")

  rbind(
    .val_row("rasch_lsat", "difficulty_item1",
             unname(fit$item_parameters$difficulty[1]), unname(b_ref[1]), 1e-2),
    .val_row("rasch_lsat", "difficulty_item3",
             unname(fit$item_parameters$difficulty[3]), unname(b_ref[3]), 1e-2),
    .val_row("rasch_lsat", "sigma_theta",
             unname(fit$ability_sd), unname(sigma_ref), 1e-2),
    .val_row("rasch_lsat", "logLik",
             fit$logLik, as.numeric(logLik(ref)), 0.5, relative = FALSE)
  )
}


#' @keywords internal
.validate_twopl_simulated <- function() {
  if (!requireNamespace("mirt", quietly = TRUE)) {
    return(NULL)
  }
  # 20 items x 1000 persons: enough information per person for the Laplace
  # approximation to be accurate. NOTE: on very short tests (e.g. 5-item
  # LSAT) Laplace-based 2PL discrimination estimates can diverge - a known
  # limitation shared with joint ML, to be resolved by the adaptive
  # quadrature option. Documented in ?fit_irt.
  set.seed(7)
  np <- 1000; ni <- 20
  theta <- rnorm(np)
  a_true <- runif(ni, 0.6, 2.0)
  b_true <- rnorm(ni, 0, 1)
  p <- plogis(outer(theta, b_true, "-") * matrix(a_true, np, ni, byrow = TRUE))
  resp <- matrix(rbinom(np * ni, 1, p), np, ni)

  fit <- fit_irt(resp, model = "2PL", se = FALSE)
  ref <- mirt::mirt(as.data.frame(resp), 1, itemtype = "2PL", verbose = FALSE)
  co <- mirt::coef(ref, simplify = TRUE)$items
  b_ref <- -co[, "d"] / co[, "a1"]

  rbind(
    .val_row("twopl_simulated", "discrimination_cor",
             cor(fit$item_parameters$discrimination, co[, "a1"]), 1, 5e-3),
    .val_row("twopl_simulated", "difficulty_cor",
             cor(fit$item_parameters$difficulty, b_ref), 1, 5e-3),
    .val_row("twopl_simulated", "mean_abs_a_diff",
             mean(abs(fit$item_parameters$discrimination - co[, "a1"])), 0,
             0.05, relative = FALSE),
    .val_row("twopl_simulated", "mean_abs_b_diff",
             mean(abs(fit$item_parameters$difficulty - b_ref)), 0,
             0.05, relative = FALSE)
  )
}


#' @keywords internal
.validate_lca_carcinoma <- function() {
  if (!requireNamespace("poLCA", quietly = TRUE)) {
    return(NULL)
  }
  data("carcinoma", package = "poLCA", envir = environment())
  resp <- as.matrix(carcinoma) - 1L   # poLCA codes 1/2; we need 0/1

  fit <- fit_lca(resp, nclass = 2, control = list(n_starts = 5))

  f <- stats::as.formula(paste0(
    "cbind(", paste(colnames(carcinoma), collapse = ","), ") ~ 1"))
  set.seed(1)
  ref <- poLCA::poLCA(f, data = carcinoma, nclass = 2, nrep = 5,
                      verbose = FALSE)

  # Same likelihood, both should reach the global optimum
  rbind(
    .val_row("lca_carcinoma", "logLik",
             fit$logLik, ref$llik, 0.1, relative = FALSE),
    .val_row("lca_carcinoma", "max_class_proportion",
             max(fit$class_probs), max(ref$P), 2e-2)
  )
}


#' @keywords internal
.validate_grm_science <- function() {
  if (!requireNamespace("mirt", quietly = TRUE)) {
    return(NULL)
  }
  data("Science", package = "mirt", envir = environment())
  resp <- as.matrix(Science)

  fit <- fit_irt(resp, model = "GRM", se = FALSE)
  ref <- mirt::mirt(Science, 1, itemtype = "graded", verbose = FALSE)

  # mirt: a*theta + d_k; threshold b_k = -d_k / a (theta ~ N(0,1) fixed).
  # gllammr estimates sigma_theta, so compare on the standardized scale:
  # b* = b / sigma_theta, a* = a * sigma_theta.
  co <- mirt::coef(ref, simplify = TRUE)$items
  b_ref_item1 <- -co[1, "d1"] / co[1, "a1"]
  a_ref_item1 <- co[1, "a1"]

  sigma <- fit$ability_sd
  b_fit_item1 <- fit$item_parameters$thresholds[[1]][1] / sigma
  a_fit_item1 <- fit$item_parameters$discrimination[1] * sigma

  rbind(
    .val_row("grm_science", "threshold1_item1_std",
             b_fit_item1, b_ref_item1, 0.10),
    .val_row("grm_science", "discrimination_item1_std",
             a_fit_item1, a_ref_item1, 0.15)
  )
}


#' @keywords internal
.validate_gamma_simulated <- function() {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    return(NULL)
  }
  set.seed(41)
  n <- 2000; g <- 50
  grp <- factor(rep(1:g, each = n %/% g))
  x <- rnorm(n)
  u <- rnorm(g, 0, 0.5)
  mu <- exp(0.5 + 0.3 * x + u[as.integer(grp)])
  d <- data.frame(y = rgamma(n, shape = 2.5, scale = mu * 0.4),
                  x = x, grp = grp)

  fit <- gllamm(y ~ x + (1 | grp), data = d,
                family = stats::Gamma(link = "log"))
  ref <- glmmTMB::glmmTMB(y ~ x + (1 | grp), data = d,
                          family = stats::Gamma(link = "log"))

  rbind(
    .val_row("gamma_simulated", "beta_x",
             unname(coef(fit)$fixed[2]),
             unname(glmmTMB::fixef(ref)$cond[2]), 1e-3),
    .val_row("gamma_simulated", "logLik",
             fit$logLik, as.numeric(logLik(ref)), 0.1, relative = FALSE),
    .val_row("gamma_simulated", "sigma_u",
             sqrt(fit$coefficients$random_var[[1]][1, 1]),
             unname(attr(glmmTMB::VarCorr(ref)$cond$grp, "stddev")), 1e-2)
  )
}


#' @keywords internal
.validate_survival_exponential <- function() {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(NULL)
  }
  # The exponential frailty model is likelihood-equivalent to a Poisson GLMM
  # on the event indicator with offset log(t): an exact cross-check.
  set.seed(51)
  n <- 1500; g <- 50
  grp <- factor(rep(1:g, each = n %/% g))
  x <- rnorm(n)
  u <- rnorm(g, 0, 0.6)
  eta <- -1 + 0.5 * x + u[as.integer(grp)]
  t_true <- rexp(n, rate = exp(eta))
  cens <- rexp(n, rate = 0.15)
  d <- data.frame(time = pmin(t_true, cens),
                  status = as.integer(t_true <= cens), x = x, grp = grp)

  fit <- fit_survival(Surv(time, status) ~ x + (1 | grp), data = d,
                      distribution = "exponential")
  ref <- lme4::glmer(status ~ x + offset(log(time)) + (1 | grp), data = d,
                     family = stats::poisson(), nAGQ = 1)

  rbind(
    .val_row("survival_exponential", "beta_x",
             unname(coef(fit)$fixed[2]), unname(lme4::fixef(ref)[2]), 1e-3),
    .val_row("survival_exponential", "sigma_frailty",
             unname(fit$coefficients$random_sd),
             unname(attr(lme4::VarCorr(ref)$grp, "stddev")), 1e-3),
    .val_row("survival_exponential", "logLik_vs_poisson_plus_constant",
             fit$logLik,
             as.numeric(logLik(ref)) - sum(d$status * log(d$time)),
             0.01, relative = FALSE)
  )
}


#' @keywords internal
.validate_sem_lavaan <- function() {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    return(NULL)
  }
  set.seed(71)
  n <- 800
  f1 <- rnorm(n)
  f2 <- 0.6 * f1 + rnorm(n, 0, 0.8)
  d <- data.frame(
    x1 = 1.0 + 1.0 * f1 + rnorm(n, 0, 0.6),
    x2 = 0.5 + 0.8 * f1 + rnorm(n, 0, 0.6),
    x3 = -0.3 + 1.2 * f1 + rnorm(n, 0, 0.6),
    y1 = 0.2 + 1.0 * f2 + rnorm(n, 0, 0.5),
    y2 = 0.0 + 0.9 * f2 + rnorm(n, 0, 0.5),
    y3 = 0.8 + 1.1 * f2 + rnorm(n, 0, 0.5)
  )
  fit <- fit_sem(measurement = list(f1 = ~ x1 + x2 + x3, f2 = ~ y1 + y2 + y3),
                 structural = list(f2 ~ f1), data = d)
  lav <- lavaan::sem("f1 =~ x1 + x2 + x3\nf2 =~ y1 + y2 + y3\nf2 ~ f1",
                     data = d)
  pe <- lavaan::parameterEstimates(lav)

  # FIML under MCAR missingness, vs lavaan missing = "fiml"
  d_na <- d
  set.seed(72)
  for (v in c("x1", "y2")) d_na[[v]][sample(nrow(d), 80)] <- NA
  fit_f <- fit_sem(measurement = list(f1 = ~ x1 + x2 + x3,
                                      f2 = ~ y1 + y2 + y3),
                   structural = list(f2 ~ f1), data = d_na,
                   missing = "fiml")
  lav_f <- lavaan::sem("f1 =~ x1 + x2 + x3\nf2 =~ y1 + y2 + y3\nf2 ~ f1",
                       data = d_na, missing = "fiml", fixed.x = FALSE)
  pe_f <- lavaan::parameterEstimates(lav_f)

  i_se <- fit$param_table$label == "f2~f1"
  j_se <- pe$lhs == "f2" & pe$op == "~" & pe$rhs == "f1"

  rbind(
    .val_row("sem_lavaan", "loading_x2",
             fit$loadings["x2", "f1"],
             pe$est[pe$lhs == "f1" & pe$op == "=~" & pe$rhs == "x2"], 5e-3),
    .val_row("sem_lavaan", "loading_y3",
             fit$loadings["y3", "f2"],
             pe$est[pe$lhs == "f2" & pe$op == "=~" & pe$rhs == "y3"], 5e-3),
    .val_row("sem_lavaan", "structural_f2_f1",
             fit$structural["f2", "f1"],
             pe$est[pe$op == "~"], 5e-3),
    .val_row("sem_lavaan", "se_structural_f2_f1",
             fit$param_table$se[i_se], pe$se[j_se], 2e-2),
    .val_row("sem_lavaan", "cfi",
             unname(fit$fit_measures["cfi"]),
             as.numeric(lavaan::fitMeasures(lav, "cfi")), 1e-3),
    .val_row("sem_lavaan", "rmsea",
             unname(fit$fit_measures["rmsea"]),
             as.numeric(lavaan::fitMeasures(lav, "rmsea")), 0.01,
             relative = FALSE),
    .val_row("sem_lavaan", "fiml_logLik",
             fit_f$logLik, as.numeric(lavaan::fitMeasures(lav_f, "logl")),
             0.01, relative = FALSE),
    .val_row("sem_lavaan", "fiml_structural_f2_f1",
             fit_f$structural["f2", "f1"],
             pe_f$est[pe_f$lhs == "f2" & pe_f$op == "~"], 5e-3)
  )
}


#' @keywords internal
.validate_lca_polytomous <- function() {
  if (!requireNamespace("poLCA", quietly = TRUE)) {
    return(NULL)
  }
  set.seed(101)
  n <- 800
  cls <- sample(1:2, n, TRUE, prob = c(0.6, 0.4))
  P1 <- list(c(.7, .2, .1), c(.1, .3, .6))
  Y <- sapply(1:4, function(j) {
    sapply(cls, function(k) sample(1:3, 1, prob = P1[[k]]))
  })
  colnames(Y) <- paste0("V", 1:4)

  fit <- fit_lca(Y, nclass = 2, control = list(n_starts = 5))
  dl <- as.data.frame(Y)
  f <- stats::as.formula(paste0("cbind(", paste(names(dl), collapse = ","),
                                ") ~ 1"))
  set.seed(1)
  ref <- poLCA::poLCA(f, data = dl, nclass = 2, nrep = 5, verbose = FALSE)

  rbind(
    .val_row("lca_polytomous", "logLik",
             fit$logLik, ref$llik, 0.1, relative = FALSE),
    .val_row("lca_polytomous", "max_class_proportion",
             max(fit$class_probs), max(ref$P), 1e-2)
  )
}


#' @keywords internal
.validate_npml_binomial <- function() {
  if (!requireNamespace("npmlreg", quietly = TRUE)) {
    return(NULL)
  }
  set.seed(121)
  g <- 80; n_per <- 10; n <- g * n_per
  grp <- factor(rep(1:g, each = n_per))
  x <- rnorm(n)
  cls <- sample(1:2, g, TRUE, prob = c(0.6, 0.4))
  locs_true <- c(-1, 1.5)
  d <- data.frame(x = x, grp = grp,
                  yb = rbinom(n, 1, plogis(locs_true[cls[as.integer(grp)]] +
                                             0.5 * x)))

  fit <- fit_npml(yb ~ x + (1 | grp), data = d, k = 2,
                  family = stats::binomial())
  ref <- suppressMessages(
    npmlreg::allvc(yb ~ x, random = ~ 1 | grp, data = d, k = 2,
                   family = binomial(), verbose = FALSE, plot.opt = 0))
  ref_locs <- sort(ref$mass.points)
  ref_masses <- ref$masses[order(ref$mass.points)]

  rbind(
    .val_row("npml_binomial", "beta_x",
             unname(coef(fit)$fixed["x"]), unname(coef(ref)["x"]), 5e-3),
    .val_row("npml_binomial", "location_1",
             fit$locations[1], unname(ref_locs[1]), 1e-2),
    .val_row("npml_binomial", "location_2",
             fit$locations[2], unname(ref_locs[2]), 1e-2),
    .val_row("npml_binomial", "mass_1",
             fit$masses[1], unname(ref_masses[1]), 1e-2),
    .val_row("npml_binomial", "logLik",
             fit$logLik, -ref$disparity / 2, 0.1, relative = FALSE)
  )
}


#' @keywords internal
.validate_aghq_binomial <- function() {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(NULL)
  }
  # Small clusters + large sigma_u: where Laplace is weakest and adaptive
  # quadrature matters. glmer with the same nAGQ is the reference.
  set.seed(131)
  g <- 100; n_per <- 6; n <- g * n_per
  grp <- factor(rep(1:g, each = n_per))
  x <- rnorm(n)
  u <- rnorm(g, 0, 2)
  d <- data.frame(x = x, grp = grp,
                  yb = rbinom(n, 1, plogis(-0.5 + 0.8 * x + u[as.integer(grp)])))

  fit <- gllamm(yb ~ x + (1 | grp), data = d, family = stats::binomial(),
                integration = aghq(15))
  ref <- lme4::glmer(yb ~ x + (1 | grp), data = d,
                     family = stats::binomial(), nAGQ = 15)

  rbind(
    .val_row("aghq_binomial", "beta_x",
             unname(coef(fit)$fixed[2]), unname(lme4::fixef(ref)[2]), 2e-3),
    .val_row("aghq_binomial", "sigma_u",
             sqrt(fit$coefficients$random_var[[1]][1, 1]),
             unname(attr(lme4::VarCorr(ref)$grp, "stddev")), 5e-3),
    .val_row("aghq_binomial", "logLik",
             fit$logLik, as.numeric(logLik(ref)), 0.05, relative = FALSE)
  )
}


#' @keywords internal
.validate_twopl_lsat_em <- function() {
  if (!requireNamespace("ltm", quietly = TRUE)) {
    return(NULL)
  }
  # The 5-item LSAT 2PL diverges under joint Laplace (documented); the EM
  # path handles it like ltm's EM does.
  data("LSAT", package = "ltm", envir = environment())
  fit <- fit_irt(as.matrix(LSAT), model = "2PL", method = "em")
  ref_coef <- coef(ltm::ltm(LSAT ~ z1))

  rbind(
    .val_row("twopl_lsat_em", "difficulty_item1",
             unname(fit$item_parameters$difficulty[1]),
             unname(ref_coef[1, "Dffclt"]), 0.02),
    .val_row("twopl_lsat_em", "discrimination_item1",
             unname(fit$item_parameters$discrimination[1]),
             unname(ref_coef[1, "Dscrmn"]), 0.03),
    .val_row("twopl_lsat_em", "difficulty_mean",
             mean(fit$item_parameters$difficulty),
             mean(ref_coef[, "Dffclt"]), 0.02)
  )
}


#' @keywords internal
.validate_eirt_verbagg <- function() {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(NULL)
  }
  # De Boeck & Wilson (2004): LLTM + error on the verbal aggression data.
  # The same model as glmer(r2 ~ btype + situ + mode + (1|id) + (1|item))
  # under the same Laplace approximation, so agreement is tight. fit_eirt
  # models difficulty (theta - b), glmer models easiness: gamma = -beta.
  data("VerbAgg", package = "lme4", envir = environment())
  VerbAgg$y <- as.integer(VerbAgg$r2 == "Y")
  resp <- with(VerbAgg, tapply(y, list(id, item), identity))
  resp <- matrix(as.integer(resp), nrow = nrow(resp),
                 dimnames = dimnames(resp))
  item_info <- unique(VerbAgg[, c("item", "btype", "situ", "mode")])
  item_info <- item_info[match(colnames(resp), as.character(item_info$item)), ]
  item_data <- data.frame(btype = factor(item_info$btype, ordered = FALSE),
                          situ = factor(item_info$situ, ordered = FALSE),
                          mode = factor(item_info$mode, ordered = FALSE))

  fit <- fit_eirt(resp, item_data,
                  difficulty_formula = ~ btype + situ + mode,
                  model = "Rasch", item_residuals = TRUE)
  ref <- .reference_fit(lme4::glmer(
    r2 ~ btype + situ + mode + (1 | id) + (1 | item),
    data = VerbAgg, family = stats::binomial(), nAGQ = 1,
    control = lme4::glmerControl(optimizer = "bobyqa")))
  beta_ref <- lme4::fixef(ref)
  gamma <- fit$regression_coefficients$difficulty

  rbind(
    .val_row("eirt_verbagg", "gamma_btype_shout",
             unname(gamma["btypeshout"]), -unname(beta_ref["btypeshout"]),
             1e-2),
    .val_row("eirt_verbagg", "gamma_situ_self",
             unname(gamma["situself"]), -unname(beta_ref["situself"]), 1e-2),
    .val_row("eirt_verbagg", "gamma_mode_do",
             unname(gamma["modedo"]), -unname(beta_ref["modedo"]), 1e-2),
    .val_row("eirt_verbagg", "sigma_theta",
             unname(fit$ability_sd),
             unname(attr(lme4::VarCorr(ref)$id, "stddev")), 1e-2),
    .val_row("eirt_verbagg", "sigma_item_residual",
             unname(fit$residual_sd$difficulty),
             unname(attr(lme4::VarCorr(ref)$item, "stddev")), 1e-2),
    .val_row("eirt_verbagg", "logLik",
             fit$logLik, as.numeric(logLik(ref)), 0.1, relative = FALSE)
  )
}


#' @keywords internal
.validate_eirt_verbagg_pcm <- function() {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(NULL)
  }
  # Kim & Wilson (2019, Measurement 151:107062): polytomous item explanatory
  # models on the 3-category verbal aggression responses. References are
  # their published Stan posterior means (Tables 5-6), so tolerances allow
  # for Bayesian-vs-ML differences. Their dummy coding: references Want
  # (mode), Self-to-blame (situ), Shout (btype).
  data("VerbAgg", package = "lme4", envir = environment())
  VerbAgg$y3 <- as.integer(VerbAgg$resp)
  resp3 <- with(VerbAgg, tapply(y3, list(id, item), identity))
  resp3 <- matrix(as.integer(resp3), nrow = nrow(resp3),
                  dimnames = dimnames(resp3))
  item_info <- unique(VerbAgg[, c("item", "btype", "situ", "mode")])
  item_info <- item_info[match(colnames(resp3),
                               as.character(item_info$item)), ]
  item_data <- data.frame(
    btype = stats::relevel(factor(item_info$btype, ordered = FALSE),
                           ref = "shout"),
    situ = stats::relevel(factor(item_info$situ, ordered = FALSE),
                          ref = "self"),
    mode = factor(item_info$mode, ordered = FALSE))

  # Their "MFRM + OIE": location-explanatory PCM with random item errors
  fit <- fit_eirt(resp3, item_data,
                  difficulty_formula = ~ mode + situ + btype,
                  model = "PCM", item_residuals = TRUE)
  gamma <- fit$regression_coefficients$difficulty

  # Their Table 5 agreement check: step difficulties calculated from the
  # explanatory model vs the directly estimated (descriptive) PCM
  pcm <- fit_irt(resp3, model = "PCM", se = FALSE)
  delta_pcm <- do.call(rbind, pcm$item_parameters$thresholds)
  pf <- fit$tmb_obj$env$last.par.best
  s1 <- fit$tmb_obj$env$parList(par = pf)$step_param[, 1]
  b <- fit$item_parameters$difficulty
  step_cor <- stats::cor(as.vector(delta_pcm), as.vector(cbind(b + s1, b - s1)))

  rbind(
    .val_row("eirt_verbagg_pcm", "gamma_intercept",
             unname(gamma["(Intercept)"]), 1.69, 0.05, relative = FALSE),
    .val_row("eirt_verbagg_pcm", "gamma_mode_do",
             unname(gamma["modedo"]), 0.49, 0.05, relative = FALSE),
    .val_row("eirt_verbagg_pcm", "gamma_situ_other",
             unname(gamma["situother"]), -0.89, 0.05, relative = FALSE),
    .val_row("eirt_verbagg_pcm", "gamma_btype_curse",
             unname(gamma["btypecurse"]), -1.38, 0.05, relative = FALSE),
    .val_row("eirt_verbagg_pcm", "gamma_btype_scold",
             unname(gamma["btypescold"]), -0.70, 0.05, relative = FALSE),
    .val_row("eirt_verbagg_pcm", "sigma_theta",
             unname(fit$ability_sd), 0.97, 0.05, relative = FALSE),
    .val_row("eirt_verbagg_pcm", "step_difficulty_cor_vs_pcm",
             step_cor, 0.99, 0.01, relative = FALSE)
  )
}


#' @keywords internal
.validate_dif_logistic <- function() {
  if (!requireNamespace("difR", quietly = TRUE)) {
    return(NULL)
  }
  # Logistic-regression DIF with score matching and purification: with
  # the same (observed-score) matching criterion, gllammr's tests are the
  # same nested-model LR tests as difR::difLogistic.
  set.seed(19)
  n <- 1200; ni <- 12
  theta <- rnorm(n)
  g <- rep(c(0, 1), length.out = n)
  b <- seq(-1.5, 1.5, length.out = ni)
  resp <- sapply(seq_len(ni), function(j) {
    eta <- theta - b[j]
    if (j %in% c(3, 7)) eta <- eta - 0.8 * g
    rbinom(n, 1, plogis(eta))
  })
  grp <- factor(ifelse(g == 1, "B", "A"))

  fit <- dif_test(resp, dif = grp, match = "score", type = "both",
                  purify = TRUE)
  ref <- difR::difLogistic(as.data.frame(resp), group = grp,
                           focal.name = "B", type = "both", purify = TRUE)
  ref_flagged <- ref$DIFitems
  if (identical(ref_flagged, "No DIF item detected")) {
    ref_flagged <- integer(0)
  }

  rbind(
    .val_row("dif_logistic", "flag_agreement",
             as.numeric(setequal(fit$flagged_items, ref_flagged)), 1,
             1e-9, relative = FALSE),
    .val_row("dif_logistic", "stat_rank_correlation",
             stats::cor(fit$dif_results$chisq, ref$Logistik,
                        method = "spearman"), 1, 0.05, relative = FALSE),
    .val_row("dif_logistic", "stat_item3",
             fit$dif_results$chisq[fit$dif_results$item == 3],
             unname(ref$Logistik[3]), 0.15)
  )
}


#' @keywords internal
.validate_dif_irt_glmm <- function() {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(NULL)
  }
  # Confirmatory IRT-LR DIF: for the Rasch model with uniform DIF, the
  # compared marginal likelihoods are exactly those of the long-format
  # GLMM y ~ 0 + item + z + item_j:z + (1|person) under the same Laplace
  # approximation (the De Boeck & Wilson formulation).
  set.seed(55)
  n <- 600; ni <- 8
  g <- rep(c(0, 1), length.out = n)
  theta <- rnorm(n, mean = -0.5 * g)
  b <- seq(-1.5, 1.5, length.out = ni)
  resp <- sapply(seq_len(ni), function(j) {
    eta <- theta - b[j]
    if (j == 4) eta <- eta - 0.8 * g
    rbinom(n, 1, plogis(eta))
  })
  grp <- factor(ifelse(g == 1, "B", "A"))

  fit <- dif_irt(resp, dif = grp, items = 4,
                 anchors = setdiff(1:ni, 4), model = "Rasch")

  long <- data.frame(y = as.vector(resp),
                     item = factor(rep(1:ni, each = n)),
                     id = factor(rep(1:n, times = ni)),
                     g = rep(g, ni))
  long$dif4 <- as.integer(long$item == 4) * long$g
  ctrl <- lme4::glmerControl(optimizer = "bobyqa")
  m1 <- .reference_fit(lme4::glmer(
    y ~ 0 + item + g + dif4 + (1 | id), data = long,
    family = stats::binomial(), nAGQ = 1, control = ctrl))
  m0 <- .reference_fit(lme4::glmer(
    y ~ 0 + item + g + (1 | id), data = long,
    family = stats::binomial(), nAGQ = 1, control = ctrl))
  lr_glmer <- 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0)))

  rbind(
    .val_row("dif_irt_glmm", "lr_statistic_item4",
             fit$dif_results$chisq[1], lr_glmer, 0.05, relative = FALSE),
    .val_row("dif_irt_glmm", "delta_item4",
             fit$dif_results$delta_groupB[1],
             unname(lme4::fixef(m1)["dif4"]), 0.02, relative = FALSE),
    .val_row("dif_irt_glmm", "impact_gamma",
             fit$impact$gamma[1], unname(lme4::fixef(m1)["g"]), 0.05,
             relative = FALSE)
  )
}


#' @keywords internal
.validate_ordinal_crossed <- function() {
  if (!requireNamespace("ordinal", quietly = TRUE)) {
    return(NULL)
  }
  # Crossed random effects in a cumulative-logit model: judge and bottle
  # effects on the wine ratings. ordinal::clmm uses the same Laplace
  # approximation, so agreement is tight.
  data("wine", package = "ordinal", envir = environment())
  wine$rating_num <- as.integer(wine$rating)

  fit <- fit_ordinal(rating_num ~ temp + (1 | judge) + (1 | bottle),
                     data = wine, link = "logit")
  ref <- ordinal::clmm(rating ~ temp + (1 | judge) + (1 | bottle),
                       data = wine, link = "logit")

  rbind(
    .val_row("ordinal_crossed", "beta_temp",
             unname(fit$coefficients$fixed["tempwarm"]),
             unname(coef(ref)["tempwarm"]), 1e-2),
    .val_row("ordinal_crossed", "sigma_judge",
             unname(sqrt(fit$coefficients$random_var$judge[1, 1])),
             unname(attr(ordinal::VarCorr(ref)$judge, "stddev")), 2e-2),
    .val_row("ordinal_crossed", "sigma_bottle",
             unname(sqrt(fit$coefficients$random_var$bottle[1, 1])),
             unname(attr(ordinal::VarCorr(ref)$bottle, "stddev")), 5e-2),
    .val_row("ordinal_crossed", "logLik",
             fit$logLik, as.numeric(logLik(ref)), 0.05, relative = FALSE)
  )
}


#' @keywords internal
.validate_cdm_fraction_dina <- function() {
  if (!requireNamespace("CDM", quietly = TRUE)) {
    return(NULL)
  }
  # de la Torre (2009) fraction-subtraction data: 536 persons, 20 items,
  # 8 attributes (256 profiles). CDM::din is the standard DINA EM.
  # The profile-prevalence surface is near-flat (many more profiles than
  # the items can separate), so compare the identified quantities: logLik
  # and item guess/slip.
  data("fraction.subtraction.data", package = "CDM", envir = environment())
  data("fraction.subtraction.qmatrix", package = "CDM",
       envir = environment())
  Y <- as.matrix(fraction.subtraction.data)
  Q <- as.matrix(fraction.subtraction.qmatrix)

  fit <- fit_cdm(Y, Q, model = "dina", control = list(n_starts = 2))
  ref <- CDM::din(Y, q.matrix = Q, rule = "DINA", progress = FALSE)

  ghat <- vapply(fit$item_params, function(e) e$guess, 0)
  shat <- vapply(fit$item_params, function(e) e$slip, 0)

  rbind(
    .val_row("cdm_fraction_dina", "logLik",
             fit$logLik, ref$loglike, 0.5, relative = FALSE),
    .val_row("cdm_fraction_dina", "mean_abs_guess_diff",
             mean(abs(ghat - ref$guess$est)), 0, 0.01, relative = FALSE),
    .val_row("cdm_fraction_dina", "mean_abs_slip_diff",
             mean(abs(shat - ref$slip$est)), 0, 0.01, relative = FALSE),
    .val_row("cdm_fraction_dina", "guess_item1",
             unname(ghat[1]), unname(ref$guess$est[1]), 0.02,
             relative = FALSE)
  )
}


#' @keywords internal
.validate_large_glmm_binomial <- function() {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    return(NULL)
  }
  set.seed(301)
  n <- 100000; g <- 1000
  grp <- factor(rep(1:g, each = n %/% g))
  x <- rnorm(n)
  u <- rnorm(g, 0, 0.8)
  d <- data.frame(x = x, grp = grp,
                  yb = rbinom(n, 1, plogis(-0.3 + 0.5 * x + u[as.integer(grp)])))

  fit <- gllamm(yb ~ x + (1 | grp), data = d, family = stats::binomial())
  ref <- glmmTMB::glmmTMB(yb ~ x + (1 | grp), data = d,
                          family = stats::binomial())

  rbind(
    .val_row("large_glmm_binomial", "beta_x",
             unname(coef(fit)$fixed[2]),
             unname(glmmTMB::fixef(ref)$cond[2]), 1e-3),
    .val_row("large_glmm_binomial", "sigma_u",
             sqrt(fit$coefficients$random_var[[1]][1, 1]),
             unname(attr(glmmTMB::VarCorr(ref)$cond$grp, "stddev")), 1e-3),
    .val_row("large_glmm_binomial", "logLik",
             fit$logLik, as.numeric(logLik(ref)), 0.5, relative = FALSE)
  )
}


#' @keywords internal
.validate_large_grm_battery <- function() {
  if (!requireNamespace("mirt", quietly = TRUE)) {
    return(NULL)
  }
  set.seed(9)
  np <- 5000; ni <- 100
  theta <- rnorm(np)
  a <- runif(ni, 0.7, 1.8)
  taus <- t(sapply(rnorm(ni), function(b0) b0 + c(-1.5, -0.5, 0.5, 1.5)))
  resp <- sapply(1:ni, function(j) {
    cum <- sapply(1:4, function(k) plogis(a[j] * (theta - taus[j, k])))
    1L + rowSums(matrix(runif(np), np, 4) < cum)
  })

  fit <- fit_irt(resp, model = "GRM", se = FALSE)
  ref <- mirt::mirt(as.data.frame(resp), 1, itemtype = "graded",
                    verbose = FALSE)
  co <- mirt::coef(ref, simplify = TRUE)$items

  rbind(
    .val_row("large_grm_battery", "discrimination_cor",
             cor(fit$item_parameters$discrimination, co[, "a1"]), 1, 1e-4),
    .val_row("large_grm_battery", "logLik",
             fit$logLik, mirt::extract.mirt(ref, "logLik"),
             0.5, relative = FALSE)
  )
}


#' @keywords internal
.validate_large_lca <- function() {
  if (!requireNamespace("poLCA", quietly = TRUE)) {
    return(NULL)
  }
  set.seed(201)
  n <- 20000
  cls <- sample(1:3, n, TRUE, prob = c(0.5, 0.3, 0.2))
  pmat <- matrix(runif(3 * 8, 0.1, 0.9), 3, 8)
  Y <- matrix(rbinom(n * 8, 1, pmat[cls, ]), n, 8)

  fit <- fit_lca(Y, nclass = 3, control = list(n_starts = 3))
  dl <- as.data.frame(Y + 1)
  f <- stats::as.formula(paste0("cbind(", paste(names(dl), collapse = ","),
                                ") ~ 1"))
  set.seed(2)
  ref <- poLCA::poLCA(f, data = dl, nclass = 3, nrep = 3, verbose = FALSE)

  rbind(
    .val_row("large_lca", "logLik",
             fit$logLik, ref$llik, 0.5, relative = FALSE),
    .val_row("large_lca", "max_class_proportion",
             max(fit$class_probs), max(ref$P), 1e-2)
  )
}


#' @keywords internal
.validate_large_sem <- function() {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    return(NULL)
  }
  set.seed(71)
  n <- 100000
  f1 <- rnorm(n); f2 <- 0.6 * f1 + rnorm(n, 0, 0.8)
  d <- data.frame(
    x1 = 1.0 * f1 + rnorm(n, 0, .6), x2 = 0.8 * f1 + rnorm(n, 0, .6),
    x3 = 1.2 * f1 + rnorm(n, 0, .6),
    y1 = 1.0 * f2 + rnorm(n, 0, .5), y2 = 0.9 * f2 + rnorm(n, 0, .5),
    y3 = 1.1 * f2 + rnorm(n, 0, .5))

  fit <- fit_sem(measurement = list(f1 = ~ x1 + x2 + x3, f2 = ~ y1 + y2 + y3),
                 structural = list(f2 ~ f1), data = d)
  lav <- lavaan::sem("f1 =~ x1 + x2 + x3\nf2 =~ y1 + y2 + y3\nf2 ~ f1",
                     data = d)
  pe <- lavaan::parameterEstimates(lav)

  rbind(
    .val_row("large_sem", "structural_f2_f1",
             fit$structural["f2", "f1"], pe$est[pe$op == "~"], 1e-3),
    .val_row("large_sem", "loading_x3",
             fit$loadings["x3", "f1"],
             pe$est[pe$lhs == "f1" & pe$op == "=~" & pe$rhs == "x3"], 1e-3)
  )
}
