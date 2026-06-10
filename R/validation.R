#' Cross-package validation of GLLAMMR estimates
#'
#' Fits canonical benchmark datasets with GLLAMMR and with established
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
#'   "lca_carcinoma", "grm_science", "gamma_simulated".
#' @param verbose Print progress messages (default TRUE)
#'
#' @return Data frame with one row per compared statistic: case, statistic,
#'   GLLAMMR value, reference value, absolute and relative difference,
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
gllammr_validate <- function(cases = "all", verbose = TRUE) {
  all_cases <- c("gaussian_sleepstudy", "binomial_toenail",
                 "poisson_grouseticks", "ordinal_wine", "rasch_lsat",
                 "twopl_simulated", "lca_carcinoma", "grm_science",
                 "gamma_simulated")
  if (identical(cases, "all")) cases <- all_cases
  unknown <- setdiff(cases, all_cases)
  if (length(unknown) > 0) {
    stop("Unknown validation case(s): ", paste(unknown, collapse = ", "))
  }

  rows <- list()
  for (case in cases) {
    if (verbose) message("Validating: ", case)
    fn <- get(paste0(".validate_", case), mode = "function")
    res <- tryCatch(fn(), error = function(e) {
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

  fit <- fit_irt(resp, model = "Rasch")

  # De Boeck-style Rasch as a GLMM: same model, same Laplace approximation
  long <- data.frame(
    y = as.vector(resp),
    item = factor(rep(colnames(resp), each = nrow(resp))),
    id = factor(rep(seq_len(nrow(resp)), times = ncol(resp)))
  )
  ref <- lme4::glmer(y ~ 0 + item + (1 | id), data = long,
                     family = stats::binomial(), nAGQ = 1)
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

  fit <- fit_irt(resp, model = "2PL")
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

  fit <- fit_irt(resp, model = "GRM")
  ref <- mirt::mirt(Science, 1, itemtype = "graded", verbose = FALSE)

  # mirt: a*theta + d_k; threshold b_k = -d_k / a (theta ~ N(0,1) fixed).
  # GLLAMMR estimates sigma_theta, so compare on the standardized scale:
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
