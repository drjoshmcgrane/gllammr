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
#'   "lca_carcinoma", "grm_science", "gamma_simulated",
#'   "survival_exponential", "sem_lavaan", "lca_polytomous",
#'   "npml_binomial", "aghq_binomial", "twopl_lsat_em".
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
                 "gamma_simulated", "survival_exponential", "sem_lavaan",
                 "lca_polytomous", "npml_binomial", "aghq_binomial",
                 "twopl_lsat_em")
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

  # Pin the Laplace path: this case validates exact agreement with the
  # identical approximation in glmer (nAGQ = 1)
  fit <- fit_irt(resp, model = "Rasch", method = "laplace")

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

  rbind(
    .val_row("sem_lavaan", "loading_x2",
             fit$loadings["x2", "f1"],
             pe$est[pe$lhs == "f1" & pe$op == "=~" & pe$rhs == "x2"], 5e-3),
    .val_row("sem_lavaan", "loading_y3",
             fit$loadings["y3", "f2"],
             pe$est[pe$lhs == "f2" & pe$op == "=~" & pe$rhs == "y3"], 5e-3),
    .val_row("sem_lavaan", "structural_f2_f1",
             fit$structural["f2", "f1"], pe$est[pe$op == "~"], 5e-3)
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
