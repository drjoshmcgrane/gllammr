# Explanatory discrimination in the polytomous context: GPCM on the
# 3-category verbal aggression data with BOTH difficulty and discrimination
# regressed on the item design facets (log-linear: a_j = exp(W delta + eps_a)).
# Cross-check: two-stage analysis (descriptive GPCM, then lm of log(a_hat)
# on the facets) should approximately reproduce the one-stage estimates.

suppressMessages({
  library(GLLAMMR)
  library(lme4)
})

data(VerbAgg, package = "lme4")
VerbAgg$y3 <- as.integer(VerbAgg$resp)
resp3 <- with(VerbAgg, tapply(y3, list(id, item), identity))
resp3 <- matrix(as.integer(resp3), nrow = 316, dimnames = dimnames(resp3))
item_info <- unique(VerbAgg[, c("item", "btype", "situ", "mode")])
item_info <- item_info[match(colnames(resp3), as.character(item_info$item)), ]
item_data <- data.frame(
  btype = relevel(factor(item_info$btype, ordered = FALSE), ref = "shout"),
  situ  = relevel(factor(item_info$situ,  ordered = FALSE), ref = "self"),
  mode  = factor(item_info$mode, ordered = FALSE))

cat("==== One-stage: GPCM, difficulty AND discrimination explanatory ====\n")
t1 <- system.time(
  m <- fit_eirt(resp3, item_data,
                difficulty_formula     = ~ mode + situ + btype,
                discrimination_formula = ~ mode + situ + btype,
                model = "GPCM", item_residuals = TRUE)
)["elapsed"]
cat("\nDifficulty effects (gamma):\n")
print(round(m$regression_coefficients$difficulty, 3))
cat("\nLog-discrimination effects (delta):\n")
print(round(m$regression_coefficients$discrimination, 3))
cat("\nMultiplicative discrimination ratios exp(delta):\n")
print(round(exp(m$regression_coefficients$discrimination), 3))
cat(sprintf("\nResidual SDs: difficulty %.3f | log-discrimination %.3f\n",
            m$residual_sd$difficulty, m$residual_sd$discrimination))
cat(sprintf("sigma_theta (fixed): %.3f | logLik %.2f | time %.2fs\n",
            m$ability_sd, m$logLik, t1))
sdr_ok <- !inherits(m$tmb_sdr, "try-error")
if (sdr_ok) {
  s <- summary(m$tmb_sdr, "fixed")
  d_se <- s[rownames(s) == "delta", "Std. Error"]
  cat("delta SEs:", paste(round(d_se, 3), collapse = " "), "\n")
}

cat("\n==== Two-stage check: descriptive GPCM, then lm(log a ~ facets) ====\n")
t2 <- system.time(mg <- fit_irt(resp3, model = "GPCM"))["elapsed"]
a_hat <- mg$item_parameters$discrimination
two_stage <- lm(log(a_hat) ~ mode + situ + btype, data = item_data)
cmp <- cbind(one_stage = m$regression_coefficients$discrimination,
             two_stage = coef(two_stage))
print(round(cmp, 3))
cat(sprintf("cor(one-stage item discriminations, descriptive GPCM): %.3f\n",
            cor(m$item_parameters$discrimination, a_hat)))
cat(sprintf("Descriptive GPCM: logLik %.2f, time %.2fs\n", mg$logLik, t2))

cat("\n==== Same idea under the cumulative-logit GRM ====\n")
mgrm <- fit_eirt(resp3, item_data,
                 difficulty_formula     = ~ mode + situ + btype,
                 discrimination_formula = ~ mode,
                 model = "GRM", item_residuals = TRUE)
cat("GRM log-discrimination effects:\n")
print(round(mgrm$regression_coefficients$discrimination, 3))
cat(sprintf("logLik %.2f\n", mgrm$logLik))

cat("\nDone.\n")
