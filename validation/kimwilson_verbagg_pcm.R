# Kim & Wilson (2019, Measurement 151:107062): polytomous item explanatory
# IRT with random item effects, on the verbal aggression data (3 categories:
# no / perhaps / yes). gllammr fit_eirt(model = "PCM") vs their published
# Stan posterior means.
#
#   K&W "MFRM"        = location-explanatory PCM, no item errors
#                       -> fit_eirt(..., model="PCM", item_residuals=FALSE)
#   K&W "MFRM + OIE"  = + overall item error on the location
#                       -> fit_eirt(..., model="PCM", item_residuals=TRUE)
#   K&W "LPCM + UISE" = step-explanatory PCM + univariate item-step errors
#                       -> fit_eirt(..., threshold_formula = ...)
#
# Their dummy coding: references are Want (mode), Self-to-blame (situ),
# Shout (btype), so relevel VerbAgg factors to match.

suppressMessages({
  library(gllammr)
  library(lme4)
})

data(VerbAgg, package = "lme4")

# ---- 3-category response matrix, coded 1/2/3 (no/perhaps/yes) ----
VerbAgg$y3 <- as.integer(VerbAgg$resp)   # ordered factor: no < perhaps < yes
resp3 <- with(VerbAgg, tapply(y3, list(id, item), identity))
resp3 <- matrix(as.integer(resp3), nrow = 316,
                dimnames = dimnames(with(VerbAgg, tapply(y3, list(id, item), identity))))
cat("Categories:", paste(sort(unique(as.vector(resp3))), collapse = "/"),
    "-", nrow(resp3), "persons x", ncol(resp3), "items\n")

# ---- Item design, releveled to Kim & Wilson's reference categories ----
item_info <- unique(VerbAgg[, c("item", "btype", "situ", "mode")])
item_info <- item_info[match(colnames(resp3), as.character(item_info$item)), ]
item_data <- data.frame(
  btype = relevel(factor(item_info$btype, ordered = FALSE), ref = "shout"),
  situ  = relevel(factor(item_info$situ,  ordered = FALSE), ref = "self"),
  mode  = factor(item_info$mode, ordered = FALSE),   # ref = want (default)
  row.names = colnames(resp3))

kw <- function(est, ours, se) {
  data.frame(KimWilson = est, gllammr = round(ours, 3), KW_SE = se)
}

# =====================================================================
# Model A: location-explanatory PCM, no item errors (K&W "MFRM")
# =====================================================================
cat("\n==== Model A: PCM location-explanatory, no item errors (K&W MFRM) ====\n")
tA <- system.time(
  mA <- fit_eirt(resp3, item_data,
                 difficulty_formula = ~ mode + situ + btype,
                 model = "PCM", item_residuals = FALSE)
)["elapsed"]
gA <- mA$regression_coefficients$difficulty
print(kw(c(1.58, 0.43, -0.82, -1.28, -0.63),
         gA[c("(Intercept)", "modedo", "situother", "btypecurse", "btypescold")],
         c(0.08, 0.04, 0.04, 0.05, 0.05)))
cat(sprintf("Person SD : K&W 0.95 | gllammr %.3f\n", mA$ability_sd))
cat(sprintf("logLik %.2f | time %.2fs\n", mA$logLik, tA))

# =====================================================================
# Model B: + overall item error (K&W "MFRM + OIE")
# =====================================================================
cat("\n==== Model B: PCM location-explanatory + item errors (K&W MFRM+OIE) ====\n")
tB <- system.time(
  mB <- fit_eirt(resp3, item_data,
                 difficulty_formula = ~ mode + situ + btype,
                 model = "PCM", item_residuals = TRUE)
)["elapsed"]
gB <- mB$regression_coefficients$difficulty
print(kw(c(1.69, 0.49, -0.89, -1.38, -0.70),
         gB[c("(Intercept)", "modedo", "situother", "btypecurse", "btypescold")],
         c(0.18, 0.15, 0.14, 0.17, 0.18)))
cat(sprintf("Item error SD : K&W 0.32 | gllammr %.3f\n",
            mB$residual_sd$difficulty))
cat(sprintf("Person SD     : K&W 0.97 | gllammr %.3f\n", mB$ability_sd))
cat(sprintf("logLik %.2f | time %.2fs\n", mB$logLik, tB))

# =====================================================================
# Descriptive PCM (saturated steps) + step-difficulty agreement
# (K&W Table 5: cor(PCM, MFRM) = 0.91, cor(PCM, MFRM+OIE) = 0.99)
# =====================================================================
cat("\n==== Descriptive PCM (fit_irt EM) and step-difficulty agreement ====\n")
tP <- system.time(mP <- fit_irt(resp3, model = "PCM"))["elapsed"]
cat(sprintf("PCM: logLik %.2f, person SD %.3f, time %.2fs\n",
            mP$logLik, mP$ability_sd, tP))

# Directly estimated step difficulties delta_im (2 per item)
delta_pcm <- do.call(rbind, mP$item_parameters$thresholds)

# Calculated step difficulties from the explanatory models:
# delta_im = b_i + s_im, s_i1 = step_param[i,1], s_i2 = -s_i1
calc_steps <- function(m) {
  pf <- m$tmb_obj$env$last.par.best
  s1 <- m$tmb_obj$env$parList(par = pf)$step_param[, 1]
  b  <- m$item_parameters$difficulty
  cbind(b + s1, b - s1)
}
dA <- calc_steps(mA); dB <- calc_steps(mB)
cat(sprintf("cor(PCM steps, MFRM calculated)     : %.3f  (K&W: 0.91)\n",
            cor(as.vector(delta_pcm), as.vector(dA))))
cat(sprintf("cor(PCM steps, MFRM+OIE calculated) : %.3f  (K&W: 0.99)\n",
            cor(as.vector(delta_pcm), as.vector(dB))))

# Step deviations vs K&W tau_i1 (Table 6, MFRM column)
pfA <- mA$tmb_obj$env$last.par.best
tauA <- mA$tmb_obj$env$parList(par = pfA)$step_param[, 1]
kw_tau <- c(-0.22, 0.00, -0.35, -0.47, -0.11, -0.11, -0.52, -0.29, -0.18,
            -0.62, -0.25, -0.19, -0.34, -0.25, -0.03, -0.17, -0.17, 0.30,
            -0.48, 0.01, 0.63, -0.60, -0.56, -0.02)
cat(sprintf("cor(gllammr tau_i1, K&W tau_i1)     : %.3f\n", cor(tauA, kw_tau)))

# =====================================================================
# Model C: step-explanatory LPCM + item-step errors (K&W "LPCM + UISE")
# =====================================================================
cat("\n==== Model C: LPCM + item-step errors (K&W LPCM+UISE) ====\n")
tC <- system.time(
  mC <- try(fit_eirt(resp3, item_data,
                     difficulty_formula = ~ 1,
                     threshold_formula = ~ mode + situ + btype,
                     model = "PCM", item_residuals = FALSE), silent = TRUE)
)["elapsed"]
if (!inherits(mC, "try-error")) {
  xi <- mC$regression_coefficients$threshold
  g0 <- mC$regression_coefficients$difficulty[["(Intercept)"]]
  # gamma0 and the per-step intercepts share a flat direction; only their
  # sum is identified -- report combined step intercepts
  step_int <- g0 + xi["(Intercept)", ]
  out <- rbind(
    `intercept`  = c(1.45, 1.90, round(step_int, 3)),
    `modedo`     = c(0.56, 0.42, round(xi["modedo", ], 3)),
    `situother`  = c(-0.70, -1.06, round(xi["situother", ], 3)),
    `btypecurse` = c(-1.71, -1.06, round(xi["btypecurse", ], 3)),
    `btypescold` = c(-0.84, -0.56, round(xi["btypescold", ], 3)))
  colnames(out) <- c("KW_step1", "KW_step2", "gllammr_step1", "gllammr_step2")
  print(out)
  cat(sprintf("Item-step error SD : K&W 0.33 | gllammr %.3f\n",
              mC$residual_sd$threshold))
  cat(sprintf("Person SD          : K&W 0.97 | gllammr %.3f\n", mC$ability_sd))
  cat(sprintf("logLik %.2f | time %.2fs\n", mC$logLik, tC))
} else {
  cat("LPCM fit failed:\n", attr(mC, "condition")$message, "\n")
}

cat("\nDone.\n")
