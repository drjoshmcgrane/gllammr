# De Boeck & Wilson (2004) verbal aggression analysis:
# GLLAMMR fit_eirt() vs the lme4 item-explanatory GLMM formulation.

suppressMessages({
  library(GLLAMMR)
  library(lme4)
})

data(VerbAgg, package = "lme4")

# ---- Reshape to persons x items ----
VerbAgg$y <- as.integer(VerbAgg$r2 == "Y")
resp <- with(VerbAgg, tapply(y, list(id, item), identity))
resp <- matrix(as.integer(resp), nrow = nrow(resp),
               dimnames = dimnames(with(VerbAgg, tapply(y, list(id, item), identity))))
cat("Response matrix:", nrow(resp), "persons x", ncol(resp), "items;",
    sum(is.na(resp)), "missing\n")

# ---- Item design (one row per item, aligned with columns of resp) ----
item_info <- unique(VerbAgg[, c("item", "btype", "situ", "mode")])
item_info <- item_info[match(colnames(resp), as.character(item_info$item)), ]
item_data <- data.frame(btype = item_info$btype,
                        situ  = item_info$situ,
                        mode  = item_info$mode,
                        row.names = colnames(resp))
stopifnot(nrow(item_data) == 24)

med3 <- function(f) {
  t <- replicate(3, system.time(f())["elapsed"])
  median(t)
}

# =====================================================================
# Model 1: LLTM ("pure" -- item difficulties fully explained by design)
#   glmer:  r2 ~ btype + situ + mode + (1|id)
#   eirt :  difficulty_formula = ~ btype + situ + mode, item_residuals = FALSE
# =====================================================================
cat("\n==== Model 1: LLTM (no item residuals) ====\n")

t_glmer1 <- system.time(
  m1_glmer <- glmer(r2 ~ btype + situ + mode + (1 | id),
                    data = VerbAgg, family = binomial)
)["elapsed"]

t_eirt1 <- system.time(
  m1_eirt <- fit_eirt(resp, item_data,
                      difficulty_formula = ~ btype + situ + mode,
                      model = "Rasch", item_residuals = FALSE)
)["elapsed"]

fe1 <- fixef(m1_glmer)
ga1 <- m1_eirt$regression_coefficients$difficulty
cmp1 <- data.frame(
  glmer_beta   = fe1,
  eirt_gamma   = ga1[names(fe1) |> sub("^\\(Intercept\\)$", "(Intercept)", x = _)],
  neg_gamma    = -ga1
)
print(round(cbind(glmer = fe1, minus_eirt_gamma = -ga1), 4))
cat(sprintf("Person SD : glmer %.4f | eirt %.4f\n",
            sqrt(unname(VarCorr(m1_glmer)$id[1])), m1_eirt$ability_sd))
cat(sprintf("logLik    : glmer %.3f | eirt %.3f\n",
            as.numeric(logLik(m1_glmer)), m1_eirt$logLik))
cat(sprintf("time      : glmer %.2fs | eirt %.2fs\n", t_glmer1, t_eirt1))

# =====================================================================
# Model 2: LLTM + error (random item residuals around the regression)
#   glmer:  r2 ~ btype + situ + mode + (1|id) + (1|item)
#   eirt :  item_residuals = TRUE  (the fit_eirt default)
# =====================================================================
cat("\n==== Model 2: LLTM + error (random item residuals) ====\n")

t_glmer2 <- system.time(
  m2_glmer <- glmer(r2 ~ btype + situ + mode + (1 | id) + (1 | item),
                    data = VerbAgg, family = binomial)
)["elapsed"]

t_eirt2 <- system.time(
  m2_eirt <- fit_eirt(resp, item_data,
                      difficulty_formula = ~ btype + situ + mode,
                      model = "Rasch", item_residuals = TRUE)
)["elapsed"]

fe2 <- fixef(m2_glmer)
ga2 <- m2_eirt$regression_coefficients$difficulty
print(round(cbind(glmer = fe2, minus_eirt_gamma = -ga2), 4))
cat(sprintf("Person SD     : glmer %.4f | eirt %.4f\n",
            sqrt(unname(VarCorr(m2_glmer)$id[1])), m2_eirt$ability_sd))
cat(sprintf("Item resid SD : glmer %.4f | eirt %.4f\n",
            sqrt(unname(VarCorr(m2_glmer)$item[1])),
            m2_eirt$residual_sd$difficulty))
cat(sprintf("logLik        : glmer %.3f | eirt %.3f\n",
            as.numeric(logLik(m2_glmer)), m2_eirt$logLik))
cat(sprintf("time          : glmer %.2fs | eirt %.2fs\n", t_glmer2, t_eirt2))

# Item parameter check: estimated difficulties vs glmer item BLUPs
b_eirt <- m2_eirt$item_parameters$difficulty
# glmer-implied item easiness: fixed part + item BLUP; convert to difficulty
Wd <- model.matrix(~ btype + situ + mode, item_data)
beta_item_glmer <- -(Wd %*% fe2[colnames(Wd)] + ranef(m2_glmer)$item[colnames(resp), 1])
cat(sprintf("\nCor(eirt difficulties, glmer-implied difficulties): %.6f\n",
            cor(b_eirt, beta_item_glmer)))

# Person ability check
th <- m2_eirt$person_abilities
u_id <- ranef(m2_glmer)$id[rownames(resp), 1]
cat(sprintf("Cor(eirt abilities, glmer person BLUPs): %.6f\n", cor(th, u_id)))

cat("\nDone.\n")
