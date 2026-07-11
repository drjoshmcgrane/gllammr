# Performance benchmarks: gllammr vs reference packages on medium-size data.
# Usage: Rscript validation/benchmarks.R
# Writes validation/BENCHMARKS.md. Development artifact (.Rbuildignore'd).

library(gllammr)

time_one <- function(label, expr) {
  gc()
  t <- system.time(force(expr))[["elapsed"]]
  cat(sprintf("  %-28s %8.2fs\n", label, t))
  data.frame(fit = label, seconds = round(t, 2))
}

set.seed(123)
results <- list()

## ---- Gaussian GLMM: n = 10,000, 100 groups, random intercept ----
cat("Gaussian GLMM (n=10k, 100 groups)\n")
g <- 100; n <- 10000
grp <- factor(rep(1:g, length.out = n))
x <- rnorm(n)
u <- rnorm(g, 0, 1)
y <- 1 + 0.5 * x + u[as.integer(grp)] + rnorm(n)
d <- data.frame(y = y, x = x, grp = grp)

results$gauss <- rbind(
  time_one("gllammr gaussian", gllamm(y ~ x + (1 | grp), data = d)),
  if (requireNamespace("lme4", quietly = TRUE))
    time_one("lme4 lmer", lme4::lmer(y ~ x + (1 | grp), data = d, REML = FALSE)),
  if (requireNamespace("glmmTMB", quietly = TRUE))
    time_one("glmmTMB", glmmTMB::glmmTMB(y ~ x + (1 | grp), data = d))
)

## ---- Binomial GLMM: n = 10,000, 100 groups ----
cat("Binomial GLMM (n=10k, 100 groups)\n")
d$yb <- rbinom(n, 1, plogis(0.5 * x + u[as.integer(grp)]))

results$binom <- rbind(
  time_one("gllammr binomial",
           gllamm(yb ~ x + (1 | grp), data = d, family = stats::binomial())),
  if (requireNamespace("lme4", quietly = TRUE))
    time_one("lme4 glmer",
             lme4::glmer(yb ~ x + (1 | grp), data = d, family = stats::binomial())),
  if (requireNamespace("glmmTMB", quietly = TRUE))
    time_one("glmmTMB",
             glmmTMB::glmmTMB(yb ~ x + (1 | grp), data = d, family = stats::binomial()))
)

## ---- Gaussian random slopes: n = 10,000 ----
cat("Gaussian slopes GLMM (n=10k, 100 groups)\n")
U <- matrix(rnorm(g * 2), g) %*% chol(matrix(c(1, .3, .3, .5), 2))
d$ys <- 1 + U[as.integer(grp), 1] + (0.5 + U[as.integer(grp), 2]) * x + rnorm(n)

results$slopes <- rbind(
  time_one("gllammr slopes", gllamm(ys ~ x + (x | grp), data = d)),
  if (requireNamespace("lme4", quietly = TRUE))
    time_one("lme4 lmer slopes",
             lme4::lmer(ys ~ x + (x | grp), data = d, REML = FALSE))
)

## ---- Rasch IRT: 1000 x 40 ----
cat("Rasch IRT (1000 persons x 40 items)\n")
np <- 1000; ni <- 40
theta <- rnorm(np); b <- rnorm(ni)
resp <- matrix(rbinom(np * ni, 1, plogis(outer(theta, b, "-"))), np, ni)

results$rasch <- rbind(
  time_one("gllammr Rasch", fit_irt(resp, model = "Rasch")),
  if (requireNamespace("mirt", quietly = TRUE))
    time_one("mirt Rasch",
             mirt::mirt(as.data.frame(resp), 1, itemtype = "Rasch", verbose = FALSE)),
  if (requireNamespace("TAM", quietly = TRUE))
    time_one("TAM Rasch",
             suppressWarnings(TAM::tam.mml(resp = resp, verbose = FALSE)))
)

## ---- GRM: 1000 x 20, 4 categories ----
cat("GRM (1000 persons x 20 items, 4 categories)\n")
ni2 <- 20
taus <- t(sapply(rnorm(ni2), function(b0) b0 + c(-1, 0, 1)))
respP <- sapply(seq_len(ni2), function(j) {
  cum <- sapply(1:3, function(k) plogis(theta - taus[j, k]))
  1L + rowSums(matrix(runif(np), np, 3) < cum)
})

results$grm <- rbind(
  time_one("gllammr GRM", fit_irt(respP, model = "GRM")),
  if (requireNamespace("mirt", quietly = TRUE))
    time_one("mirt graded",
             mirt::mirt(as.data.frame(respP), 1, itemtype = "graded", verbose = FALSE))
)

## ---- LCA: 1000 x 8, 3 classes ----
cat("LCA (1000 x 8 items, 3 classes)\n")
cls <- sample(1:3, np, replace = TRUE, prob = c(.5, .3, .2))
pmat <- matrix(runif(3 * 8, .1, .9), 3, 8)
respL <- matrix(rbinom(np * 8, 1, pmat[cls, ]), np, 8)

results$lca <- rbind(
  time_one("gllammr LCA", fit_lca(respL, nclass = 3, control = list(n_starts = 3))),
  if (requireNamespace("poLCA", quietly = TRUE)) {
    dl <- as.data.frame(respL + 1)
    f <- stats::as.formula(paste0("cbind(", paste(names(dl), collapse = ","), ") ~ 1"))
    time_one("poLCA", poLCA::poLCA(f, data = dl, nclass = 3, nrep = 3, verbose = FALSE))
  }
)

## ---- Marginal predictions: Monte Carlo integrator (n_sim = 1000/5000) ----
# Population-averaged predict(type = "marginal") integrates over the random
# effects by Monte Carlo. The integrator is vectorized (all draws reduced
# column-wise in a couple of matrix ops), replacing the former per-replicate
# R loop. Timed here on binomial random-intercept and random-slope fits.
cat("Marginal predictions (Monte Carlo integrator)\n")
time_predict <- function(label, fit, n_sim) {
  gc()
  ts <- replicate(5, {
    set.seed(1)
    system.time(predict(fit, type = "marginal", n_sim = n_sim,
                        se.fit = TRUE))[["elapsed"]]
  })
  t <- median(ts)
  cat(sprintf("  %-28s %8.3fs\n", label, t))
  data.frame(fit = label, seconds = round(t, 3))
}

set.seed(42)
gm <- 200; nm <- 5000
grpm <- factor(rep(1:gm, length.out = nm)); xm <- rnorm(nm); um <- rnorm(gm)
ybm <- rbinom(nm, 1, plogis(0.3 + 0.5 * xm + um[as.integer(grpm)]))
fit_marg_ri <- gllamm(ybm ~ xm + (1 | grpm),
                      data = data.frame(ybm = ybm, xm = xm, grpm = grpm),
                      family = stats::binomial())

gs <- 100; ns <- 2000
grps <- factor(rep(1:gs, length.out = ns)); xs <- rnorm(ns)
Us <- matrix(rnorm(gs * 2), gs) %*% chol(matrix(c(1, .3, .3, .5), 2))
etas <- 0.2 + Us[as.integer(grps), 1] + (0.5 + Us[as.integer(grps), 2]) * xs
ybs <- rbinom(ns, 1, plogis(etas))
fit_marg_rs <- gllamm(ybs ~ xs + (xs | grps),
                      data = data.frame(ybs = ybs, xs = xs, grps = grps),
                      family = stats::binomial())

results$marginal <- rbind(
  time_predict("marginal RI n=5k nsim=1k", fit_marg_ri, 1000),
  time_predict("marginal RI n=5k nsim=5k", fit_marg_ri, 5000),
  time_predict("marginal RS n=2k nsim=1k", fit_marg_rs, 1000),
  time_predict("marginal RS n=2k nsim=5k", fit_marg_rs, 5000)
)

## ---- Write report ----
all <- do.call(rbind, results)
lines <- c(
  "# gllammr Performance Benchmarks",
  "",
  paste0("Generated ", format(Sys.Date()), " | R ", getRversion(),
         " | gllammr ", as.character(utils::packageVersion("gllammr"))),
  "",
  "| Section | Fit | Seconds |",
  "|---|---|---|"
)
for (sec in names(results)) {
  r <- results[[sec]]
  for (i in seq_len(nrow(r))) {
    lines <- c(lines, paste0("| ", sec, " | ", r$fit[i], " | ", r$seconds[i], " |"))
  }
}
writeLines(lines, "validation/BENCHMARKS.md")
cat("\nWrote validation/BENCHMARKS.md\n")
