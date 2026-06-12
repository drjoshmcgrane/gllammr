# Package-wide audit probe: fit every model class including its newest
# variants, call every applicable S3 method, and check basic invariants.
# Findings feed the audit fix list; this script is a dev artifact.

suppressMessages(library(gllammr))
set.seed(2026)

probe <- function(label, fit, methods = c("print", "summary", "fitted",
                                          "residuals", "predict",
                                          "simulate", "icc", "VarCorr",
                                          "ranef", "fixef", "vcov",
                                          "coef", "logLik")) {
  for (m in methods) {
    res <- tryCatch({
      out <- suppressMessages(switch(m,
        print = { capture.output(print(fit)); "ok" },
        summary = { capture.output(summary(fit)); "ok" },
        simulate = { s <- simulate(fit, nsim = 2, seed = 1)
                     sprintf("ok [%dx%d]", NROW(s), NCOL(s)) },
        predict = { p <- predict(fit)
                    sprintf("ok [%s]", paste(class(p)[1], length(p))) },
        { do.call(m, list(fit)); "ok" }))
      out
    }, error = function(e) paste("ERROR:", trimws(strsplit(
         conditionMessage(e), "\n")[[1]][1])),
       warning = function(w) paste("WARN:", trimws(conditionMessage(w))))
    if (!startsWith(res, "ok")) {
      cat(sprintf("%-28s %-10s %s\n", label, m, res))
    }
  }
  invisible(NULL)
}

n <- 600
g1 <- factor(rep(1:30, each = 20))
g2 <- factor(rep(1:20, 30))
x <- rnorm(n)
u1 <- rnorm(30, 0, 0.8); u2 <- rnorm(20, 0, 0.5)

cat("=== GLMM variants ===\n")
d <- data.frame(g1 = g1, g2 = g2, x = x,
                yg = 1 + 0.5 * x + u1[g1] + rnorm(n),
                yb = rbinom(n, 1, plogis(0.5 * x + u1[g1])),
                yp = rpois(n, exp(0.3 * x + 0.5 * u1[g1])))
d$yg2 <- d$yg + u2[g2]
probe("gaussian single", gllamm(yg ~ x + (1 | g1), data = d))
probe("gaussian crossed", gllamm(yg2 ~ x + (1 | g1) + (1 | g2), data = d))
probe("gaussian slopes", gllamm(yg ~ x + (x | g1), data = d))
probe("binomial crossed", gllamm(yb ~ x + (1 | g1) + (1 | g2), data = d,
                                 family = binomial()))
probe("poisson single", gllamm(yp ~ x + (1 | g1), data = d,
                               family = poisson()))
probe("binomial aghq", gllamm(yb ~ x + (1 | g1), data = d,
                              family = binomial(), integration = aghq(7)))
probe("binomial npml", gllamm(yb ~ x + (1 | g1), data = d,
                              family = binomial(), integration = npml(2)))

cat("=== multinomial ===\n")
e2 <- exp(0.4 + 0.6 * x + u1[g1]); e3 <- exp(-0.2 + 1.0 * x + u1[g1])
den <- 1 + e2 + e3; r <- runif(n)
d$ym <- factor(ifelse(r < 1 / den, "a", ifelse(r < (1 + e2) / den, "b", "c")))
probe("multinomial single", fit_multinomial(ym ~ x + (1 | g1), data = d))
probe("multinomial crossed", fit_multinomial(ym ~ x + (1 | g1) + (1 | g2),
                                             data = d))

cat("=== survival ===\n")
tt <- rexp(n, exp(-1 + 0.5 * x + u1[g1])); cc <- rexp(n, 0.2)
ds <- data.frame(time = pmin(tt, cc), status = as.integer(tt <= cc),
                 x = x, g = g1)
probe("survival exponential", fit_survival(Surv(time, status) ~ x + (1 | g),
                                           data = ds,
                                           distribution = "exponential"))
probe("survival weibull", fit_survival(Surv(time, status) ~ x + (1 | g),
                                       data = ds,
                                       distribution = "weibull"))

cat("=== rank ===\n")
dr <- expand.grid(alt = 1:4, chooser = 1:120)
dr$price <- rnorm(nrow(dr))
util <- -0.8 * dr$price + rlogis(nrow(dr))
dr$rank <- ave(-util, dr$chooser, FUN = rank)
probe("rank", fit_rank(rank ~ price, case = ~ chooser, data = dr))

cat("=== mixed responses ===\n")
dm <- data.frame(g = g1, x = x,
                 y1 = 1 + 0.3 * x + u1[g1] + rnorm(n, 0, 0.7),
                 y2 = rbinom(n, 1, plogis(0.8 * u1[g1])))
probe("mixed responses", fit_mixed(list(gaussian = y1 ~ x,
                                        binomial = y2 ~ x),
                                   random = ~ 1 | g, data = dm))

cat("=== IRT variants ===\n")
np <- 500; ni <- 10
theta <- rnorm(np)
b <- seq(-1.5, 1.5, length.out = ni)
resp <- sapply(1:ni, function(j) rbinom(np, 1, plogis(theta - b[j])))
resp3 <- sapply(1:8, function(j) {
  1L + rowSums(outer(theta - (j - 4) / 2 + rlogis(np), c(-0.8, 0.8), ">"))
})
probe("irt rasch EM", fit_irt(resp, model = "Rasch"))
probe("irt rasch laplace", fit_irt(resp, model = "Rasch",
                                   method = "laplace"))
probe("irt 2pl EM", fit_irt(resp, model = "2PL"))
probe("irt GRM EM", fit_irt(resp3, model = "GRM"))
pd <- data.frame(sch = factor(rep(1:25, each = 20)))
probe("irt multilevel", fit_irt(resp, model = "Rasch", person_data = pd,
                                random = ~ (1 | sch)))

cat("=== EIRT variants ===\n")
idata <- data.frame(z = rnorm(ni))
probe("eirt rasch", fit_eirt(resp, idata, difficulty_formula = ~ z,
                             model = "Rasch"))
idata3 <- data.frame(z = rnorm(8))
probe("eirt pcm", fit_eirt(resp3, idata3, difficulty_formula = ~ z,
                           model = "PCM"))
probe("eirt lpcm", fit_eirt(resp3, idata3, difficulty_formula = ~ 1,
                            threshold_formula = ~ z, model = "PCM"))
probe("eirt multilevel", fit_eirt(resp, idata, difficulty_formula = ~ z,
                                  model = "Rasch", person_data = pd,
                                  random = ~ (1 | sch)))
probe("eirt weighted", fit_eirt(resp, idata, difficulty_formula = ~ z,
                                model = "Rasch",
                                weights = rep(1:2, length.out = np)))

cat("=== LCA / CDM / SEM ===\n")
cls <- rbinom(np, 1, 0.45) + 1
Yl <- sapply(1:5, function(j) rbinom(np, 1, c(0.2, 0.8)[cls]))
probe("lca", fit_lca(Yl, nclass = 2))
probe("lca ordered", fit_lca(Yl, nclass = 2, ordering = "increasing"))
probe("lca rasch struct", fit_lca(resp[, 1:6], nclass = 3,
                                  structure = "rasch"))
Q <- rbind(diag(2), diag(2), c(1, 1))
alpha <- matrix(rbinom(np * 2, 1, 0.5), np, 2)
Yc <- sapply(1:5, function(j) {
  m <- which(Q[j, ] == 1)
  rbinom(np, 1, 0.15 + 0.7 * (rowSums(alpha[, m, drop = FALSE]) ==
                                length(m)))
})
probe("cdm dina", fit_cdm(Yc, Q, model = "dina"))
f1 <- rnorm(np)
dsem <- data.frame(x1 = f1 + rnorm(np, 0, .6),
                   x2 = 0.8 * f1 + rnorm(np, 0, .6),
                   x3 = 1.2 * f1 + rnorm(np, 0, .6))
probe("sem ml", fit_sem(measurement = list(f = ~ x1 + x2 + x3),
                        data = dsem))
dsem$x1[1:50] <- NA
probe("sem fiml", fit_sem(measurement = list(f = ~ x1 + x2 + x3),
                          data = dsem, missing = "fiml"))

cat("=== done ===\n")
