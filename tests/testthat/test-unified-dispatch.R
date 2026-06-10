# Stage D: full family dispatch through gllamm(), multi-term random effects
# (crossed and nested), and stats-interop of the binomial constructor.

test_that("crossed random effects match lme4 exactly", {
  skip_if_not_installed("lme4")
  set.seed(31)
  n <- 1200
  g1 <- factor(sample(1:25, n, TRUE)); g2 <- factor(sample(1:15, n, TRUE))
  x <- rnorm(n)
  u1 <- rnorm(25, 0, 0.9); u2 <- rnorm(15, 0, 0.6)
  d <- data.frame(x = x, g1 = g1, g2 = g2,
                  y = 1 + 0.5 * x + u1[as.integer(g1)] + u2[as.integer(g2)] + rnorm(n))

  fit <- gllamm(y ~ x + (1 | g1) + (1 | g2), data = d)
  ref <- lme4::lmer(y ~ x + (1 | g1) + (1 | g2), data = d, REML = FALSE)

  expect_equal(unname(coef(fit)$fixed), unname(lme4::fixef(ref)), tolerance = 1e-4)
  expect_equal(fit$logLik, as.numeric(logLik(ref)), tolerance = 1e-4)
  expect_equal(sqrt(fit$coefficients$random_var$g1[1, 1]),
               unname(attr(lme4::VarCorr(ref)$g1, "stddev")), tolerance = 1e-3)
  expect_equal(sqrt(fit$coefficients$random_var$g2[1, 1]),
               unname(attr(lme4::VarCorr(ref)$g2, "stddev")), tolerance = 1e-3)
})

test_that("nested (1 | a/b) expands and matches lme4", {
  skip_if_not_installed("lme4")
  set.seed(32)
  n <- 1500
  a <- factor(rep(1:10, each = 150)); b <- factor(rep(1:50, each = 30))
  x <- rnorm(n)
  ua <- rnorm(10); ub <- rnorm(50, 0, 0.5)
  d <- data.frame(x = x, a = a, b = b,
                  y = 2 + 0.3 * x + ua[as.integer(a)] + ub[as.integer(b)] + rnorm(n))

  fit <- gllamm(y ~ x + (1 | a/b), data = d)
  ref <- lme4::lmer(y ~ x + (1 | a/b), data = d, REML = FALSE)

  expect_equal(fit$logLik, as.numeric(logLik(ref)), tolerance = 1e-4)
  expect_named(fit$coefficients$random_var, c("a", "a:b"))
})

test_that("crossed binomial random effects match lme4", {
  skip_if_not_installed("lme4")
  set.seed(33)
  n <- 1200
  g1 <- factor(sample(1:25, n, TRUE)); g2 <- factor(sample(1:15, n, TRUE))
  x <- rnorm(n)
  u1 <- rnorm(25, 0, 0.9); u2 <- rnorm(15, 0, 0.6)
  d <- data.frame(x = x, g1 = g1, g2 = g2)
  d$yb <- rbinom(n, 1, plogis(0.5 * x + u1[as.integer(g1)] + u2[as.integer(g2)]))

  fit <- gllamm(yb ~ x + (1 | g1) + (1 | g2), data = d, family = stats::binomial())
  ref <- lme4::glmer(yb ~ x + (1 | g1) + (1 | g2), data = d,
                     family = stats::binomial())

  expect_equal(fit$logLik, as.numeric(logLik(ref)), tolerance = 1e-3)
})

test_that("gllamm() dispatches irt(), lca(), multinomial() families", {
  set.seed(34)
  np <- 150; ni <- 8
  resp <- matrix(rbinom(np * ni, 1, plogis(outer(rnorm(np), rnorm(ni), "-"))),
                 np, ni)

  fit_i <- gllamm(resp, family = irt("Rasch"))
  expect_s3_class(fit_i, "gllamm_irt")
  expect_true(fit_i$convergence$converged)

  # Multi-level IRT through the unified interface
  pd <- data.frame(cls = factor(rep(1:10, each = 15)))
  fit_ml <- gllamm(resp, data = pd, family = irt("Rasch"), random = ~ (1 | cls))
  expect_s3_class(fit_ml, "gllamm_irt_multilevel")

  fit_l <- gllamm(matrix(rbinom(200 * 4, 1, 0.5), 200, 4), family = lca(2))
  expect_s3_class(fit_l, "gllamm_lca")

  d <- data.frame(x = rnorm(300), g = factor(rep(1:15, each = 20)))
  u <- rnorm(15)
  p <- cbind(1, exp(0.5 * d$x + u[as.integer(d$g)]), exp(-0.3 * d$x))
  d$ym <- apply(p / rowSums(p), 1, function(pr) sample(1:3, 1, prob = pr))
  fit_m <- gllamm(ym ~ x + (1 | g), data = d, family = multinomial())
  expect_s3_class(fit_m, "gllamm_multinomial")
})

test_that("GLLAMMR binomial() stays usable by glm and lme4", {
  set.seed(35)
  d <- data.frame(x = rnorm(200), g = factor(rep(1:10, each = 20)),
                  yb = rbinom(200, 1, 0.5))

  fam <- binomial()
  expect_s3_class(fam, "binomial_family")
  expect_true(is.function(fam$linkinv))   # full stats family object

  g <- glm(yb ~ x, data = d, family = binomial(link = "probit"))
  expect_s3_class(g, "glm")

  skip_if_not_installed("lme4")
  m <- suppressMessages(lme4::glmer(yb ~ x + (1 | g), data = d, family = binomial()))
  expect_s4_class(m, "glmerMod")
})

test_that("matrix response with wrong family errors clearly", {
  expect_error(gllamm(matrix(0, 5, 2), family = stats::gaussian()),
               "data")
  d <- data.frame(y = rnorm(10), x = rnorm(10))
  expect_error(gllamm(y ~ x, data = d, family = irt("Rasch")),
               "response matrix")
})
