# Category-probability predictions for Laplace-fitted polytomous IRT models.
# Regression tests for the defect where predict(type = "probs") only worked
# for EM fits because the Laplace path never assigned the gllamm_irt_poly
# S3 class.

# Simulate polytomous responses (persons x items) from a given model.
.sim_poly <- function(model, n = 300, J = 5, K = 3, seed = 11) {
  set.seed(seed)
  theta <- rnorm(n)
  th <- lapply(seq_len(J), function(j) sort(rnorm(K - 1, 0, 0.8)))
  a <- runif(J, 0.9, 1.4)
  Y <- matrix(0L, n, J)
  pm <- if (model == "PCM") "GPCM" else model
  for (j in seq_len(J)) {
    pr <- gllammr:::irt_category_probs(pm, theta, th[[j]], a[j])
    cum <- t(apply(pr, 1, cumsum))
    r <- runif(n)
    Y[, j] <- 1L + rowSums(r > cum[, -ncol(cum), drop = FALSE])
  }
  Y
}

.check_probs <- function(pr, n, J, K) {
  expect_type(pr, "list")
  expect_length(pr, J)
  for (j in seq_len(J)) {
    expect_equal(dim(pr[[j]]), c(n, K))
    expect_false(anyNA(pr[[j]]))
    expect_true(all(pr[[j]] >= 0 & pr[[j]] <= 1))
    expect_true(all(abs(rowSums(pr[[j]]) - 1) < 1e-8))
  }
}

test_that("Laplace polytomous fits carry the gllamm_irt_poly class", {
  for (m in c("PCM", "GPCM", "GRM", "NRM")) {
    fit <- fit_irt(.sim_poly(m), model = m, method = "laplace")
    expect_s3_class(fit, "gllamm_irt_poly")
    # gllamm_irt_poly must precede gllamm_irt so predict() dispatches to the
    # polytomous method and NextMethod() still reaches the general one.
    cl <- class(fit)
    expect_lt(match("gllamm_irt_poly", cl), match("gllamm_irt", cl))
  }
})

test_that("predict(type = 'probs') works for each Laplace polytomous model", {
  n <- 300; J <- 5; K <- 3
  for (m in c("PCM", "GPCM", "GRM", "NRM")) {
    fit <- fit_irt(.sim_poly(m, n = n, J = J, K = K), model = m,
                   method = "laplace")
    pr <- predict(fit, type = "probs")
    .check_probs(pr, n, J, K)
    # Expected-score shape matches the EM path (persons x items).
    es <- predict(fit, type = "expected")
    expect_equal(dim(es), c(n, J))
    expect_false(anyNA(es))
    expect_true(all(es >= 1 & es <= K))
  }
})

test_that("EM and Laplace category probabilities agree closely", {
  # Same data, two estimators: predicted category probabilities should be
  # near-identical (both maximize essentially the same marginal likelihood).
  for (m in c("PCM", "GPCM", "GRM")) {
    Y <- .sim_poly(m, n = 400, J = 5, K = 3, seed = 21)
    fl <- fit_irt(Y, model = m, method = "laplace")
    fe <- fit_irt(Y, model = m, method = "em")
    pl <- predict(fl, type = "probs")
    pe <- predict(fe, type = "probs")
    mad <- mean(abs(unlist(pl) - unlist(pe)))
    expect_lt(mad, 0.02)
  }
})

test_that("category-probability prediction is reachable via gllamm()/irt()", {
  Y <- .sim_poly("GPCM", n = 300, J = 5, K = 3, seed = 31)
  fit <- gllamm(Y, family = irt("GPCM"))   # se = TRUE default -> Laplace
  expect_s3_class(fit, "gllamm_irt_poly")
  pr <- predict(fit, type = "probs")
  .check_probs(pr, 300, 5, 3)
})
