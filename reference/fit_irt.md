# Fit Item Response Theory Models

Fit dichotomous (Rasch, 2PL, 3PL) or polytomous (GRM, PCM, GPCM, NRM)
IRT models. Optionally include multi-level random effects for
hierarchical or clustered data.

## Usage

``` r
fit_irt(
  response_matrix,
  model = c("Rasch", "2PL", "3PL", "GRM", "PCM", "GPCM", "NRM"),
  person_data = NULL,
  random = NULL,
  weights = NULL,
  mc_items = NULL,
  method = c("auto", "em", "laplace"),
  quad_points = 61,
  se = TRUE,
  start = NULL,
  control = list()
)
```

## Arguments

- response_matrix:

  Matrix of item responses (persons x items). For dichotomous models:
  coded 0/1. For polytomous models: coded 1, 2, ..., K (ordered
  categories).

- model:

  Type of IRT model: Dichotomous: "Rasch", "2PL", "3PL" Polytomous:
  "GRM" (Graded Response), "PCM" (Partial Credit), "GPCM" (Generalized
  Partial Credit), "NRM" (Nominal Response)

- person_data:

  Optional data frame with person-level variables for multi-level
  models. Must have one row per person (rows correspond to rows in
  response_matrix). Used to specify grouping variables in the random
  effects formula.

- random:

  Optional random effects formula using lme4-style syntax. Examples:
  `~ (1 | class)`, `~ (1 | school/class)`,
  `~ (1 | student) + (1 | time)`. Requires person_data to contain the
  grouping variables.

- weights:

  Optional vector of person-level weights (length = number of persons).
  Under `method = "em"` arbitrary non-negative weights are supported
  (each person's log marginal likelihood is weighted directly). Under
  `method = "laplace"` only integer frequency weights are supported;
  they are implemented by exact replication of weighted persons, so
  results are identical to fitting the duplicated data.

- mc_items:

  For 3PL model only: which items have guessing parameters. Can be: NULL
  (default, all items have guessing), logical vector (length = n_items),
  or integer vector (indices of MC items). Non-MC items use 2PL
  likelihood (no guessing).

- method:

  Estimation method. "auto" (default) uses "laplace" whenever
  multi-level structure (`random`) or standard errors (`se = TRUE`, the
  default) require it, and "em" otherwise – i.e. for single-level fits
  with `se = FALSE`, or whenever non-integer person weights (an EM-only
  feature) are supplied. "em" is Bock-Aitkin marginal maximum likelihood
  with fixed Gauss-Hermite quadrature (the mirt/TAM algorithm; typically
  20-50x faster and evaluates the marginal likelihood more exactly than
  the Laplace approximation). "laplace" is the TMB path, required for
  multi-level models. EM person abilities are EAP scores; Laplace
  abilities are posterior modes.

- quad_points:

  Number of quadrature nodes for method = "em" (default 61)

- se:

  Compute parameter standard errors via TMB::sdreport (default TRUE,
  consistent with the rest of the package). SEs require the Laplace
  path, so under `method = "auto"` the default selects "laplace"; pass
  `se = FALSE` to get the faster EM path for single-level models and
  skip SE computation (which roughly doubles the fitting time for large
  person samples). SEs are not yet available under `method = "em"`: an
  explicit `se = TRUE` is then ignored with a warning.

- start:

  Optional starting values

- control:

  Control parameters for optimization

## Value

An object of class `gllamm_irt`

## Examples

``` r
if (FALSE) { # \dontrun{
# Dichotomous example (Rasch)
set.seed(123)
n_persons <- 500
n_items <- 20
theta <- rnorm(n_persons, 0, 1)
difficulty <- rnorm(n_items, 0, 1)

# Generate binary responses
responses <- matrix(NA, n_persons, n_items)
for (i in 1:n_persons) {
  for (j in 1:n_items) {
    p <- plogis(theta[i] - difficulty[j])
    responses[i, j] <- rbinom(1, 1, p)
  }
}

# Fit Rasch model
fit_rasch <- fit_irt(responses, model = "Rasch")
summary(fit_rasch)

# Polytomous example (GRM)
# Generate 5-category responses
responses_poly <- matrix(NA, n_persons, n_items)
thresholds <- matrix(seq(-2, 2, length.out = 4), n_items, 4, byrow = TRUE)
for (i in 1:n_persons) {
  for (j in 1:n_items) {
    probs <- c(plogis(theta[i] - thresholds[j, 1]),
               diff(plogis(theta[i] - thresholds[j, ])),
               1 - plogis(theta[i] - thresholds[j, 4]))
    responses_poly[i, j] <- sample(1:5, 1, prob = probs)
  }
}

# Fit GRM model
fit_grm <- fit_irt(responses_poly, model = "GRM")
summary(fit_grm)

# 3PL with selective guessing (mixed MC and non-MC items)
# Assessment: 20 items, first 15 are MC, last 5 are open-ended
fit_3pl <- fit_irt(responses, model = "3PL", mc_items = 1:15)
# Only items 1-15 get guessing parameters
# Items 16-20 use 2PL likelihood (no guessing)

# Multi-level IRT: students nested in classes
person_data <- data.frame(
  person_id = 1:n_persons,
  class_id = rep(1:10, each = 50)
)
fit_multilevel <- fit_irt(responses, model = "2PL",
                           person_data = person_data,
                           random = ~ (1 | class_id))
# theta_i = theta_0i + u_class[class[i]]
} # }
```
