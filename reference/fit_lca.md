# Fit Latent Class Analysis Models

Fit latent class models to categorical data

## Usage

``` r
fit_lca(
  formula,
  data = NULL,
  nclass = 2,
  weights = NULL,
  method = c("em", "tmb"),
  ordering = "none",
  item_ordering = "none",
  structure = c("free", "rasch"),
  start = NULL,
  control = list()
)
```

## Arguments

- formula:

  Formula specifying manifest variables (can be a matrix)

- data:

  Data frame containing variables

- nclass:

  Number of latent classes

- weights:

  Optional vector of case weights (one per observation)

- method:

  Estimation method: "em" (default; closed-form EM with Ramsay
  acceleration, the poLCA algorithm) or "tmb" (direct quasi-Newton on
  the marginal likelihood via TMB)

- ordering:

  Order restriction on the classes. One of:

  - `"none"` (default): unrestricted LCA.

  - `"increasing"`: Croon's (1990) ordered latent class model - a total
    order, where every binary item probability and every gaussian
    indicator mean is nondecreasing from class 1 to class K.

  - A *partial order*: a list of class index pairs (or a two-column
    matrix), each pair `c(a, b)` constraining class `a` \\\preceq\\
    class `b` (item probabilities and gaussian means of `a` no greater
    than those of `b`). Classes not connected by any chain of pairs are
    unconstrained relative to each other. Example - a diamond with two
    incomparable intermediate profiles between a low and a high class:
    `ordering = list(c(1, 2), c(1, 3), c(2, 4), c(3, 4))`.

  The constrained M-step is a weighted isotonic regression over the
  class poset (pool-adjacent-violators on a chain, Dykstra's projection
  algorithm on a general DAG), so estimation remains closed-form EM. A
  total order resolves label switching by construction; a partial order
  resolves it up to the automorphisms of the poset (e.g. the two
  incomparable middle classes of a diamond can swap). Requires
  `method = "em"`; not available with categorical (\> 2 category)
  indicators. Note that likelihood-ratio tests against the unrestricted
  model have a non-standard (chi-bar-square) null distribution, and
  AIC/BIC are reported with the nominal parameter count.

- item_ordering:

  Item-monotonicity restriction (invariant item ordering; Croon 1991):
  "none" (default), "increasing" (success probabilities nondecreasing
  across items in column order, within every class), or a list/matrix of
  item index pairs defining a partial order over items. Combined with
  `ordering` this gives the double monotonicity model. Requires
  all-binary indicators and `method = "em"`.

- structure:

  "free" (default) or "rasch": the located latent class model (latent
  class Rasch; Lindsay, Clogg & Grego 1991), where
  `logit P(x_ic = 1) = theta_c - delta_i` - classes lie on an interval
  scale shared with the items. Classes are reported sorted by location;
  `class_locations` and `item_difficulties` are returned. Implies both
  monotonicities, so the ordering arguments must be left at their
  defaults. Requires all-binary indicators.

- start:

  Optional starting values

- control:

  Control parameters

## Value

An object of class `gllamm_lca`

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate 2-class data
set.seed(123)
n <- 500

# Class 1: high probability of yes
class1_probs <- c(0.8, 0.7, 0.9, 0.75)
# Class 2: low probability of yes
class2_probs <- c(0.2, 0.3, 0.1, 0.25)

# Generate data
true_class <- sample(1:2, n, replace = TRUE, prob = c(0.6, 0.4))
data <- matrix(NA, n, 4)
for (i in 1:n) {
  probs <- if (true_class[i] == 1) class1_probs else class2_probs
  data[i, ] <- rbinom(4, 1, probs)
}
colnames(data) <- paste0("Item", 1:4)

# Fit 2-class model
fit <- fit_lca(data, nclass = 2)
summary(fit)
} # }
```
