# Cognitive Diagnosis Models

## From a continuous trait to attribute profiles

IRT summarizes a person with one number. Cognitive diagnosis models
(CDMs) replace that with a *profile* of binary skills: can this student
subtract fractions, find common denominators, simplify? The latent
classes are the $`2^A`$ attribute profiles, and a **Q-matrix** declares
which skills each item requires. The payoff is diagnostic: instead of a
score, each person gets a posterior probability of mastering each skill.

[`fit_cdm()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_cdm.md)
(equivalently `gllamm(Y, family = cdm(Q, ...))`) implements the
saturated **G-DINA** family with **DINA** and **DINO** as special cases,
estimated by a compiled, accelerated EM and validated against the CDM
package (de la Torre’s fraction-subtraction data: identical item
parameters; simulated G-DINA: log-likelihoods equal to 1e-4).

## Simulating a diagnostic assessment

Three skills, 15 items: single-skill items plus items requiring two or
all three skills (conjunctively — the DINA rule).

``` r

library(gllammr)
#> 
#> Attaching package: 'gllammr'
#> The following object is masked from 'package:stats':
#> 
#>     binomial
set.seed(2026)
A <- 3
Q <- rbind(diag(3), diag(3), diag(3),
           c(1, 1, 0), c(0, 1, 1), c(1, 0, 1),
           c(1, 1, 0), c(0, 1, 1), c(1, 1, 1))
colnames(Q) <- c("subtract", "denominator", "simplify")
n <- 1000
alpha <- matrix(rbinom(n * A, 1, 0.55), n, A)   # true skill profiles
guess <- runif(nrow(Q), 0.05, 0.20)
slip  <- runif(nrow(Q), 0.05, 0.20)
eta <- sapply(seq_len(nrow(Q)), function(j) {
  as.integer(rowSums(alpha[, Q[j, ] == 1, drop = FALSE]) == sum(Q[j, ]))
})
Y <- sapply(seq_len(nrow(Q)), function(j) {
  rbinom(n, 1, ifelse(eta[, j] == 1, 1 - slip[j], guess[j]))
})
```

## Fitting the DINA model

DINA is conjunctive: you answer correctly (beyond guessing) only if you
master *every* required skill. Each item gets a guessing and a slip
parameter.

``` r

fit <- fit_cdm(Y, Q, model = "dina")
print(fit)
#> Cognitive Diagnosis Model (DINA)
#> 
#> Attributes: 3 | Latent profiles: 8  
#> Persons: 1000 | Items: 15 | Monotone: TRUE 
#> 
#> Profile prevalences:
#>   000   100   010   110   001   101   011   111 
#> 0.077 0.115 0.130 0.131 0.114 0.131 0.141 0.162 
#> 
#> Attribute mastery prevalences:
#>    subtract denominator    simplify 
#>       0.538       0.563       0.547 
#> 
#> Item guess/slip:
#>        guess  slip
#> Item1  0.044 0.079
#> Item2  0.168 0.091
#> Item3  0.095 0.141
#> Item4  0.143 0.058
#> Item5  0.072 0.222
#> Item6  0.155 0.207
#> Item7  0.111 0.138
#> Item8  0.207 0.119
#> Item9  0.144 0.085
#> Item10 0.128 0.171
#> Item11 0.133 0.173
#> Item12 0.132 0.176
#> Item13 0.128 0.173
#> Item14 0.136 0.130
#> Item15 0.173 0.085
#> 
#> Log-likelihood: -7637.32 | AIC: 15348.64 | BIC: 15530.22
```

Parameter recovery is direct to check here:

``` r

ghat <- sapply(fit$item_params, function(e) e$guess)
shat <- sapply(fit$item_params, function(e) e$slip)
round(rbind(guess_error = max(abs(ghat - guess)),
            slip_error = max(abs(shat - slip))), 3)
#>              [,1]
#> guess_error 0.038
#> slip_error  0.041
mean(fit$modal_attributes == alpha)   # per-skill classification accuracy
#> [1] 0.9773333
```

## The diagnostic output

The clinically useful quantity is the person-by-skill mastery posterior:

``` r

head(round(fit$attribute_posteriors, 3))
#>      subtract denominator simplify
#> [1,]    0.001       0.158    0.005
#> [2,]    0.000       0.188    0.950
#> [3,]    0.995       0.192    0.045
#> [4,]    1.000       0.985    1.000
#> [5,]    0.000       0.998    0.005
#> [6,]    1.000       0.000    1.000
```

Each row is one student; each entry is P(mastered that skill \| their
responses). Profile prevalences (`fit$profile_probs`) describe the
population.

## G-DINA: when conjunction is too strict

The saturated G-DINA (the default, `model = "gdina"`) frees one response
probability per *combination* of required skills — partial mastery can
partially help. Because DINA is nested in G-DINA, a likelihood-ratio
comparison asks whether the conjunctive restriction holds:

``` r

fit_g <- fit_cdm(Y, Q, model = "gdina")
c(dina = fit$logLik, gdina = fit_g$logLik,
  lr = 2 * (fit_g$logLik - fit$logLik))
#>         dina        gdina           lr 
#> -7637.317987 -7633.619303     7.397368
```

A small statistic relative to the extra parameters (here the data are
truly DINA) supports the simpler model. `model = "dino"` gives the
disjunctive counterpart (any one skill suffices).

Item kernels under G-DINA read directly as P(correct) per reduced skill
pattern — for a two-skill item, four probabilities:

``` r

fit_g$item_params[[10]]
#> $measured
#> [1] "subtract"    "denominator"
#> 
#> $prob
#>        00        01        10        11 
#> 0.1079932 0.1366552 0.1366170 0.8299876
```

## Monotonicity

By default (`monotone = TRUE`) mastering more of an item’s skills can
never *lower* the success probability — enforced by isotonic regression
over the skill lattice in the M-step, the same machinery behind
gllammr’s order-restricted latent class models. Set `monotone = FALSE`
for the unconstrained saturated model; the constrained fit can never
beat it on likelihood, and a large gap flags items whose kernels violate
the ordering.

## Attribute hierarchies

If skills have prerequisites — you cannot simplify fractions before you
can subtract them — the profile space shrinks. Declare prerequisite
pairs and impossible profiles are removed:

``` r

fit_h <- fit_cdm(Y, Q, model = "dina",
                 hierarchy = list(c("subtract", "simplify")))
fit_h$profile_labels    # profiles with simplify-but-not-subtract are gone
#> [1] "000" "100" "010" "110" "101" "111"
```

Linear, convergent, and divergent hierarchies are all expressible as
pair lists; cycles are rejected.

## Connection to the GLLAMM framework

A CDM is a structured latent class model: with one attribute, the
unconstrained G-DINA *is* a two-class LCA (gllammr’s test suite verifies
the log-likelihoods agree to 1e-11). The Q-matrix supplies the equality
structure (items depend only on their measured skills) and monotonicity
the order structure — both restrictions of the GLLAMM parameter space,
not departures from it. See
[`vignette("latent-class")`](https://drjoshmcgrane.github.io/gllammr/articles/latent-class.md)
for the unstructured and order-restricted latent class models this
builds on.

## References

- de la Torre, J. (2011). The generalized DINA model framework.
  *Psychometrika*, 76, 179–199.
- Junker, B. W., & Sijtsma, K. (2001). Cognitive assessment models with
  few assumptions. *Applied Psychological Measurement*, 25, 258–272.
- Rupp, A. A., Templin, J., & Henson, R. A. (2010). *Diagnostic
  Measurement: Theory, Methods, and Applications*. Guilford.
