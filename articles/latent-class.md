# Latent Class Analysis with gllammr

## Overview

Latent class analysis (LCA) explains associations among categorical
indicators through membership in a small number of unobserved classes.
In the GLLAMM framework it is the discrete-latent-variable member of the
family: replace the normal latent trait with K mass points and the model
is an LCA.

[`fit_lca()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_lca.md)
(equivalently `gllamm(Y, family = lca(...))`) handles binary,
polytomous, and continuous (gaussian) indicators in any mix, missing
responses under MAR, case weights, and — beyond standard LCA —
**order-restricted** and **partially ordered** class structures.
Estimation is accelerated closed-form EM, validated against poLCA
(identical log-likelihoods on the classic carcinoma data and 20,000-case
simulations, at 2–3x poLCA’s speed).

## A basic two-class model

``` r

library(gllammr)
#> 
#> Attaching package: 'gllammr'
#> The following object is masked from 'package:stats':
#> 
#>     binomial
set.seed(2026)
n <- 800
true_class <- rbinom(n, 1, 0.4) + 1
p_yes <- rbind(c(.80, .70, .85, .75, .65),   # class 1: high endorsers
               c(.20, .25, .15, .30, .25))   # class 2: low endorsers
Y <- sapply(1:5, function(j) rbinom(n, 1, p_yes[true_class, j]))
colnames(Y) <- paste0("sym", 1:5)

fit2 <- fit_lca(Y, nclass = 2)
print(fit2)
#> Latent Class Analysis
#> 
#> Number of classes: 2 
#> Number of observations: 800 
#> Number of items: 5 
#> 
#> Class probabilities:
#> Class1 Class2 
#>  0.592  0.408 
#> 
#> Item response probabilities by class:
#>      Class1 Class2
#> sym1  0.799  0.204
#> sym2  0.693  0.253
#> sym3  0.877  0.121
#> sym4  0.787  0.243
#> sym5  0.653  0.275
#> 
#> Model fit:
#> Log-likelihood: -2498.58 
#> AIC: 5019.16 
#> BIC: 5070.69
```

Class-membership posteriors and modal assignments support downstream
profiling; entropy summarizes how cleanly the classes separate:

``` r

post <- fit2$posterior
entropy <- 1 - sum(-post * log(pmax(post, 1e-12))) /
  (nrow(post) * log(fit2$nclass))
round(c(entropy = entropy,
        accuracy = max(mean(fit2$modal_class == true_class),
                       mean(fit2$modal_class == 3 - true_class))), 3)
#>  entropy accuracy 
#>    0.741    0.917
```

## How many classes?

Fit a sequence and compare information criteria (BIC is the usual
default for LCA):

``` r

fits <- lapply(2:4, function(k) fit_lca(Y, nclass = k,
                                        control = list(n_starts = 5)))
data.frame(classes = 2:4,
           logLik = sapply(fits, `[[`, "logLik"),
           BIC = sapply(fits, `[[`, "BIC"))
#>   classes    logLik      BIC
#> 1       2 -2498.582 5070.694
#> 2       3 -2492.952 5099.543
#> 3       4 -2488.575 5130.895
```

BIC correctly favours the two-class structure that generated the data.

## Ordered latent classes

When classes represent *grades* of a single phenomenon — severity,
proficiency, engagement — an unrestricted LCA leaves their ordering to
chance and label switching. `ordering = "increasing"` fits Croon’s
ordered latent class model: every item probability (and every gaussian
indicator mean) is constrained nondecreasing from class 1 to class K.
The constrained M-step is a weighted isotonic regression, so estimation
stays closed-form.

``` r

set.seed(7)
cls3 <- sample(1:3, 900, TRUE, prob = c(.4, .35, .25))
pmat <- rbind(c(.10, .50, .90), c(.20, .40, .80), c(.15, .55, .85),
              c(.30, .50, .70), c(.05, .45, .90), c(.25, .60, .95))
Y3 <- sapply(1:6, function(j) rbinom(900, 1, pmat[j, cls3]))

ord <- fit_lca(Y3, nclass = 3, ordering = "increasing")
round(ord$item_probs, 3)
#>       Class1 Class2 Class3
#> Item1  0.000  0.393  0.842
#> Item2  0.149  0.343  0.762
#> Item3  0.121  0.380  0.808
#> Item4  0.242  0.436  0.701
#> Item5  0.004  0.280  0.887
#> Item6  0.286  0.384  0.942
```

Every row increases left to right by construction: class labels are
pinned, and the model is the formal statement of “these classes are
ordered”. Comparing its log-likelihood with the unrestricted fit tests
whether the ordering story is tenable (note the likelihood-ratio null is
chi-bar-square, so use the comparison descriptively):

``` r

unr <- fit_lca(Y3, nclass = 3, control = list(n_starts = 5))
c(unrestricted = unr$logLik, ordered = ord$logLik)
#> unrestricted      ordered 
#>    -3325.765    -3325.765
```

A negligible gap (as here, where the truth is ordered) supports the
restriction.

## Partially ordered classes

Sometimes classes are ordered but not *linearly* — two intermediate
profiles may each dominate a low class and be dominated by a high class
while being incomparable with each other (mastered skill A vs mastered
skill B). Pass the comparable pairs:

``` r

set.seed(61)
cls4 <- sample(1:4, 1200, TRUE)
# items 1-3 respond to profile 2; items 4-6 to profile 3
pmat4 <- rbind(
  c(.10, .80, .15, .90), c(.15, .85, .20, .90), c(.10, .75, .10, .85),
  c(.10, .15, .80, .90), c(.15, .20, .85, .90), c(.10, .10, .75, .85))
Y4 <- sapply(1:6, function(j) rbinom(1200, 1, pmat4[j, cls4]))

diamond <- fit_lca(Y4, nclass = 4,
                   ordering = list(c(1, 2), c(1, 3), c(2, 4), c(3, 4)))
round(diamond$item_probs, 3)
#>       Class1 Class2 Class3 Class4
#> Item1  0.127  0.150  0.764  0.925
#> Item2  0.168  0.223  0.879  0.916
#> Item3  0.096  0.145  0.798  0.818
#> Item4  0.117  0.749  0.148  0.894
#> Item5  0.145  0.855  0.241  0.903
#> Item6  0.072  0.746  0.114  0.868

chain <- fit_lca(Y4, nclass = 4, ordering = "increasing",
                 control = list(n_starts = 5))
c(diamond = diamond$logLik, total_order = chain$logLik)
#>     diamond total_order 
#>   -4216.864   -4355.846
```

The diamond fits essentially as well as an unrestricted model, while
*no* total order can accommodate the crossing intermediate profiles —
which is precisely the diagnostic logic: if `"increasing"` costs a lot
but your hypothesized lattice does not, the classes are ordered but not
linearly. Fully Q-matrix-structured class spaces are the next step up:
see
[`vignette("cognitive-diagnosis")`](https://drjoshmcgrane.github.io/gllammr/articles/cognitive-diagnosis.md).

## Categorization, ordering or quantification: comparing latent structures

Torres Irribarra & Diakow propose selecting *the structure of the latent
variable itself* by comparing six progressively constrained models:
unconstrained classes (UN), class-monotone ordered classes (MON),
invariant item ordering (IIO), double monotonicity (DM), located latent
classes (LCR — the latent class Rasch model, where
$`\mathrm{logit}\,\pi_{ic} = \theta_c - \delta_i`$ puts classes on an
interval scale), and the Rasch model (RM). Successive comparisons
decompose the order and scale assumptions: UN vs MON/IIO asks whether
ordering is tenable; the single-monotonicity models vs DM whether
persons and items share one progression; DM vs LCR isolates the
interval-scale assumption; LCR vs RM asks whether a continuum beats
located classes. All six are available here — the item-side constraints
through `item_ordering`, the located classes through
`structure = "rasch"` — and one call runs the whole framework:

``` r

d_loc <- local({
  set.seed(9)
  theta <- c(-1.5, -0.5, 0.5, 1.5)
  delta <- seq(1.4, -1.4, length.out = 8)
  cls <- sample(1:4, 1000, TRUE)
  sapply(1:8, function(j) rbinom(1000, 1, plogis(theta[cls] - delta[j])))
})
cmp <- latent_structure_comparison(d_loc, nclass = 4, n_starts = 3)
print(cmp)
#> Latent structure comparison (Torres Irribarra & Diakow framework)
#> Classes: 4 | IIO/DM item order: 1 < 2 < 3 < 4 < 5 < 6 < 7 < 8 
#> 
#>  model                 structure   logLik n_params     AIC  dAIC     BIC   dBIC
#>     UN               qualitative -4721.16       35 9512.31  0.00 9684.09 108.33
#>    MON                   ordinal -4726.43       35 9522.85 10.54 9694.63 118.87
#>    IIO                   ordinal -4729.50       35 9529.01 16.69 9700.78 125.03
#>     DM                   ordinal -4730.54       35 9531.08 18.76 9702.85 127.09
#>    LCR   quantitative (discrete) -4749.37       14 9526.73 14.42 9595.44  19.69
#>     RM quantitative (continuous) -4756.79        9 9531.58 19.27 9575.75   0.00
#> 
#> Lowest BIC: RM 
#> Read successively: UN vs MON/IIO (is ordering tenable?),
#> single vs double monotonicity (one shared progression?),
#> DM vs LCR (interval scale?), LCR vs RM (continuum vs grain).
```

These data were generated from four *located* classes, and the framework
finds quantitative structure: the qualitative and ordinal models buy
almost no likelihood over LCR despite spending far more parameters. The
located class model itself reports the interval-scale quantities:

``` r

lcr <- attr(cmp, "fits")$LCR
round(lcr$class_locations, 2)
#> Class1 Class2 Class3 Class4 
#>  -7.42  -1.27  -0.21   1.10
round(lcr$item_difficulties, 2)
#> Item1 Item2 Item3 Item4 Item5 Item6 Item7 Item8 
#>  1.41  0.90  0.52  0.29 -0.15 -0.52 -0.96 -1.49
```

(A neat identity worth knowing: the located latent class model is
exactly the Rasch model with a nonparametric, mass-point ability
distribution — the package’s test suite verifies its likelihood against
[`fit_npml()`](https://drjoshmcgrane.github.io/gllammr/reference/fit_npml.md)
on the long-format GLMM.)

## Mixed indicator types and missing data

Indicators are auto-detected: 0/1 columns are binary, integer 1..K
columns categorical, anything else gaussian. Missing responses are
handled by the EM under MAR — each person contributes whatever they
answered:

``` r

Ym <- cbind(Y3[, 1:4],
            score = rnorm(900, mean = c(-1, 0, 1)[cls3], sd = .8))
Ym[sample(length(Ym), 300)] <- NA
mix <- fit_lca(Ym, nclass = 3, ordering = "increasing")
round(mix$gaussian_params$means, 2)
#>       Class1 Class2 Class3
#> score  -1.02  -0.05      1
```

The ordering restriction applies to the gaussian means as well — the
class with the higher symptom probabilities also has the higher mean
score, by construction.

## Practical guidance

- Use several random starts (`control = list(n_starts = ...)`) for
  unrestricted models; ordering reduces the label-switching part of the
  multimodality problem but not the rest.
- Profile classes against external variables using `posterior` weights
  rather than modal assignment when entropy is modest.
- Inequality-restricted models report AIC/BIC with the nominal parameter
  count; treat comparisons against unrestricted models as descriptive.

## References

- Croon, M. (1990). Latent class analysis with ordered latent classes.
  *British Journal of Mathematical and Statistical Psychology*, 43,
  171–192.
- Linzer, D. A., & Lewis, J. B. (2011). poLCA: an R package for
  polytomous variable latent class analysis. *Journal of Statistical
  Software*, 42(10).
- Skrondal, A., & Rabe-Hesketh, S. (2004). *Generalized Latent Variable
  Modeling*. Chapman & Hall/CRC.
