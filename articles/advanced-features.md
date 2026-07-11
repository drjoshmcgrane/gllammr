# Advanced Features in gllammr

This vignette tours the specialist corners of the GLLAMM framework, all
reachable through the unified
[`gllamm()`](https://drjoshmcgrane.github.io/gllammr/reference/gllamm.md)
interface: joint mixed-type outcomes, frailty survival, rank-ordered
choices, nonparametric latent distributions, adaptive quadrature, and
cluster-robust standard errors.

``` r

library(gllammr)
#> 
#> Attaching package: 'gllammr'
#> The following object is masked from 'package:stats':
#> 
#>     binomial
set.seed(2026)
```

## Mixed responses: joint outcomes sharing a random effect

A continuous severity score and a binary dropout indicator from the same
clinics, linked by a shared clinic effect. The first argument is the
shared random-effects formula; the outcome formulas live in the family:

``` r

n_clinic <- 40; n_per <- 15
clinic <- factor(rep(1:n_clinic, each = n_per))
u <- rnorm(n_clinic, 0, 0.8)
d <- data.frame(
  clinic = clinic,
  age = rnorm(n_clinic * n_per),
  severity = 1 + 0.3 * rnorm(n_clinic * n_per) + u[clinic] +
    rnorm(n_clinic * n_per, 0, 0.7))
d$dropout <- rbinom(nrow(d), 1, plogis(-1 + 0.8 * u[clinic]))

joint <- gllamm(~ 1 | clinic, data = d,
                family = mixed_response(gaussian = severity ~ age,
                                        binomial = dropout ~ age))
joint$coefficients
#> $gaussian
#> (Intercept)         age 
#>  0.94646256 -0.01898972 
#> 
#> $binomial
#> (Intercept)         age 
#> -1.06730557  0.09986005
```

The loading on the shared effect for the second outcome captures how
strongly clinic-level severity and dropout propensity travel together.

## Parametric frailty survival

Exponential or Weibull survival with a shared log-normal frailty per
cluster, `Surv(time, event)` on the left-hand side. The exponential
frailty model is likelihood-equivalent to a Poisson GLMM with a log-time
offset, which the validation suite uses as an exact cross-check.

``` r

g <- factor(rep(1:40, each = 10))
x <- rnorm(400)
u_s <- rnorm(40, 0, 0.6)
t_true <- rexp(400, rate = exp(-1 + 0.5 * x + u_s[g]))
cens <- rexp(400, 0.15)
ds <- data.frame(time = pmin(t_true, cens),
                 status = as.integer(t_true <= cens), x = x, g = g)

surv <- gllamm(Surv(time, status) ~ x + (1 | g), data = ds,
               family = survival_family("exponential"))
coef(surv)$fixed
#> (Intercept)           x 
#>  -1.0340678   0.5865287
```

## Rank-ordered (exploded) logit

When respondents *rank* alternatives rather than pick one, the exploded
logit decomposes each ranking into successive choices. Random effects
act as taste shifters (intercepts cancel within a ranking). Partial
rankings — unranked alternatives — are handled natively.

``` r

n_cases <- 150; n_alt <- 4
dr <- expand.grid(alt = 1:n_alt, chooser = 1:n_cases)
dr$price <- rnorm(nrow(dr))
util <- -0.9 * dr$price + rlogis(nrow(dr))
dr$rank <- ave(-util, dr$chooser, FUN = rank)

rk <- gllamm(rank ~ price, data = dr, family = ranking(case = ~ chooser))
coef(rk)$fixed
#>      price 
#> -0.6111276
```

## Nonparametric latent distributions (NPML)

If normality of the random effects is in doubt, `integration = npml(k)`
replaces the normal with k estimated mass points and masses (Aitkin
1999) — the GLLAMM framework’s nonparametric option, validated against
the npmlreg package.

``` r

g2 <- factor(rep(1:80, each = 8))
cls <- sample(1:2, 80, TRUE, prob = c(0.6, 0.4))
locs <- c(-1, 1.5)
dn <- data.frame(g = g2, x = rnorm(640))
dn$y <- rbinom(640, 1, plogis(locs[cls[g2]] + 0.5 * dn$x))

np <- gllamm(y ~ x + (1 | g), data = dn, family = binomial(),
             integration = npml(2))
np$locations
#> [1] -1.212472  1.577401
np$masses
#> [1] 0.6054162 0.3945838
```

The two mass points recover the bimodal cluster structure that a normal
random effect would have smoothed over.

## Adaptive quadrature

Laplace integration (the default) is fast and accurate for most designs,
but with small clusters and large random-effect variances it can bias
variance components. `integration = aghq(k)` runs k-node adaptive
Gauss-Hermite quadrature — the Stata gllamm algorithm — matching
`glmer(nAGQ = k)`:

``` r

g3 <- factor(rep(1:100, each = 6))
u3 <- rnorm(100, 0, 2)
da <- data.frame(g = g3, x = rnorm(600))
da$y <- rbinom(600, 1, plogis(-0.5 + 0.8 * da$x + u3[g3]))

fit_lap <- gllamm(y ~ x + (1 | g), data = da, family = binomial())
fit_agh <- gllamm(y ~ x + (1 | g), data = da, family = binomial(),
                  integration = aghq(15))
c(laplace = sqrt(fit_lap$coefficients$random_var[[1]][1, 1]),
  aghq = sqrt(fit_agh$coefficients$random_var[[1]][1, 1]))
#>  laplace     aghq 
#> 1.782626 1.861647
```

The adaptive-quadrature estimate of the random-effect SD is the more
trustworthy one in this small-cluster, large-variance regime.

## Cluster-robust (sandwich) standard errors

Model-based standard errors assume the variance structure is right.
Sandwich standard errors from per-cluster scores stay valid under
misspecification:

``` r

fit_g <- gllamm(severity ~ age + (1 | clinic), data = d)
cbind(model = sqrt(diag(vcov(fit_g)))[1:2],
      sandwich = sqrt(diag(vcov(fit_g, type = "sandwich")))[1:2])
#>                 model   sandwich
#> (Intercept) 0.1362311 0.13587748
#> age         0.0320627 0.02928522
```

## Survey weights at multiple levels

Complex surveys weight clusters and respondents differently;
`weights = list(level1 = ..., level2 = ...)` implements the
pseudo-likelihood approach of Rabe-Hesketh & Skrondal (2006). See
[`vignette("weights")`](https://drjoshmcgrane.github.io/gllammr/articles/weights.md)
for a full treatment.

## References

- Aitkin, M. (1999). A general maximum likelihood analysis of variance
  components in generalized linear models. *Biometrics*, 55, 117–128.
- Rabe-Hesketh, S., & Skrondal, A. (2006). Multilevel modelling of
  complex survey data. *JRSS A*, 169, 805–827.
- Rabe-Hesketh, S., Skrondal, A., & Pickles, A. (2005). Maximum
  likelihood estimation of limited and discrete dependent variable
  models with nested random effects. *Journal of Econometrics*, 128,
  301–323.
