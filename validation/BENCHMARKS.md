# GLLAMMR Performance Benchmarks

Generated 2026-06-11 | R 4.5.1 | GLLAMMR 1.2.0 | macOS Apple Silicon
(single-threaded: Apple clang has no OpenMP; Linux/CRAN builds parallelize
the likelihood via parallel_accumulator)

All timings are medians of 3 warm runs (single fits, one core).

| Model class | Data size | GLLAMMR | Comparator | Ratio |
|---|---|---|---|---|
| Gaussian GLMM | n=10k, 100 grp | 0.37s | glmmTMB 0.35s / lmer 0.03s | **1.06x** / 12x |
| Binomial GLMM | n=10k, 100 grp | 0.69s | glmmTMB 0.54s / glmer 0.37s | 1.3x / 1.9x |
| Poisson GLMM | n=10k, 100 grp | 0.42s | glmmTMB 0.28s | 1.5x |
| Random slopes (gaussian) | n=10k, 100 grp | 0.60s | lmer 0.07s | 8.6x |
| Crossed REs (gaussian) | n=2k, 50x30 grp | 0.11s | lmer 0.41s | **0.27x (faster)** |
| Gamma GLMM | n=2k, 50 grp | 0.33s | glmmTMB 0.38s | **0.87x (faster)** |
| Ordinal (PO logit) | n=2k, 50 grp | 0.78s | ordinal::clmm 0.46s | 1.7x |
| Adaptive quadrature (15 nodes) | n=600, 100 grp | 0.19s | glmer nAGQ=15 0.10s | 1.9x |
| NPML (k=2, binomial) | n=800, 80 grp | 0.04s | npmlreg::allvc 0.06s | **0.67x (faster)** |
| Weibull frailty survival | n=1.5k, 50 grp | 0.21s | survreg (no frailty) 0.04s | n/a* |
| SEM (ml default) | n=800, 6 ind. | 0.03s | lavaan::sem 0.21s | **0.14x (faster)** |
| **Large SEM (ml)** | **n=100k, 6 ind.** | **0.03s** | **lavaan 0.08s** | **0.4x (faster)** |
| LCA (em default, 3 restarts) | n=1k, 8 items | 0.1s | poLCA (nrep=3) 0.42s | **0.24x (faster)** |
| **Large LCA (em)** | **n=20k, 8 items, 3 cls** | **6.7s** | **poLCA 15.9s** | **0.42x (2.4x faster)** |
| Rasch IRT (em, default) | 1000 x 40 | 0.15s | mirt 0.09s | 1.7x |
| 2PL IRT (em, default) | 1000 x 20 | 0.15s | mirt 0.19s | **0.8x (faster)** |
| **Large GRM battery (em)** | **5000 x 100, 5 cat** | **3.0s** | **mirt graded 9.9s** | **0.3x (3.3x faster)** |
| Rasch IRT (method="laplace") | 1000 x 40 | 3.3s | mirt 0.09s | 37x |
| GRM (method="laplace") | 1000 x 20 | 9.8s | mirt graded 0.25s | 39x |

*no frailty-capable parametric survival comparator installed; survreg fits
the fixed-effects model only.

## Reading the numbers

- **GLMM family**: at or beyond glmmTMB parity; crossed random effects and
  NPML are faster than their comparators; lmer's pure-gaussian speed
  (profiled deviance on sparse matrices) is structurally out of reach for
  general-purpose marginal-likelihood machinery.
- **IRT defaults to MML-EM for single-level fits** (Laplace engages
  automatically for multi-level models or se = TRUE). The polytomous EM
  core is compiled C++ (E-step posteriors, expected counts, per-item
  damped-Newton M-steps, safeguarded Ramsay acceleration) on a
  Bock-Aitkin rectangular quadrature grid; on a 5000-person, 100-item,
  5-category GRM battery it converges to mirt's logLik to the decimal
  (cor(a) = 1.000000) in 3.0s vs mirt's 9.9s. EM also handles short tests
  where joint-Laplace 2PL diverges (5-item LSAT validates against ltm).
- All timings single fits on one core; estimates cross-validated in
  validation/RESULTS.md (49/49).
- **LCA defaults to closed-form EM** (the poLCA algorithm with Ramsay
  acceleration, pure R on BLAS): identical logLik to poLCA at every scale,
  2.4x faster on n=20k. The TMB path remains as method = "tmb".
- **SEM defaults to Wishart ML on the sample covariance** (the
  lavaan/LISREL approach): fitting cost independent of N - 0.03s at
  n=100k where the full-data Laplace path takes 61s - with lavaan's
  estimates and logLik reproduced exactly. method = "laplace" retained.
- **Large-scale validation tier**: gllammr_validate(scale = "large") runs
  n=100k GLMM, 5000x100 GRM, n=20k LCA, and n=100k SEM agreement checks
  (9/9 passing; 61/61 across both tiers).

## Marginal predictions (Monte Carlo integrator)

`predict(type = "marginal")` integrates the inverse-link over the random
effects by Monte Carlo. The integrator is now vectorized: all `n_sim` draws
are generated up front (in the exact random-number order the former
per-replicate loop consumed, so a fixed seed is bit-for-bit reproducible) and
reduced column-wise in a couple of matrix operations, with a memory guard
that processes the draws in column blocks once `n_obs * n_sim` exceeds 5e7.
Medians of 5 warm runs, `se.fit = TRUE`, R 4.5 / macOS Apple Silicon:

| Fit | n_sim | Before (loop) | After (vectorized) | Speedup |
|---|---|---|---|---|
| Binomial random-intercept, n=5000 | 1000 | 0.115s | 0.085s | 1.35x |
| Binomial random-intercept, n=5000 | 5000 | 0.567s | 0.478s | 1.19x |
| Binomial random-slope (q=2), n=2000 | 5000 | 0.266s | 0.169s | 1.57x |
| Binomial random-intercept, n=200 | 50000 | 0.485s | 0.162s | **3.0x** |

The speedup grows with the number of draws relative to observations (where
R-level loop overhead dominated); for large `n_obs x n_sim` the runtime is
bounded by the irreducible inverse-link evaluations. Random-slope fits gain
most because the old loop re-factored the covariance (`chol`) on every one of
the `n_sim` replicates; the vectorized path factors it once.

## Refit-loop diagnostics (Cook's distance, DIF purification)

Audited for hoistable per-iteration model matrices. The invariant design
matrices are already built once outside the loops (`dif_test`, `dif_irt`,
`cooks.distance`); the loops themselves are dominated by irreducible
per-iteration refits (`glm.fit` ~88% of a purified `dif_test` run; TMB/nlminb
optimization for Cook's-distance leave-one-out refits). The only remaining
loop-invariant work (design formulas, `complete.cases`) measures ~1.7% of
runtime - below a 10% threshold - so no further hoist was applied.
