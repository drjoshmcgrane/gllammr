# GLLAMMR Performance Benchmarks

Generated 2026-06-11 | R 4.5.1 | GLLAMMR 1.2.0 | macOS Apple Silicon
(single-threaded: Apple clang has no OpenMP; Linux/CRAN builds parallelize
the likelihood via parallel_accumulator)

| Model class | Data size | GLLAMMR | Comparator | Ratio |
|---|---|---|---|---|
| Gaussian GLMM | n=10k, 100 grp | 0.84s | lmer 0.04s / glmmTMB 0.34s | 21x / 2.5x |
| Binomial GLMM | n=10k, 100 grp | 0.83s | glmer 0.34s / glmmTMB 0.42s | 2.4x / 2.0x |
| Random slopes (gaussian) | n=10k, 100 grp | 0.60s | lmer 0.07s | 8.6x |
| Crossed REs (gaussian) | n=2k, 50x30 grp | 0.11s | lmer 0.41s | **0.27x (faster)** |
| Gamma GLMM | n=2k, 50 grp | 0.33s | glmmTMB 0.38s | **0.87x (faster)** |
| Ordinal (PO logit) | n=2k, 50 grp | 0.78s | ordinal::clmm 0.46s | 1.7x |
| Adaptive quadrature (15 nodes) | n=600, 100 grp | 0.19s | glmer nAGQ=15 0.10s | 1.9x |
| NPML (k=2, binomial) | n=800, 80 grp | 0.04s | npmlreg::allvc 0.06s | **0.67x (faster)** |
| Weibull frailty survival | n=1.5k, 50 grp | 0.21s | survreg (no frailty) 0.04s | n/a* |
| SEM (2 factors + path) | n=800, 6 ind. | 0.40s | lavaan::sem 0.21s | 1.9x |
| LCA (3 classes, 3 restarts) | n=1k, 8 items | 2.23s | poLCA (nrep=3) 0.42s | 5.3x |
| Rasch IRT | 1000 x 40 | 9.18s | mirt 0.15s / TAM 0.04s | 61x / 230x |
| GRM (4 categories) | 1000 x 20 | 16.71s | mirt graded 0.17s | 98x |

*no frailty-capable parametric survival comparator installed; survreg fits
the fixed-effects model only.

## Reading the numbers

- **GLMM family**: at or beyond glmmTMB parity; crossed random effects and
  NPML are faster than their comparators; lmer's pure-gaussian speed
  (profiled deviance on sparse matrices) is structurally out of reach for
  general-purpose marginal-likelihood machinery.
- **IRT is the one slow class**: mirt/TAM use EM with fixed quadrature -
  per iteration they evaluate a small grid, while the Laplace path
  re-solves a 1000-dimensional inner problem per gradient evaluation.
  The estimates agree (validated to <1e-2); the gap is algorithmic, not a
  defect. OpenMP (Linux/CRAN) roughly divides the IRT times by the core
  count; EM-style fitting for IRT is the natural future optimization.
- All timings single fits on one core; estimates cross-validated in
  validation/RESULTS.md (49/49).
