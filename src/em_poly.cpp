// C++ EM core for polytomous IRT (GRM / PCM / GPCM / NRM).
//
// Bock-Aitkin MML-EM with fixed Gauss-Hermite quadrature. The whole EM
// loop runs in C++: the E-step accumulates person-by-node posteriors and
// expected category counts; M-steps take a few damped Newton steps per
// item with finite-difference derivatives (any improvement preserves the
// generalized-EM convergence guarantee; the outer loop monitors the true
// marginal log-likelihood). Parameterizations mirror R/irt_em.R exactly:
//   GRM : [tau_1, log-spacings (K-2), log a]
//   PCM : [delta_1..delta_{K-1}]            (common slope, updated globally)
//   GPCM: [delta_1..delta_{K-1}, log a]
//   NRM : [c_1..c_{K-1}, log a]
//
// PCM runs on the internal N(0,1) scale as GPCM with a shared slope (the
// free-ability-SD equivalence); the R wrapper back-transforms.

#include <cmath>
#include <vector>
#include <algorithm>

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

namespace {

enum Model { GRM = 1, GPCM = 3, NRM = 4 };   // PCM arrives as GPCM+common_a

// Category probabilities for one item at one node. K <= 32.
inline void cat_probs(int model, double a, const double *th, int K,
                      double z, double *out) {
  if (model == GRM) {
    // thresholds th are the ordered taus (already constructed)
    double prev = 1.0;
    for (int k = 0; k < K - 1; k++) {
      double cum = 1.0 / (1.0 + std::exp(-a * (z - th[k])));
      out[k] = prev - cum;
      prev = cum;
    }
    out[K - 1] = prev;
    for (int k = 0; k < K; k++) out[k] = std::max(out[k], 1e-12);
  } else if (model == GPCM) {
    double eta[32];
    eta[0] = 0.0;
    double mx = 0.0;
    for (int m = 1; m < K; m++) {
      eta[m] = eta[m - 1] + a * (z - th[m - 1]);
      if (eta[m] > mx) mx = eta[m];
    }
    double denom = 0.0;
    for (int m = 0; m < K; m++) {
      out[m] = std::exp(eta[m] - mx);
      denom += out[m];
    }
    for (int m = 0; m < K; m++) out[m] = std::max(out[m] / denom, 1e-12);
  } else {                                   // NRM
    double eta[32];
    eta[0] = 0.0;
    double mx = 0.0;
    for (int k = 1; k < K; k++) {
      eta[k] = a * z + th[k - 1];
      if (eta[k] > mx) mx = eta[k];
    }
    double denom = 0.0;
    for (int k = 0; k < K; k++) {
      out[k] = std::exp(eta[k] - mx);
      denom += out[k];
    }
    for (int k = 0; k < K; k++) out[k] = std::max(out[k] / denom, 1e-12);
  }
}

// Map an item's unconstrained parameter vector to (a, taus/deltas)
inline void unpack(int model, bool common_a, double a_common,
                   const double *par, int K,
                   double &a, double *th) {
  int p_th = K - 1;
  if (model == GRM) {
    th[0] = par[0];
    for (int k = 1; k < p_th; k++) th[k] = th[k - 1] + std::exp(par[k]);
  } else {
    for (int k = 0; k < p_th; k++) th[k] = par[k];
  }
  if (common_a) {
    a = a_common;
  } else {
    a = std::min(std::exp(par[p_th]), 10.0);
  }
}

// Expected complete-data negative log-likelihood for one item:
// -(sum_q sum_k R[q,k] log P_k(z_q))
inline double item_nll(int model, bool common_a, double a_common,
                       const double *par, int K,
                       const double *nodes, int Q, const double *Rj) {
  double a, th[32], P[32];
  unpack(model, common_a, a_common, par, K, a, th);
  double total = 0.0;
  for (int q = 0; q < Q; q++) {
    cat_probs(model, a, th, K, nodes[q], P);
    for (int k = 0; k < K; k++) {
      double r = Rj[q + Q * k];
      if (r > 0) total -= r * std::log(P[k]);
    }
  }
  return total;
}

// A few damped Newton steps with finite-difference derivatives (GEM step)
void mstep_item(int model, bool common_a, double a_common,
                double *par, int n_par, int K,
                const double *nodes, int Q, const double *Rj) {
  const double h = 1e-4;
  std::vector<double> g(n_par), trial(n_par);
  std::vector<double> H(n_par * n_par);

  for (int newton = 0; newton < 4; newton++) {
    double f0 = item_nll(model, common_a, a_common, par, K, nodes, Q, Rj);

    // Gradient and (diagonal-dominant) Hessian by central differences
    for (int i = 0; i < n_par; i++) {
      std::copy(par, par + n_par, trial.begin());
      trial[i] = par[i] + h;
      double fp = item_nll(model, common_a, a_common, trial.data(), K, nodes, Q, Rj);
      trial[i] = par[i] - h;
      double fm = item_nll(model, common_a, a_common, trial.data(), K, nodes, Q, Rj);
      g[i] = (fp - fm) / (2 * h);
      H[i + n_par * i] = std::max((fp - 2 * f0 + fm) / (h * h), 1e-4);
    }
    // Off-diagonals (forward differences of the gradient) - only for small p
    if (n_par <= 6) {
      for (int i = 0; i < n_par; i++) {
        for (int j = i + 1; j < n_par; j++) {
          std::copy(par, par + n_par, trial.begin());
          trial[i] += h; trial[j] += h;
          double fpp = item_nll(model, common_a, a_common, trial.data(), K, nodes, Q, Rj);
          trial[j] = par[j] - h;
          double fpm = item_nll(model, common_a, a_common, trial.data(), K, nodes, Q, Rj);
          trial[i] = par[i] - h; trial[j] = par[j] + h;
          double fmp = item_nll(model, common_a, a_common, trial.data(), K, nodes, Q, Rj);
          trial[j] = par[j] - h;
          double fmm = item_nll(model, common_a, a_common, trial.data(), K, nodes, Q, Rj);
          double hij = (fpp - fpm - fmp + fmm) / (4 * h * h);
          H[i + n_par * j] = H[j + n_par * i] = hij;
        }
      }
    }

    // Solve H step = g (Cholesky with diagonal loading on failure)
    std::vector<double> A(H), step(g);
    double load = 0.0;
    bool ok = false;
    for (int attempt = 0; attempt < 4 && !ok; attempt++) {
      A = H;
      for (int i = 0; i < n_par; i++) A[i + n_par * i] += load;
      ok = true;
      // Cholesky
      for (int i = 0; i < n_par && ok; i++) {
        for (int j = 0; j <= i; j++) {
          double sum = A[i + n_par * j];
          for (int k = 0; k < j; k++) sum -= A[i + n_par * k] * A[j + n_par * k];
          if (i == j) {
            if (sum <= 0) { ok = false; break; }
            A[i + n_par * i] = std::sqrt(sum);
          } else {
            A[i + n_par * j] = sum / A[j + n_par * j];
          }
        }
      }
      if (!ok) load = (load == 0.0) ? 1e-2 : load * 10;
    }
    if (!ok) return;                          // keep current params (GEM-safe)

    // forward/back substitution
    for (int i = 0; i < n_par; i++) {
      double sum = g[i];
      for (int k = 0; k < i; k++) sum -= A[i + n_par * k] * step[k];
      step[i] = sum / A[i + n_par * i];
    }
    for (int i = n_par - 1; i >= 0; i--) {
      double sum = step[i];
      for (int k = i + 1; k < n_par; k++) sum -= A[k + n_par * i] * step[k];
      step[i] = sum / A[i + n_par * i];
    }

    // Damped update: halve until improvement (GEM property)
    double damp = 1.0;
    bool improved = false;
    for (int half = 0; half < 8; half++) {
      for (int i = 0; i < n_par; i++) trial[i] = par[i] - damp * step[i];
      double f1 = item_nll(model, common_a, a_common, trial.data(), K, nodes, Q, Rj);
      if (f1 < f0) {
        std::copy(trial.begin(), trial.end(), par);
        improved = true;
        break;
      }
      damp *= 0.5;
    }
    if (!improved) break;
  }
}

} // namespace

extern "C" {

// Y: N x J integer (1..K, NA), Kvec: J, wts: N, nodes/logA: Q,
// model: 1 GRM, 3 GPCM, 4 NRM; common_a: 1 -> shared slope (PCM/free sigma)
SEXP C_em_poly(SEXP Y_, SEXP Kvec_, SEXP wts_, SEXP nodes_, SEXP logA_,
               SEXP model_, SEXP common_a_, SEXP maxiter_, SEXP tol_) {
  const int *Y = INTEGER(Y_);
  const int *Kvec = INTEGER(Kvec_);
  const double *w = REAL(wts_);
  const double *nodes = REAL(nodes_);
  const double *logA = REAL(logA_);
  const int model = INTEGER(model_)[0];
  const bool common_a = INTEGER(common_a_)[0] == 1;
  const int maxiter = INTEGER(maxiter_)[0];
  const double tol = REAL(tol_)[0];

  const int N = Rf_nrows(Y_);
  const int J = Rf_ncols(Y_);
  const int Q = Rf_length(nodes_);
  int maxK = 2;
  for (int j = 0; j < J; j++) maxK = std::max(maxK, Kvec[j]);
  if (maxK > 31) Rf_error("More than 31 response categories are not supported");

  // Per-item parameter vectors
  std::vector<int> n_par(J);
  std::vector<std::vector<double>> par(J);
  double a_common = 1.0;

  for (int j = 0; j < J; j++) {
    int K = Kvec[j];
    n_par[j] = (K - 1) + (common_a ? 0 : 1);
    par[j].assign(n_par[j], 0.0);

    // Initialize from marginal cumulative proportions
    std::vector<double> counts(K, 0.5);       // smoothing
    double tot = 0.5 * K;
    for (int i = 0; i < N; i++) {
      int y = Y[i + N * j];
      if (y != NA_INTEGER) { counts[y - 1] += 1.0; tot += 1.0; }
    }
    double cum = 0.0;
    std::vector<double> tau0(K - 1);
    for (int k = 0; k < K - 1; k++) {
      cum += counts[k] / tot;
      double p = std::min(std::max(cum, 0.05), 0.95);
      tau0[k] = std::log(p / (1.0 - p));
    }
    if (model == GRM) {
      par[j][0] = tau0[0];
      for (int k = 1; k < K - 1; k++) {
        par[j][k] = std::log(std::max(tau0[k] - tau0[k - 1], 1e-2));
      }
    } else if (model == GPCM) {
      for (int k = 0; k < K - 1; k++) par[j][k] = tau0[k];
    } // NRM: zeros
    if (!common_a) par[j][K - 1] = 0.0;       // log a = 0
  }

  // Work arrays
  std::vector<double> L(N * Q);               // person x node loglik
  std::vector<double> W(N * Q);               // posteriors
  std::vector<double> Pj(Q * maxK);
  std::vector<std::vector<double>> Rj(J);
  for (int j = 0; j < J; j++) Rj[j].assign(Q * Kvec[j], 0.0);

  double loglik = R_NegInf, loglik_old = R_NegInf;
  bool converged = false;
  int iter = 0;

  // Ramsay-style safeguarded acceleration over the flattened parameters
  int n_total = common_a ? 1 : 0;
  for (int j = 0; j < J; j++) n_total += n_par[j];
  auto flatten = [&](std::vector<double> &v) {
    v.resize(n_total);
    int k = 0;
    for (int j = 0; j < J; j++)
      for (int i = 0; i < n_par[j]; i++) v[k++] = par[j][i];
    if (common_a) v[k++] = std::log(a_common);
  };
  auto unflatten = [&](const std::vector<double> &v) {
    int k = 0;
    for (int j = 0; j < J; j++)
      for (int i = 0; i < n_par[j]; i++) par[j][i] = v[k++];
    if (common_a) a_common = std::exp(v[k++]);
  };
  std::vector<double> th_before, th_after, step_prev, th_accel;
  double step_prev_norm = -1.0;
  bool accelerated = false;

  for (iter = 1; iter <= maxiter; iter++) {
    // ---- log-likelihood matrix ----
    std::fill(L.begin(), L.end(), 0.0);
    for (int j = 0; j < J; j++) {
      int K = Kvec[j];
      double a, th[32];
      unpack(model, common_a, a_common, par[j].data(), K, a, th);
      for (int q = 0; q < Q; q++) {
        cat_probs(model, a, th, K, nodes[q], Pj.data() + q * maxK);
      }
      for (int i = 0; i < N; i++) {
        int y = Y[i + N * j];
        if (y == NA_INTEGER) continue;
        for (int q = 0; q < Q; q++) {
          L[i + N * q] += std::log(Pj[q * maxK + (y - 1)]);
        }
      }
    }

    // ---- E-step: posteriors and marginal loglik ----
    loglik = 0.0;
    for (int i = 0; i < N; i++) {
      double mx = R_NegInf;
      for (int q = 0; q < Q; q++) {
        double v = L[i + N * q] + logA[q];
        W[i + N * q] = v;
        if (v > mx) mx = v;
      }
      double s = 0.0;
      for (int q = 0; q < Q; q++) {
        W[i + N * q] = std::exp(W[i + N * q] - mx);
        s += W[i + N * q];
      }
      loglik += w[i] * (mx + std::log(s));
      for (int q = 0; q < Q; q++) W[i + N * q] *= w[i] / s;
    }

    if (accelerated && loglik < loglik_old) {
      // Safeguard: the extrapolated point lost likelihood - revert to the
      // plain EM update and recompute from there next iteration
      unflatten(th_after);
      accelerated = false;
      continue;
    }
    accelerated = false;

    if (std::abs(loglik - loglik_old) < tol) {   // absolute, mirt-style
      converged = true;
      break;
    }
    loglik_old = loglik;

    // ---- expected category counts ----
    for (int j = 0; j < J; j++) {
      int K = Kvec[j];
      std::fill(Rj[j].begin(), Rj[j].end(), 0.0);
      for (int i = 0; i < N; i++) {
        int y = Y[i + N * j];
        if (y == NA_INTEGER) continue;
        double *col = Rj[j].data() + Q * (y - 1);
        for (int q = 0; q < Q; q++) col[q] += W[i + N * q];
      }
    }

    // ---- M-step: per-item Newton ----
    flatten(th_before);
    for (int j = 0; j < J; j++) {
      mstep_item(model, common_a, a_common, par[j].data(), n_par[j],
                 Kvec[j], nodes, Q, Rj[j].data());
    }

    // ---- common slope update (PCM): golden-section on log a ----
    if (common_a) {
      double lo = std::log(0.05), hi = std::log(10.0);
      const double phi = 0.6180339887498949;
      double x1 = hi - phi * (hi - lo), x2 = lo + phi * (hi - lo);
      auto total_nll = [&](double log_a) {
        double a = std::exp(log_a), total = 0.0;
        for (int j = 0; j < J; j++) {
          total += item_nll(model, true, a, par[j].data(), Kvec[j],
                            nodes, Q, Rj[j].data());
        }
        return total;
      };
      double f1 = total_nll(x1), f2 = total_nll(x2);
      for (int it = 0; it < 40 && (hi - lo) > 1e-6; it++) {
        if (f1 < f2) {
          hi = x2; x2 = x1; f2 = f1;
          x1 = hi - phi * (hi - lo); f1 = total_nll(x1);
        } else {
          lo = x1; x1 = x2; f1 = f2;
          x2 = lo + phi * (hi - lo); f2 = total_nll(x2);
        }
      }
      a_common = std::exp((lo + hi) / 2);
    }

    // ---- acceleration: extrapolate along successive EM steps ----
    // step_n = theta_after_mstep - theta_before_mstep; with EM's linear
    // convergence rate r = ||step_n||/||step_{n-1}||, the remaining
    // trajectory sums to ~ step_n * r/(1-r). Safeguarded above.
    flatten(th_after);
    double step_norm = 0.0;
    std::vector<double> step_now(n_total);
    if ((int)th_before.size() == n_total) {
      for (int k = 0; k < n_total; k++) {
        step_now[k] = th_after[k] - th_before[k];
        step_norm += step_now[k] * step_now[k];
      }
      step_norm = std::sqrt(step_norm);
      if (step_prev_norm > 0 && step_norm > 0) {
        double r = step_norm / step_prev_norm;
        if (r > 0.1 && r < 0.98) {
          double gain = std::min(r / (1.0 - r), 20.0);
          th_accel = th_after;
          for (int k = 0; k < n_total; k++) {
            th_accel[k] += gain * step_now[k];
          }
          unflatten(th_accel);
          accelerated = true;
        }
      }
      step_prev_norm = step_norm;
    }
  }
  if (iter > maxiter) iter = maxiter;

  // ---- EAP abilities on the internal scale ----
  SEXP theta_ = PROTECT(Rf_allocVector(REALSXP, N));
  {
    // Recompute posteriors at the final parameters
    std::fill(L.begin(), L.end(), 0.0);
    for (int j = 0; j < J; j++) {
      int K = Kvec[j];
      double a, th[32];
      unpack(model, common_a, a_common, par[j].data(), K, a, th);
      for (int q = 0; q < Q; q++) {
        cat_probs(model, a, th, K, nodes[q], Pj.data() + q * maxK);
      }
      for (int i = 0; i < N; i++) {
        int y = Y[i + N * j];
        if (y == NA_INTEGER) continue;
        for (int q = 0; q < Q; q++) {
          L[i + N * q] += std::log(Pj[q * maxK + (y - 1)]);
        }
      }
    }
    for (int i = 0; i < N; i++) {
      double mx = R_NegInf;
      for (int q = 0; q < Q; q++) {
        double v = L[i + N * q] + logA[q];
        W[i + N * q] = v;
        if (v > mx) mx = v;
      }
      double s = 0.0, m1 = 0.0;
      for (int q = 0; q < Q; q++) {
        double e = std::exp(W[i + N * q] - mx);
        s += e;
        m1 += e * nodes[q];
      }
      REAL(theta_)[i] = m1 / s;
    }
  }

  // ---- assemble result ----
  SEXP thresholds_ = PROTECT(Rf_allocVector(VECSXP, J));
  SEXP a_ = PROTECT(Rf_allocVector(REALSXP, J));
  for (int j = 0; j < J; j++) {
    int K = Kvec[j];
    double a, th[32];
    unpack(model, common_a, a_common, par[j].data(), K, a, th);
    SEXP tj = PROTECT(Rf_allocVector(REALSXP, K - 1));
    for (int k = 0; k < K - 1; k++) REAL(tj)[k] = th[k];
    SET_VECTOR_ELT(thresholds_, j, tj);
    UNPROTECT(1);
    REAL(a_)[j] = a;
  }

  const char *names[] = {"thresholds", "a", "logLik", "iterations",
                         "converged", "theta", ""};
  SEXP out = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(out, 0, thresholds_);
  SET_VECTOR_ELT(out, 1, a_);
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(loglik));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(iter));
  SET_VECTOR_ELT(out, 4, Rf_ScalarLogical(converged));
  SET_VECTOR_ELT(out, 5, theta_);
  UNPROTECT(4);
  return out;
}

} // extern "C"
