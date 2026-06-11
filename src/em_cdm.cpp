// C++ EM core for cognitive diagnosis models (G-DINA / DINA / DINO).
//
// The whole EM loop runs in C++ on BLAS level-3 kernels: with binary
// data, y*log(p) + (1-y)*log(1-p) = y*(log p - log(1-p)) + log(1-p), so
// the E-step is one dgemm L = Y %*% D (D the per-class logit matrix)
// plus a per-class base term, and the expected success counts are one
// dgemm num = t(Y) %*% W. Under missingness the mask matrix supplies the
// base term and denominators through two further dgemm calls. M-steps
// are closed-form weighted proportions per reduced-profile group,
// isotonically projected over the group lattice (Dykstra) when monotone.
//
// Acceleration is SQUAREM (Varadhan & Roland 2008, scheme S3): each
// cycle takes two EM steps, extrapolates theta0 - 2*alpha*r + alpha^2*v
// on the unconstrained scale [log pi, logit items], projects back into
// the monotone cone, and falls back to the plain second EM iterate if
// the true marginal log-likelihood decreased. This cuts the long
// EM tail on near-flat profile-prevalence ridges (2^A profiles) by an
// order of magnitude relative to per-step extrapolation.
//
// One start per call: R draws the random starting values (keeping RNG
// use in R) and picks the best of n_starts calls.

#include <cmath>
#include <vector>
#include <algorithm>

#define R_NO_REMAP
#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

namespace {

// Weighted isotonic regression over a DAG via Dykstra's algorithm.
// Mirrors .isotonic_poset() in R/lca_em.R (edges are 0-based here).
void isotonic_poset(double *x, const double *w, int G,
                    const int *ea, const int *eb, int E) {
  if (E == 0) return;
  std::vector<double> inc_a(E, 0.0), inc_b(E, 0.0);
  for (int it = 0; it < 10000; it++) {
    double delta = 0.0;
    for (int e = 0; e < E; e++) {
      int a = ea[e], b = eb[e];
      double xa = x[a] + inc_a[e];
      double xb = x[b] + inc_b[e];
      double na, nb;
      if (xa > xb) {
        double m = (w[a] * xa + w[b] * xb) / (w[a] + w[b]);
        na = m; nb = m;
      } else {
        na = xa; nb = xb;
      }
      inc_a[e] = xa - na;
      inc_b[e] = xb - nb;
      delta += std::fabs(na - x[a]) + std::fabs(nb - x[b]);
      x[a] = na; x[b] = nb;
    }
    if (delta < 1e-12) break;
  }
  for (int e = 0; e < E; e++) {           // feasibility to precision
    int a = ea[e], b = eb[e];
    if (x[a] > x[b]) {
      double m = (w[a] * x[a] + w[b] * x[b]) / (w[a] + w[b]);
      x[a] = m; x[b] = m;
    }
  }
}

inline double clamp_prob(double p) {
  return std::min(std::max(p, 1e-6), 1.0 - 1e-6);
}

} // namespace

extern "C" {

// Y_      : N x J real matrix, 0/1 with NA for missing
// w_      : N person weights
// gid_    : J x K integer matrix of 1-based group ids
// ng_     : J integers, groups per item
// edges_  : list over items of E_j x 2 integer matrices (1-based)
// mono_   : scalar integer flag
// pi0_    : K starting prevalences
// items0_ : list over items of G_j starting probabilities
// maxit_, tol_ : EM controls (absolute logLik tolerance)
SEXP C_em_cdm(SEXP Y_, SEXP w_, SEXP gid_, SEXP ng_, SEXP edges_,
              SEXP mono_, SEXP pi0_, SEXP items0_, SEXP maxit_, SEXP tol_) {
  const int N = Rf_nrows(Y_);
  const int J = Rf_ncols(Y_);
  const int K = Rf_ncols(gid_);
  const double *Y = REAL(Y_);
  const double *w = REAL(w_);
  const int *gid = INTEGER(gid_);
  const int *ng = INTEGER(ng_);
  const int monotone = Rf_asInteger(mono_);
  const int max_iter = Rf_asInteger(maxit_);
  const double tol = Rf_asReal(tol_);

  // Response copy with missing zeroed; mask only kept under missingness
  std::vector<double> yb((size_t)N * J);
  std::vector<double> mb;
  bool complete = true;
  for (size_t t = 0; t < (size_t)N * J; t++) {
    if (ISNAN(Y[t])) { yb[t] = 0.0; complete = false; }
    else             { yb[t] = Y[t]; }
  }
  if (!complete) {
    mb.resize((size_t)N * J);
    for (size_t t = 0; t < (size_t)N * J; t++) {
      mb[t] = ISNAN(Y[t]) ? 0.0 : 1.0;
    }
  }

  // Group bookkeeping: offsets into the stacked item-parameter vector
  std::vector<int> goff(J + 1, 0);
  for (int j = 0; j < J; j++) goff[j + 1] = goff[j] + ng[j];
  const int P = goff[J];

  // Edges (0-based, stacked per item)
  std::vector<int> ea, eb, eoff(J + 1, 0);
  for (int j = 0; j < J; j++) {
    SEXP ej = VECTOR_ELT(edges_, j);
    int E = Rf_nrows(ej);
    const int *e = INTEGER(ej);
    for (int r = 0; r < E; r++) {
      ea.push_back(e[r] - 1);
      eb.push_back(e[r + E] - 1);
    }
    eoff[j + 1] = eoff[j] + E;
  }

  // Parameters
  std::vector<double> pi(K), items(P);
  {
    const double *p0 = REAL(pi0_);
    for (int k = 0; k < K; k++) pi[k] = p0[k];
    for (int j = 0; j < J; j++) {
      const double *i0 = REAL(VECTOR_ELT(items0_, j));
      for (int g = 0; g < ng[j]; g++) items[goff[j] + g] = clamp_prob(i0[g]);
    }
  }

  // Work arrays
  std::vector<double> logp(P), log1mp(P);
  std::vector<double> D((size_t)J * K);          // logit(p) per item/class
  std::vector<double> B((size_t)J * K);          // log(1-p), missing case
  std::vector<double> base(K), logpi(K);
  std::vector<double> L((size_t)N * K);          // loglik -> posterior W
  std::vector<double> rowm(N), rowsum(N);
  std::vector<double> Nk(K), numfull((size_t)J * K), denfull;
  std::vector<double> num(P), den(P);
  std::vector<double> th0(K + P), th1(K + P), th2(K + P);
  std::vector<double> rvec(K + P), vvec(K + P);
  std::vector<double> pi2(K), items2(P);
  int Gmax = 1;
  for (int j = 0; j < J; j++) Gmax = std::max(Gmax, ng[j]);
  std::vector<double> unit_w(Gmax, 1.0), dgv(Gmax);
  if (!complete) denfull.resize((size_t)J * K);

  const double one = 1.0, zero = 0.0;

  double loglik = NA_REAL;
  bool converged = false;
  int esteps = 0;

  // E-step into L/W (+ loglik); also fills Nk
  auto estep = [&](void) -> double {
    for (int k = 0; k < K; k++) logpi[k] = std::log(pi[k]);
    for (int p = 0; p < P; p++) {
      logp[p] = std::log(items[p]);
      log1mp[p] = std::log1p(-items[p]);
    }
    for (int k = 0; k < K; k++) {
      double bk = 0.0;
      double *Dk = D.data() + (size_t)k * J;
      double *Bk = complete ? nullptr : B.data() + (size_t)k * J;
      for (int j = 0; j < J; j++) {
        int g = goff[j] - 1 + gid[(size_t)k * J + j];
        Dk[j] = logp[g] - log1mp[g];
        if (complete) bk += log1mp[g];
        else Bk[j] = log1mp[g];
      }
      base[k] = bk + logpi[k];
    }

    // L = Yb %*% D  (N x K), then + base (and + Mb %*% B under missingness)
    F77_CALL(dgemm)("N", "N", &N, &K, &J, &one, yb.data(), &N,
                    D.data(), &J, &zero, L.data(), &N FCONE FCONE);
    if (!complete) {
      F77_CALL(dgemm)("N", "N", &N, &K, &J, &one, mb.data(), &N,
                      B.data(), &J, &one, L.data(), &N FCONE FCONE);
      for (int k = 0; k < K; k++) base[k] = logpi[k];
    }

    std::fill(rowm.begin(), rowm.end(), -INFINITY);
    for (int k = 0; k < K; k++) {
      const double bk = base[k];
      double *Lk = L.data() + (size_t)k * N;
      for (int i = 0; i < N; i++) {
        Lk[i] += bk;
        if (Lk[i] > rowm[i]) rowm[i] = Lk[i];
      }
    }
    std::fill(rowsum.begin(), rowsum.end(), 0.0);
    for (int k = 0; k < K; k++) {
      double *Lk = L.data() + (size_t)k * N;
      for (int i = 0; i < N; i++) {
        Lk[i] = std::exp(Lk[i] - rowm[i]);
        rowsum[i] += Lk[i];
      }
    }
    double ll = 0.0;
    for (int i = 0; i < N; i++) ll += w[i] * (rowm[i] + std::log(rowsum[i]));
    // L -> weighted posterior W
    std::fill(Nk.begin(), Nk.end(), 0.0);
    for (int k = 0; k < K; k++) {
      double *Lk = L.data() + (size_t)k * N;
      double nk = 0.0;
      for (int i = 0; i < N; i++) {
        Lk[i] *= w[i] / rowsum[i];
        nk += Lk[i];
      }
      Nk[k] = nk;
    }
    return ll;
  };

  // M-step from the current E-step quantities (Nk + posterior in L)
  auto mstep = [&](void) {
    double nsum = 0.0;
    for (int k = 0; k < K; k++) nsum += Nk[k];
    double psum = 0.0;
    for (int k = 0; k < K; k++) {
      pi[k] = std::max(Nk[k] / nsum, 1e-12);
      psum += pi[k];
    }
    for (int k = 0; k < K; k++) pi[k] /= psum;

    // num = t(Yb) %*% W (J x K), pooled to groups; same for denominators
    F77_CALL(dgemm)("T", "N", &J, &K, &N, &one, yb.data(), &N,
                    L.data(), &N, &zero, numfull.data(), &J FCONE FCONE);
    if (!complete) {
      F77_CALL(dgemm)("T", "N", &J, &K, &N, &one, mb.data(), &N,
                      L.data(), &N, &zero, denfull.data(), &J FCONE FCONE);
    }
    std::fill(num.begin(), num.end(), 0.0);
    std::fill(den.begin(), den.end(), 0.0);
    for (int k = 0; k < K; k++) {
      for (int j = 0; j < J; j++) {
        int g = goff[j] - 1 + gid[(size_t)k * J + j];
        num[g] += numfull[(size_t)k * J + j];
        den[g] += complete ? Nk[k] : denfull[(size_t)k * J + j];
      }
    }
    for (int j = 0; j < J; j++) {
      double *pj = items.data() + goff[j];
      double *dg = dgv.data();
      for (int g = 0; g < ng[j]; g++) {
        dg[g] = std::max(den[goff[j] + g], 1e-12);
        pj[g] = num[goff[j] + g] / dg[g];
      }
      if (monotone && eoff[j + 1] > eoff[j]) {
        isotonic_poset(pj, dg, ng[j],
                       ea.data() + eoff[j], eb.data() + eoff[j],
                       eoff[j + 1] - eoff[j]);
      }
      for (int g = 0; g < ng[j]; g++) pj[g] = clamp_prob(pj[g]);
    }
  };

  auto flatten = [&](std::vector<double> &th) {
    for (int k = 0; k < K; k++) th[k] = std::log(pi[k]);
    for (int p = 0; p < P; p++) {
      th[K + p] = std::log(items[p] / (1.0 - items[p]));
    }
  };

  // Set (pi, items) from an unconstrained vector, projecting items back
  // into the monotone cone
  auto set_from = [&](const std::vector<double> &th) {
    double s = 0.0;
    for (int k = 0; k < K; k++) {
      pi[k] = std::exp(th[k]);
      s += pi[k];
    }
    for (int k = 0; k < K; k++) pi[k] = std::max(pi[k] / s, 1e-12);
    for (int j = 0; j < J; j++) {
      double *pj = items.data() + goff[j];
      for (int g = 0; g < ng[j]; g++) {
        pj[g] = 1.0 / (1.0 + std::exp(-th[K + goff[j] + g]));
      }
      if (monotone && eoff[j + 1] > eoff[j]) {
        isotonic_poset(pj, unit_w.data(), ng[j],
                       ea.data() + eoff[j], eb.data() + eoff[j],
                       eoff[j + 1] - eoff[j]);
      }
      for (int g = 0; g < ng[j]; g++) pj[g] = clamp_prob(pj[g]);
    }
  };

  // ---- SQUAREM (S3) main loop ----
  double ll_prev = -INFINITY;
  while (esteps < max_iter) {
    double ll0 = estep(); esteps++;
    loglik = ll0;
    if (std::fabs(ll0 - ll_prev) < tol) { converged = true; break; }
    ll_prev = ll0;

    flatten(th0);
    mstep();                                       // theta1 = EM(theta0)
    double ll1 = estep(); esteps++;
    loglik = ll1;
    if (std::fabs(ll1 - ll0) < tol) { converged = true; break; }
    ll_prev = ll1;

    flatten(th1);
    mstep();                                       // theta2 = EM(theta1)
    flatten(th2);
    std::copy(pi.begin(), pi.end(), pi2.begin());
    std::copy(items.begin(), items.end(), items2.begin());

    double rTr = 0.0, vTv = 0.0;
    for (int t = 0; t < K + P; t++) {
      rvec[t] = th1[t] - th0[t];
      vvec[t] = (th2[t] - th1[t]) - rvec[t];
      rTr += rvec[t] * rvec[t];
      vTv += vvec[t] * vvec[t];
    }
    if (vTv > 1e-300) {
      double alpha = -std::sqrt(rTr / vTv);
      if (alpha > -1.0) alpha = -1.0;              // never short of one EM step
      if (alpha < -100.0) alpha = -100.0;
      for (int t = 0; t < K + P; t++) {
        th0[t] = th0[t] - 2.0 * alpha * rvec[t] + alpha * alpha * vvec[t];
      }
      set_from(th0);
      double ll_acc = estep(); esteps++;
      loglik = ll_acc;
      if (ll_acc >= ll1) {
        mstep();                                   // free polish step
        ll_prev = ll_acc;
      } else {
        // Reject: fall back to the plain second EM iterate (monotone)
        std::copy(pi2.begin(), pi2.end(), pi.begin());
        std::copy(items2.begin(), items2.end(), items.begin());
        ll_prev = ll1;
      }
    }
  }

  // ---- Final posterior at the returned parameters ----
  estep();   // L now holds w-scaled posteriors; rescale rows to sum 1
  for (int i = 0; i < N; i++) {
    rowsum[i] = 0.0;
  }
  for (int k = 0; k < K; k++) {
    double *Lk = L.data() + (size_t)k * N;
    for (int i = 0; i < N; i++) rowsum[i] += Lk[i];
  }
  SEXP post_ = PROTECT(Rf_allocMatrix(REALSXP, N, K));
  double *post = REAL(post_);
  for (int k = 0; k < K; k++) {
    double *Lk = L.data() + (size_t)k * N;
    for (int i = 0; i < N; i++) {
      post[(size_t)k * N + i] = Lk[i] / rowsum[i];
    }
  }

  // ---- Return ----
  SEXP pi_ = PROTECT(Rf_allocVector(REALSXP, K));
  for (int k = 0; k < K; k++) REAL(pi_)[k] = pi[k];
  SEXP items_ = PROTECT(Rf_allocVector(VECSXP, J));
  for (int j = 0; j < J; j++) {
    SEXP ij = PROTECT(Rf_allocVector(REALSXP, ng[j]));
    for (int g = 0; g < ng[j]; g++) REAL(ij)[g] = items[goff[j] + g];
    SET_VECTOR_ELT(items_, j, ij);
    UNPROTECT(1);
  }

  const char *names[] = {"pi", "items", "loglik", "posterior",
                         "converged", "iterations", ""};
  SEXP out = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(out, 0, pi_);
  SET_VECTOR_ELT(out, 1, items_);
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(loglik));
  SET_VECTOR_ELT(out, 3, post_);
  SET_VECTOR_ELT(out, 4, Rf_ScalarLogical(converged ? 1 : 0));
  SET_VECTOR_ELT(out, 5, Rf_ScalarInteger(esteps));
  UNPROTECT(4);
  return out;
}

} // extern "C"
