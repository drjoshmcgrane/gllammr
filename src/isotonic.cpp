// Weighted isotonic regression over a partial order (.Call interface).
//
// Minimizes sum(w * (y - x)^2) subject to x[a] <= x[b] for each edge:
// Dykstra's cyclic projection onto the half-spaces, each projection a
// weighted two-point pool. Shared backend for the order-restricted
// latent class M-steps (class posets, item-by-class grids); the chain
// special case is handled by PAVA on the R side.

#include <cmath>
#include <vector>

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

extern "C" {

// y_ : numeric values; w_ : positive weights (same length);
// edges_ : E x 2 integer matrix, 1-based (a precedes b)
SEXP C_isotonic_poset(SEXP y_, SEXP w_, SEXP edges_) {
  const int n = Rf_length(y_);
  const int E = Rf_nrows(edges_);
  const double *y = REAL(y_);
  const double *w = REAL(w_);
  const int *ed = INTEGER(edges_);

  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double *x = REAL(out);
  for (int i = 0; i < n; i++) x[i] = y[i];
  if (E == 0) { UNPROTECT(1); return out; }

  std::vector<double> inc_a(E, 0.0), inc_b(E, 0.0);
  for (int it = 0; it < 10000; it++) {
    double delta = 0.0;
    for (int e = 0; e < E; e++) {
      int a = ed[e] - 1;
      int b = ed[e + E] - 1;
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
  for (int e = 0; e < E; e++) {            // feasibility to precision
    int a = ed[e] - 1;
    int b = ed[e + E] - 1;
    if (x[a] > x[b]) {
      double m = (w[a] * x[a] + w[b] * x[b]) / (w[a] + w[b]);
      x[a] = m; x[b] = m;
    }
  }
  UNPROTECT(1);
  return out;
}

} // extern "C"
