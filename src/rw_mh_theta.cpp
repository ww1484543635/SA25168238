#include <Rcpp.h>
using namespace Rcpp;

inline double log_post_theta(double theta, int k1, int k2, int k3, int k4) {
  if (!(theta > 0.0 && theta < 1.0)) return R_NegInf;
  const double p1 = 0.5 + theta / 4.0;
  const double p2 = (1.0 - theta) / 4.0;
  const double p3 = (1.0 - theta) / 4.0;
  const double p4 = theta / 4.0;
  if (p1 <= 0.0 || p2 <= 0.0 || p3 <= 0.0 || p4 <= 0.0) return R_NegInf;
  return k1*std::log(p1) + k2*std::log(p2) + k3*std::log(p3) + k4*std::log(p4);
}

// [[Rcpp::export]]
List rw_mh_theta_cpp_impl(IntegerVector k, int m = 10000, double w = 0.10, double x0 = 0.50) {
  if (k.size() != 4) stop("k must be length 4.");
  if (m < 2) stop("m must be >= 2.");
  if (!(w > 0.0) || !R_finite(w)) stop("w must be positive and finite.");
  if (!(x0 > 0.0 && x0 < 1.0) || !R_finite(x0)) stop("x0 must be in (0,1) and finite.");

  const int k1 = k[0], k2 = k[1], k3 = k[2], k4 = k[3];

  NumericVector x(m);
  x[0] = x0;

  int accept = 0;
  for (int i = 1; i < m; ++i) {
    const double cur = x[i-1];
    const double y = cur + R::runif(-w, w);

    const double log_old = log_post_theta(cur, k1, k2, k3, k4);
    const double log_new = log_post_theta(y,   k1, k2, k3, k4);

    if (!R_finite(log_new)) {
      x[i] = cur;
    } else {
      const double logr = log_new - log_old;
      const double u = std::log(R::runif(0.0, 1.0));
      if (u <= std::min(0.0, logr)) {
        x[i] = y;
        accept++;
      } else {
        x[i] = cur;
      }
    }
  }

  return List::create(
    _["chain"] = x,
    _["accept_rate"] = (double)accept / (double)(m - 1)
  );
}
