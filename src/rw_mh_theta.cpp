#include <Rcpp.h>
using namespace Rcpp;

// 数值稳定的 inverse-logit
static inline double inv_logit(double z) {
  if (z >= 0.0) {
    double ez = std::exp(-z);
    return 1.0 / (1.0 + ez);
  } else {
    double ez = std::exp(z);
    return ez / (1.0 + ez);
  }
}

// 计算 log target(z) 及其梯度 d/dz log target(z)
// target = loglik(theta(z)) + log|dtheta/dz|
// theta = inv_logit(z) ∈ (0,1)
static inline void logp_and_grad_z(double z, int k1, int k2, int k3, int k4,
                                   double &logp, double &grad) {
  const double th = inv_logit(z);
  
  // 极端 z 可能导致 th 变成 0/1（数值上），直接判不合法
  if (!(th > 0.0 && th < 1.0) || !R_finite(th)) {
    logp = R_NegInf;
    grad = 0.0;
    return;
  }
  
  const double p1 = 0.5 + th / 4.0;
  const double p2 = (1.0 - th) / 4.0;
  const double p3 = p2;
  const double p4 = th / 4.0;
  
  if (p1 <= 0.0 || p2 <= 0.0 || p3 <= 0.0 || p4 <= 0.0) {
    logp = R_NegInf;
    grad = 0.0;
    return;
  }
  
  // log-likelihood（去掉多项分布常数项）
  const double loglik =
    (double)k1 * std::log(p1) +
    (double)k2 * std::log(p2) +
    (double)k3 * std::log(p3) +
    (double)k4 * std::log(p4);
  
  // Jacobian: log|dtheta/dz| = log(theta*(1-theta))
  const double logjac = std::log(th) + std::log(1.0 - th);
  
  logp = loglik + logjac;
  // d/dtheta loglik
  const double dloglik_dth =
    ( (double)k1 ) / (4.0 * p1) -
    ( (double)k2 ) / (4.0 * p2) -
    ( (double)k3 ) / (4.0 * p3) +
    ( (double)k4 ) / (4.0 * p4);
  
  // dtheta/dz = th*(1-th)
  const double dth_dz = th * (1.0 - th);
  // d/dz logjac = 1 - 2*theta
  const double dlogjac_dz = 1.0 - 2.0 * th;
  grad = dloglik_dth * dth_dz + dlogjac_dz;
  if (!R_finite(grad)) grad = 0.0;
}

// [[Rcpp::export]]
List rw_mh_theta_cpp_impl(IntegerVector k,
                          int m = 10000,
                          double w = 0.02,   // 现在 w 当作 HMC step size (eps)
                          double x0 = 0.50,
                          int L = 10) {
  RNGScope scope;
  // 参数检查
  if (k.size() != 4) stop("k must be length 4.");
  if (m < 2) stop("m must be >= 2.");
  if (!R_finite(w) || w <= 0.0) stop("w (HMC step size) must be positive and finite.");
  if (!R_finite(x0) || !(x0 > 0.0 && x0 < 1.0)) stop("x0 must be in (0,1) and finite.");
  if (L < 1) stop("L must be >= 1.");
  
  for (int i = 0; i < 4; ++i) {
    if (IntegerVector::is_na(k[i])) stop("k contains NA.");
    if (k[i] < 0) stop("k must be nonnegative.");
  }
  
  const int k1 = k[0], k2 = k[1], k3 = k[2], k4 = k[3];
  
  // 初始 z = logit(x0)
  double z = std::log(x0 / (1.0 - x0));
  
  NumericVector chain(m);
  chain[0] = x0;
  
  int accept = 0;
  
  for (int i = 1; i < m; ++i) {
    
    double logp_cur, grad_cur;
    logp_and_grad_z(z, k1, k2, k3, k4, logp_cur, grad_cur);
    
    // 动量 r ~ N(0,1)
    double r = R::rnorm(0.0, 1.0);
    
    // 保存候选
    double z_prop = z;
    double r_prop = r;
    
    r_prop += 0.5 * w * grad_cur;
    double logp_prop = logp_cur, grad_prop = grad_cur;
    // L 次 leapfrog
    for (int t = 0; t < L; ++t) {
      // 一步位置
      z_prop += w * r_prop;
      // 更新 logp/grad
      logp_and_grad_z(z_prop, k1, k2, k3, k4, logp_prop, grad_prop);
      if (!R_finite(logp_prop)) break;
      // 动量更新（最后一次留给收尾半步）
      if (t != L - 1) r_prop += w * grad_prop;
    }
    // 收尾半步动量
    r_prop += 0.5 * w * grad_prop;
    // 反转动量保证可逆
    r_prop = -r_prop;
    bool take = false;
    if (R_finite(logp_prop)) {
      const double H_cur  = -logp_cur  + 0.5 * r * r;
      const double H_prop = -logp_prop + 0.5 * r_prop * r_prop;
      const double loga = H_cur - H_prop;
      const double u = std::log(R::runif(0.0, 1.0));
      if (u <= std::min(0.0, loga)) take = true;
    }
    
    if (take) {
      z = z_prop;
      accept++;
    }
    chain[i] = inv_logit(z);
  }
  
  return List::create(
    _["chain"] = chain,
    _["accept_rate"] = (double)accept / (double)(m - 1)
  );
}
