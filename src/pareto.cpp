#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;



// -------------------------------
// - Helper functions ------------
// -------------------------------


static inline double choose2(int n) {
  return (static_cast<double>(n) * (n - 1)) / 2.0;
}

static inline void ensure_sorted_copy(NumericVector x, std::vector<double>& v) {
  v.assign(x.begin(), x.end());
  std::sort(v.begin(), v.end());
}


static inline bool finite_dbl(double z) { return R_finite(z); }

// sum over i<j within v_u (already filtered >= u)
static inline double pair_sum_t(const std::vector<double>& v_u) {
  const int m = static_cast<int>(v_u.size());
  double s = 0.0;
  for (int i = 0; i < m - 1; ++i) {
    const double xi = v_u[i];
    for (int j = i + 1; j < m; ++j) {
      const double xj = v_u[j];
      s += (xj - xi) / (xi + xj);
    }
  }
  return s;
}

// compute t_hat (s/choose(n_u,2)), u1 = s/choose(n,2), u2 = choose(n_u,2)/choose(n,2)
static inline void t_components(const std::vector<double>& x_sorted, double u,
                                double& t_hat, double& u1, double& u2) {
  const int n = static_cast<int>(x_sorted.size());
  std::vector<double> v_u; v_u.reserve(n);
  for (int i = 0; i < n; ++i) if (x_sorted[i] >= u) v_u.push_back(x_sorted[i]);
  const int n_u = static_cast<int>(v_u.size());
  const double c2_nu = choose2(n_u);
  const double c2_n  = choose2(n);
  if (n_u < 2 || n < 2 || c2_nu == 0.0 || c2_n == 0.0) {
    t_hat = NA_REAL; u1 = NA_REAL; u2 = NA_REAL; return;
  }
  const double s = pair_sum_t(v_u);
  t_hat = s / c2_nu;
  u1    = s / c2_n;
  u2    = c2_nu / c2_n;
}


// ----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix t_rec(NumericVector x) {
  std::vector<double> v; ensure_sorted_copy(x, v);
  const int n = static_cast<int>(v.size());
  if (n < 2) {
    return NumericMatrix(0, 2);
  }
  NumericVector t(n - 1);
  t[n - 2] = (v[n - 1] - v[n - 2]) / (v[n - 1] + v[n - 2]);

  std::vector<double> xj; xj.reserve(n);
  xj.push_back(v[n - 2]);
  xj.push_back(v[n - 1]);

  for (int k = n - 1; k >= 2; --k) {
    const double xm = v[k - 1];
    double sum_term = 0.0;
    for (double y : xj) sum_term += (y - xm) / (xm + y);

    const double a = static_cast<double>(n - k) / (n - k + 2.0);
    const double b = 2.0 / ((n - k + 1.0) * (n - k));
    const double next_t = a * ( t[k - 1] + b * sum_term );
    t[k - 2] = next_t;

    xj.insert(xj.begin(), xm);
  }

  NumericMatrix out(n - 1, 2);
  for (int i = 0; i < n - 1; ++i) {
    out(i, 0) = v[i];
    out(i, 1) = t[i];
  }
  return out;
}


// [[Rcpp::export]]
double t_scalar(NumericVector x, double u) {
  std::vector<double> v; ensure_sorted_copy(x, v);
  std::vector<double> v_u; v_u.reserve(v.size());
  for (double xi : v) if (xi >= u) v_u.push_back(xi);
  const int n_u = static_cast<int>(v_u.size());
  if (n_u < 2) return NA_REAL;
  const double s = pair_sum_t(v_u);
  return 2.0 * s / (static_cast<double>(n_u) * (n_u - 1));
}

// [[Rcpp::export]]
NumericVector t_comp(NumericVector x, double u) {
  std::vector<double> v; ensure_sorted_copy(x, v);
  double t_hat, u1, u2; t_components(v, u, t_hat, u1, u2);
  return NumericVector::create(t_hat, u1, u2);
}

// [[Rcpp::export]]
SEXP var_pareto(NumericVector x, double u, std::string method, int R = 1000) {
  std::vector<double> v; ensure_sorted_copy(x, v);
  const int n = static_cast<int>(v.size());
  if (n < 2) return wrap(NA_REAL);

  double t_hat = NA_REAL, u1 = NA_REAL, u2 = NA_REAL;
  t_components(v, u, t_hat, u1, u2);
  if (!finite_dbl(u2) || u2 <= 0.0) return wrap(NA_REAL);

  if (method == "unbiased") {
    long double c1aa = 0.0L, c1ab = 0.0L, c1bb = 0.0L;
    long double c2aa = 0.0L, c2ab = 0.0L, c2bb = 0.0L;

    for (int i = 0; i < n; ++i) {
      const double xi = v[i];
      long double sum_term = 0.0L;   // sum( I * term )
      long double sum_term2 = 0.0L;  // sum( I * term^2 )
      long double sum_ind = 0.0L;    // sum( I )

      for (int k = 0; k < n; ++k) if (k != i) {
        const double xk = v[k];
        if (std::min(xi, xk) > u) {
          const double term = std::fabs(xk - xi) / (xi + xk);
          sum_term  += term;
          sum_term2 += term * term;
          sum_ind   += 1.0L;
        }
      }

      c1aa += sum_term * sum_term;
      c1ab += sum_term * sum_ind;
      c1bb += sum_ind   * sum_ind;

      c2aa += sum_term2;
      c2ab += sum_term;
      c2bb += sum_ind;
    }

    const long double nld = static_cast<long double>(n);
    const long double denom_big = nld * (nld - 1) * (nld - 2) * (nld - 3);
    if (denom_big <= 0) return wrap(NA_REAL);

    const long double A = (4.0L*c1aa - 2.0L*c2aa) / denom_big;
    const long double B = (4.0L*c1ab - 2.0L*c2ab) / denom_big;
    const long double C = (4.0L*c1bb - 2.0L*c2bb) / denom_big;

    const long double var1 = A - (4*nld - 6) * (static_cast<long double>(u1)*u1) / ((nld-2)*(nld-3));
    const long double cov  = B - (4*nld - 6) * (static_cast<long double>(u1)*u2) / ((nld-2)*(nld-3));
    const long double var2 = C - (4*nld - 6) * (static_cast<long double>(u2)*u2) / ((nld-2)*(nld-3));

    const long double v_u = (static_cast<long double>(n) / u2) * (var1 - 2.0L*t_hat*cov + (t_hat*t_hat)*var2);
    return wrap(static_cast<double>(std::max(0.0L, v_u)));
  }

  if (method == "bootstrap") {
    RNGScope scope;
    std::vector<double> values; values.reserve(R);
    for (int b = 0; b < R; ++b) {
      std::vector<double> xb; xb.reserve(n);
      for (int t = 0; t < n; ++t) {
        int idx = static_cast<int>(std::floor(unif_rand() * n));
        if (idx == n) idx = n - 1;
        xb.push_back(v[idx]);
      }
      std::sort(xb.begin(), xb.end());
      std::vector<double> v_u; v_u.reserve(n);
      for (double xi : xb) if (xi >= u) v_u.push_back(xi);
      const int n_u = static_cast<int>(v_u.size());
      double ghat_b = NA_REAL;
      if (n_u >= 2) {
        const double s = pair_sum_t(v_u);
        ghat_b = 2.0 * s / (static_cast<double>(n_u) * (n_u - 1));
      }
      if (finite_dbl(ghat_b)) values.push_back(ghat_b);
    }
    if (values.size() <= 1) return wrap(NA_REAL);
    double sum = 0.0, sumsq = 0.0; int m = values.size();
    for (double z : values) { sum += z; sumsq += z*z; }
    double mean = sum / m;
    double varvals = (sumsq - m*mean*mean) / (m - 1);
    return wrap(n * u2 * varvals);
  }

  if (method == "jackknife") {
    std::vector<double> values; values.reserve(n);
    for (int i = 0; i < n; ++i) {
      std::vector<double> xj; xj.reserve(n - 1);
      xj.insert(xj.end(), v.begin(), v.begin() + i);
      xj.insert(xj.end(), v.begin() + i + 1, v.end());
      std::vector<double> v_u; v_u.reserve(n - 1);
      for (double xi : xj) if (xi >= u) v_u.push_back(xi);
      const int n_u = static_cast<int>(v_u.size());
      double ghat_i = NA_REAL;
      if (n_u >= 2) {
        const double s = pair_sum_t(v_u);
        ghat_i = 2.0 * s / (static_cast<double>(n_u) * (n_u - 1));
      }
      if (finite_dbl(ghat_i)) values.push_back(ghat_i);
    }
    if (values.size() <= 1) return wrap(NA_REAL);
    double sum = 0.0, sumsq = 0.0; int m = values.size();
    for (double z : values) { sum += z; sumsq += z*z; }
    double mean = sum / m;
    double varvals = (sumsq - m*mean*mean) / (m - 1);
    return wrap((n - 1.0)*(n - 1.0) * u2 * varvals);
  }

  stop("Unknown method. Use 'unbiased', 'bootstrap', or 'jackknife'.");
  return wrap(NA_REAL);
}


// [[Rcpp::export]]
NumericMatrix ci_pareto(NumericVector x, NumericVector u_vec,
                            std::string method, int R = 1000,
                            double conf_level = 0.95) {
  std::vector<double> v; ensure_sorted_copy(x, v);
  const int n = static_cast<int>(v.size());
  const int m = u_vec.size();
  NumericMatrix out(2, m); 
  const double alpha = 1.0 - conf_level;
  const double z = R::qnorm(1.0 - alpha/2.0, 0.0, 1.0, 1, 0);

  for (int t = 0; t < m; ++t) {
    const double u = u_vec[t];
    double t_hat = NA_REAL, u1 = NA_REAL, u2 = NA_REAL;
    t_components(v, u, t_hat, u1, u2);
    if (!finite_dbl(t_hat) || !finite_dbl(u2) || u2 <= 0.0) {
      out(0, t) = NA_REAL; out(1, t) = NA_REAL; continue;
    }

    SEXP vs_sexp = var_pareto(x, u, method, R);
    double vs = as<double>(vs_sexp);
    if (!finite_dbl(vs)) { out(0, t) = NA_REAL; out(1, t) = NA_REAL; continue; }

    const double halfwidth = z * std::sqrt(vs) / std::sqrt(n * u2);
    const double lower = std::max(0.0, t_hat - halfwidth);
    const double upper = std::min(1.0, t_hat + halfwidth);
    out(0, t) = lower; out(1, t) = upper;
  }
  return out;
}
