#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;



// --------------------------------
// - Helper functions -------------
// --------------------------------

static inline void u1_u2_from_sorted(const std::vector<double>& v, double d,
                                     double& u1, double& u2) {
  u1 = 0.0;
  u2 = 0.0;
  int i = 0, j = static_cast<int>(v.size()) - 1;
  while (i < j) {
    double xj = v[j];
    if (v[i] + xj > d) {
      for (int k = i; k < j; ++k) {
        double vk = v[k];
        u1 += (xj - vk) / (xj + vk);
      }
      u2 += (j - i);
      --j;
    } else {
      ++i;
    }
  }
}

static inline void ensure_sorted_copy(NumericVector x, std::vector<double>& v) {
  v.assign(x.begin(), x.end());
  std::sort(v.begin(), v.end());
}

static inline double choose2(int n) {
  return (static_cast<double>(n) * (n - 1)) / 2.0;
}

static inline bool is_finite(double z) {
  return R_finite(z);
}

static inline double sample_variance(const std::vector<double>& vals) {
  double sum = 0.0, sumsq = 0.0;
  int m = 0;
  for (double z : vals) {
    if (is_finite(z)) {
      sum += z;
      sumsq += z * z;
      ++m;
    }
  }
  if (m <= 1) return NA_REAL;
  double mean = sum / m;
  return (sumsq - m * mean * mean) / (m - 1);
}

static inline double g_scalar_from_sorted(const std::vector<double>& v, double d) {
  double u1 = 0.0, u2 = 0.0;
  u1_u2_from_sorted(v, d, u1, u2);
  return (u2 == 0.0) ? NA_REAL : (u1 / u2);
}

// -------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector g_comp(NumericVector x, double d) {
    // Copy to plain vector to avoid mutating the R object and to speed sort
    std::vector<double> v(x.begin(), x.end());
    std::sort(v.begin(), v.end());

    const int n = static_cast<int>(v.size());
    if (n < 2) {
        return NumericVector::create(NA_REAL, NA_REAL, NA_REAL);
    }

    double u1 = 0.0;
    double u2 = 0.0;
    int i = 0;
    int j = n - 1;

    while (i < j) {
        const double xj = v[j];
        if (v[i] + xj > d) {
            // sum_{k=i}^{j-1} (xj - v[k]) / (xj + v[k])
            for (int k = i; k < j; ++k) {
                const double vk = v[k];
                u1 += (xj - vk) / (xj + vk);
            }
            u2 += (j - i);
            --j;
        } else {
            ++i;
        }
    }

    const double denom = (static_cast<double>(n) * (n - 1)) / 2.0;
    return NumericVector::create(u1 / u2, u1 / denom, u2 / denom);
}

// [[Rcpp::export]]
double g_scalar(NumericVector x, double d) {
    std::vector<double> v(x.begin(), x.end());
    std::sort(v.begin(), v.end());

    const int n = static_cast<int>(v.size());
    if (n < 2) return NA_REAL;

    double u1 = 0.0;
    double u2 = 0.0;
    int i = 0;
    int j = n - 1;

    while (i < j) {
        const double xj = v[j];
        if (v[i] + xj > d) {
            for (int k = i; k < j; ++k) {
                const double vk = v[k];
                u1 += (xj - vk) / (xj + vk);
            }
            u2 += (j - i);
            --j;
        } else {
            ++i;
        }
    }
    return u1 / u2;
}

// [[Rcpp::export]]
SEXP var_gamma(NumericVector x, double d, std::string method = "unbiased", int R = 1000) {
  std::vector<double> v;
  ensure_sorted_copy(x, v);
  int n = v.size();
  if (n < 2) return wrap(NA_REAL);

  double U1 = 0.0, U2 = 0.0;
  u1_u2_from_sorted(v, d, U1, U2);
  double denom = choose2(n);
  double u1 = U1 / denom;
  double u2 = U2 / denom;
  double g = (U2 == 0) ? NA_REAL : (U1 / U2);

  if (method == "unbiased") {
    long double c1aa = 0.0L, c1ab = 0.0L, c1bb = 0.0L;
    long double c2aa = 0.0L, c2ab = 0.0L, c2bb = 0.0L;

    for (int i = 0; i < n; ++i) {
      double xi = v[i];
      long double sum_term = 0.0L, sum_term2 = 0.0L, sum_ind = 0.0L;
      for (int k = 0; k < n; ++k) {
        if (i == k) continue;
        double xk = v[k];
        if (xi + xk > d) {
          double term = std::fabs(xk - xi) / (xi + xk);
          sum_term += term;
          sum_term2 += term * term;
          sum_ind += 1.0L;
        }
      }

      c1aa += sum_term * sum_term;
      c1ab += sum_term * sum_ind;
      c1bb += sum_ind * sum_ind;

      c2aa += sum_term2;
      c2ab += sum_term;
      c2bb += sum_ind;
    }

    long double nld = n;
    long double denom_big = nld * (nld - 1) * (nld - 2) * (nld - 3);
    if (denom_big <= 0) return wrap(NA_REAL);

    long double A = (4.0L * c1aa - 2.0L * c2aa) / denom_big;
    long double B = (4.0L * c1ab - 2.0L * c2ab) / denom_big;
    long double C = (4.0L * c1bb - 2.0L * c2bb) / denom_big;

    long double var1 = A - (4*nld - 6) * (u1 * u1) / ((nld - 2)*(nld - 3));
    long double cov = B - (4*nld - 6) * (u1 * u2) / ((nld - 2)*(nld - 3));
    long double var2 = C - (4*nld - 6) * (u2 * u2) / ((nld - 2)*(nld - 3));

    if (u2 == 0.0 || !R_finite(g)) return wrap(NA_REAL);
    long double v_u = (n / u2) * (var1 - 2.0L * g * cov + g*g * var2);
    return wrap((double)v_u);

  } else if (method == "bootstrap") {
    RNGScope scope;
    std::vector<double> values; values.reserve(R);
    for (int b = 0; b < R; ++b) {
      std::vector<double> xb; xb.reserve(n);
      for (int t = 0; t < n; ++t) {
        int idx = std::floor(unif_rand() * n);
        if (idx == n) idx = n - 1;
        xb.push_back(v[idx]);
      }
      std::sort(xb.begin(), xb.end());
      values.push_back(g_scalar_from_sorted(xb, d));
    }
    double varvals = sample_variance(values);
    if (!R_finite(varvals)) return wrap(NA_REAL);
    return wrap(n * u2 * varvals);

  } else if (method == "jackknife") {
    std::vector<double> values; values.reserve(n);
    for (int i = 0; i < n; ++i) {
      std::vector<double> xj;
      xj.insert(xj.end(), v.begin(), v.begin() + i);
      xj.insert(xj.end(), v.begin() + i + 1, v.end());
      values.push_back(g_scalar_from_sorted(xj, d));
    }
    double varvals = sample_variance(values);
    if (!R_finite(varvals)) return wrap(NA_REAL);
    return wrap((n - 1.0) * (n - 1.0) * u2 * varvals);

  } else {
    stop("Unknown method. Use 'unbiased', 'bootstrap', or 'jackknife'.");
    return wrap(NA_REAL);
  }
}

// [[Rcpp::export]]
NumericMatrix ci_gamma(NumericVector x, NumericVector d_vec,
                           std::string method = "unbiased",
                           int R = 1000, double conf_level = 0.95) {
  std::vector<double> v;
  ensure_sorted_copy(x, v);
  int n = v.size(), m = d_vec.size();
  NumericMatrix out(2, m);
  double alpha = 1.0 - conf_level;
  double z = R::qnorm(1.0 - alpha / 2.0, 0.0, 1.0, 1, 0);

  for (int t = 0; t < m; ++t) {
    double d = d_vec[t];
    double U1 = 0.0, U2 = 0.0;
    u1_u2_from_sorted(v, d, U1, U2);
    double ghat = (U2 == 0.0) ? NA_REAL : (U1 / U2);
    double denom = choose2(n);
    double u2_scaled = U2 / denom;

    SEXP vs_sexp = var_gamma(x, d, method, R);
    double vs = as<double>(vs_sexp);
    if (!R_finite(vs) || !R_finite(ghat) || u2_scaled <= 0.0) {
      out(0, t) = NA_REAL;
      out(1, t) = NA_REAL;
      continue;
    }

    double halfwidth = z * std::sqrt(vs) / std::sqrt(n * u2_scaled);
    out(0, t) = std::max(0.0, ghat - halfwidth);
    out(1, t) = std::min(1.0, ghat + halfwidth);
  }

  return out;
}
