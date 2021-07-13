#include <RcppArmadillo.h>

// Threshold function

void threshold(arma::vec& beta, const unsigned& k) {
  unsigned p = beta.size();
  arma::uvec pos_beta = arma::linspace<arma::uvec>(0, p - k - 1, p - k);
  arma::uvec order = arma::sort_index(arma::abs(beta), "ascend");
  beta(order(pos_beta)).fill(0);
}

// Gradient descent function

void gd(const arma::mat& x, const arma::mat& xt, const arma::vec& y, arma::vec& beta,
        arma::vec& eta, double& objval, const int k, const int h, const double step,
        const int max_iter, const double eps) {

  // Preliminary preparation

  int n = x.n_rows;
  double obj = 0.5 * std::pow(arma::norm(y - x * beta - eta, 2), 2);

  // Take gradient descent steps

  for (int m = 1; m <= max_iter; m++) {
    double old_obj = obj;
    beta = beta + step * xt * (y - x * beta - eta);
    threshold(beta, k);
    eta = y - x * beta;
    threshold(eta, n - h);
    obj = 0.5 * std::pow(arma::norm(y - x * beta - eta, 2), 2);
    if (old_obj - obj <= eps) break;
  }

  // Polish coefficients

  arma::uvec nz_j = arma::find(beta != 0);
  arma::uvec z_i = arma::find(eta == 0);
  arma::uvec nz_i = arma::find(eta != 0);
  beta(nz_j) = arma::solve(xt(nz_j, z_i) * x(z_i, nz_j), xt(nz_j, z_i) * y(z_i));
  eta = y - x * beta;
  eta(z_i).fill(0);
  objval = 0.5 * std::pow(arma::norm(y - x * beta - eta, 2), 2);

}

// Neighbourhood search function

// [[Rcpp::export]]
Rcpp::List ns(const arma::mat& x, const arma::vec& y, const arma::vec& k, const arma::vec& h,
              const unsigned max_ns_iter, const unsigned max_gd_iter, const double eps) {

  // Preliminaries

  double step = 1 / std::pow(arma::norm(x), 2);
  unsigned n = x.n_rows;
  unsigned p = x.n_cols;
  int nk = k.size();
  int nh = h.size();
  arma::mat xt = x.t();
  arma::cube beta = arma::cube(p, nk, nh);
  arma::cube eta = arma::cube(n, nk, nh);
  arma::mat objval = arma::mat(nk, nh);

  if (max_ns_iter > 0) {

    //  Initialise coefficient vectors

    for (int j = 0; j < nh; j++) {
      for (int i = 0; i < nk; i++) {
        arma::vec beta_i_j = arma::vec(p);
        arma::vec eta_i_j = arma::vec(n);
        double objval_i_j;
        gd(x, xt, y, beta_i_j, eta_i_j, objval_i_j, k(i), h(j), step, max_gd_iter, eps);
        beta.subcube(0, i, j, p - 1, i, j) = beta_i_j;
        eta.subcube(0, i, j, n - 1, i, j) = eta_i_j;
        objval.submat(i, j, i, j) = objval_i_j;
      }
    }

    // Perform neighbourhood swaps

    for (unsigned iter = 0; iter < max_ns_iter; iter++) {
      double objval_old_sum = arma::accu(objval);

      // Loop over h_1,h_2,...

      for (int j = 0; j < nh; j++) {

        // Loop over k_1,k_2,...

        for (int i = 0; i < nk; i++) {

          // Loop over neighbours of k_i

          for (int a = i - 1; a <= i + 1; a++) {
            if ((a < 0) | (a > nk - 1)) continue;

            // Loop over neighbours of h_j

            for (int b = j - 1; b <= j + 1; b++) {
              if ((std::abs<int>(i - a) + std::abs<int>(j - b) != 1) | (b < 0) | (b > nh - 1)) {
                continue;
              }

              // Neighbour of beta(k_i,h_j) and eta(k_i,h_j)

              arma::vec beta_a_b = beta.subcube(0, a, b, p - 1, a, b);
              arma::vec eta_a_b = eta.subcube(0, a, b, n - 1, a , b);

              // Threshold neighbour so that it is feasible

              threshold(beta_a_b, k(i));
              threshold(eta_a_b, h(j));

              // Run gradient descent and update solution if improved

              double objval_i_j;
              gd(x, xt, y, beta_a_b, eta_a_b, objval_i_j, k(i), h(j), step, max_gd_iter, eps);
              if (objval_i_j < arma::as_scalar(objval.submat(i, j, i, j))) {
                beta.subcube(0, i, j, p - 1, i, j) = beta_a_b;
                eta.subcube(0, i, j, n - 1, i, j) = eta_a_b;
                objval.submat(i, j, i, j) = objval_i_j;
              }

            }

          }

        }

      }

      // Exit if converged

      double objval_sum = arma::accu(objval);
      if (objval_old_sum - objval_sum <= eps) break;

    }

  }

  // Return fitted coefficients

  return Rcpp::List::create(Rcpp::Named("beta") = beta, Rcpp::Named("eta") = eta,
                            Rcpp::Named("objval") = objval);

}
