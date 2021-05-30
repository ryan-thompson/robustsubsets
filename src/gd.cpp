#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// This is a C++ implementation of the projected block-coordinate gradient descent (gd) algorithm
// described in the paper.

// [[Rcpp::export]]
Rcpp::List gd(const arma::mat& X, const arma::vec& y, arma::vec beta, arma::vec eta, const int k,
              const int h, const double step, const int max_iter, const double tol) {

  // Preliminary preparation

  int n = X.n_rows;
  int p = X.n_cols;
  arma::uvec order;
  arma::uvec pos_beta = arma::linspace<arma::uvec>(0, p - k - 1, p - k);
  arma::uvec pos_eta = arma::linspace<arma::uvec>(0, h - 1, h);
  arma::mat tX = trans(X);
  double old_obj;
  double obj = 0.5 * std::pow(arma::norm(y - X * beta - eta, 2), 2);

// Take gradient descent steps

  for (int m = 1; m <= max_iter; m++) {
    old_obj = obj;
    beta = beta + step * tX * (y - X * beta - eta);
    order = arma::sort_index(arma::abs(beta), "ascend");
    beta(order(pos_beta)).fill(0);
    eta = y - X * beta;
    order = arma::sort_index(arma::abs(eta), "ascend");
    eta(order(pos_eta)).fill(0);
    obj = 0.5 * std::pow(arma::norm(y - X * beta - eta, 2), 2);
    if (old_obj - obj <= tol) break;
  }

  // Polish the coefficients

  arma::uvec nz_j = arma::find(beta != 0);
  arma::uvec z_i = arma::find(eta == 0);
  arma::uvec nz_i = arma::find(eta != 0);
  beta(nz_j) = arma::solve(tX(nz_j, z_i) * X(z_i, nz_j), tX(nz_j, z_i) * y(z_i));
  eta = y - X * beta;
  eta(z_i).fill(0);

  // Prepare and return the result

  arma::vec beta_nz = arma::zeros(p);
  beta_nz(nz_j).fill(1);
  arma::vec eta_nz = arma::zeros(n);
  eta_nz(nz_i).fill(1);
  arma::vec x = arma::join_vert(beta, eta, beta_nz, eta_nz);
  obj = 0.5 * std::pow(arma::norm(y - X * beta - eta, 2), 2);
  return Rcpp::List::create(Named("x") = x, Named("objval") = obj);

}
