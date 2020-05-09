#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

/* This is a C++ implementation of the projected block-coordinate gradient descent (pbgd) algorithm
 * described in the paper.
 */

// [[Rcpp::export]]
List pbgd(arma::mat X, arma::vec y, arma::vec beta, arma::vec eta, int k, int h, double step,
          int max_iter, double tol) {

  /* Preliminary preparation */

  int n = X.n_rows;
  int p = X.n_cols;
  arma::uvec order;
  arma::uvec pos_beta = arma::linspace<arma::uvec>(0, p - k - 1, p - k);
  arma::uvec pos_eta = arma::linspace<arma::uvec>(0, h - 1, h);
  arma::mat tX = trans(X);
  double old_obj;
  double obj = 0.5 * dot(y - X * beta - eta, y - X * beta - eta);

  /* Take gradient descent steps */

  for (int m = 1; m <= max_iter; m++) {
    old_obj = obj;
    beta = beta + step * tX * (y - X * beta - eta);
    order = sort_index(abs(beta), "ascend");
    beta.elem(order(pos_beta)).fill(0);
    eta = y - X * beta;
    order = sort_index(abs(eta), "ascend");
    eta.elem(order(pos_eta)).fill(0);
    obj = 0.5 * dot(y - X * beta - eta, y - X * beta - eta);
    if (old_obj - obj <= tol) {
      break;
    }
  }

  /* Polish the coefficients */

  arma::uvec nz_j = find(beta != 0);
  arma::uvec z_i = find(eta == 0);
  arma::uvec nz_i = find(eta != 0);
  beta(nz_j) = solve(tX(nz_j, z_i) * X(z_i, nz_j), tX(nz_j, z_i) * y(z_i));
  eta = y - X * beta;
  eta(z_i).fill(0);

  /* Finally, prepare and return the result */

  arma::vec beta_nz = arma::zeros(p);
  beta_nz(nz_j).fill(1);
  arma::vec eta_nz = arma::zeros(n);
  eta_nz(nz_i).fill(1);
  arma::vec x = join_vert(beta, eta, beta_nz, eta_nz);
  obj = 0.5 * dot(y - X * beta - eta, y - X * beta - eta);
  return List::create(Named("x") = x, Named("objval") = obj);

}
