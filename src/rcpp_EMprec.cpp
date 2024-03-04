// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

/*
 Copyright 2019-2024 Carl F. Falk
 
 This program is free software: you can redistribute it and/or
 modify it under the terms of the GNU General Public License as
 published by the Free Software Foundation, either version 3 of
 the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 <http://www.gnu.org/licenses/>
 */

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// First part of E step for an entire data matrix
// [[Rcpp::export]]
arma::mat imp1matprec(Rcpp::NumericMatrix D, const arma::colvec& muest, const arma::mat& kest){
    
  arma::mat d(D.begin(), D.rows(), D.cols(), true); // I think false means no copy; here we want one to avoid replacing original data
  
  int N = d.n_rows;

  for(int i = 0; i<N; i++){
    arma::uvec io = find_finite(d.row(i));
    arma::uvec im = find_nonfinite(d.row(i));
    arma::uvec idx = {static_cast<arma::uword>(i)};

    d.submat(idx,im) = (muest(im) - inv(kest.submat(im,im)) * kest.submat(im,io) * (d.submat(idx,io).t() - muest(io))).t();        
  }
  return(d);
}
      
// Second part of E step for an entire data matrix
// operates directly on T2 to avoid making a copy
// [[Rcpp::export]]
void imp2matprec(Rcpp::NumericMatrix D, const arma::mat& kest, arma::mat& t2){
  
  arma::mat d(D.begin(), D.rows(), D.cols(), true);
  
  int N = d.n_rows;

  for(int i = 0; i<N; i++){
    arma::uvec im = find_nonfinite(d.row(i));
    t2.submat(im,im) += inv(kest.submat(im,im));
  }
  
}

// EM cycle all in one shot, covariance matrix parameterization
// [[Rcpp::export]]
Rcpp::List EMcycleprec(const Rcpp::NumericMatrix& D, const arma::colvec& muest, const arma::mat& kest){
  
  // First part of imputation
  arma::mat dimp = imp1matprec(D, muest, kest);

  int N = dimp.n_rows;
  
  // Vector of 1's  
  arma::colvec ones = arma::ones<arma::vec>(N);
  
  // Compute T1 matrix
  arma::mat t1 = dimp.t() * ones;

  // Compute T2 matrix, then add stuff to T2
  arma::mat t2 = dimp.t() * dimp;

  // Then, add stuff to T2
  imp2matprec(D, kest, t2);

  // Now, compute mu
  arma::colvec newmu = t1/N;

  // update S and optionally K
  arma::mat newS = t2/N - newmu*newmu.t();
  arma::mat newK = inv(newS);
  
  Rcpp::List out;
  out["mu"] = newmu;
  out["S"] = newS;
  out["K"] = newK;

  return out;
}
        
// Negative log-likelihood; precision matrix parameterization
// [[Rcpp::export]]
double nllprec(const arma::mat d, const arma::colvec muest, const arma::mat kest){
    
    int N = d.n_rows;
    int J = d.n_cols;
    
    arma::mat kinv = inv(kest);
    
    double nll = 0;
    
    for(int i = 0; i<N; i++){
      arma::uvec io = find_finite(d.row(i));
      arma::uvec idx = {static_cast<arma::uword>(i)};
      nll += (0.5*(log(det(kinv.submat(io,io))) + (d.submat(idx,io).t() - muest(io)).t() * inv(kinv.submat(io,io)) * (d.submat(idx,io).t() - muest(io)) + J*log(2*arma::datum::pi))).eval()(0,0);
    }
    
    return(nll);
}
