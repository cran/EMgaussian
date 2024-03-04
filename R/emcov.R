######################################################################
##
## Copyright 2019-2024 Carl F. Falk
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
## <http://www.gnu.org/licenses/>

#' EM algorithm for multivariate normal, covariance matrix parameterization
#' 
#' @param dat Data frame or matrix that contains the raw data.
#' @param max.iter Max number of EM cycles.
#' @param tol Tolerance for change in parameter estimates across EM Cycles. If
#'   all changes are less than \code{tol}, the algorithm terminates.
#' @param start Starting value method (see details).
#' @param debug (Experimental) set an integer value > 1 for some information as the algorithm runs.
#' @param ... Space for additional arguments, not currently used.
#' @details
#' This function computes all means and covariances among a set of variables using
#' the Expectation-Maximization (EM) algorithm to handle missing values, and assuming
#' multivariate normality. The EM code was originally developed for the precision
#' matrix parameterization (\code{\link{em.prec}}), i.e., the parameters are the
#' means and the inverse of the covariance matrix. But, this is easily modifiable
#' to handle a covariance matrix parameterization such that means and covariances
#' are the model parameters.
#' 
#' For starting values, the function accepts either a list that has \code{mu} and
#' \code{S} slots corresponding to the starting mean and covariance matrix. This
#' is useful if the user would like to use custom starting values. Otherwise, a
#' character corresponding to any of the options available in the
#' \code{\link{startvals.cov}} function will be used to take a guess at starting values.
#' 
#' @return
#' A list with the following:
#' \itemize{
#'  \item{\code{p.est}: all parameter estimates as a vector (means followed by unique elements of precision matrix).}
#'  \item{\code{mu}: estimated means.}
#'  \item{\code{S}: estimated covariance matrix.}
#'  \item{\code{it}: number of EM cycles completed.}
#'  \item{\code{conv}: boolean value indicating convergence (TRUE) or not (FALSE).}
#' }
#' 
#' @importFrom stats cov
#' @importFrom lavaan lav_matrix_vechr lav_matrix_vechr_reverse
#' @importFrom matrixcalc is.positive.definite is.symmetric.matrix
#' @importFrom Matrix nearPD forceSymmetric
#' @importFrom Rcpp evalCpp
#' @useDynLib EMgaussian
#' @export
#' @examples
#' \donttest{
#'   library(psych)
#'   data(bfi)
#'   test <- em.cov(bfi[,1:25])
#' }
#' 
##########################################################################################
# EM algorithm originally in Stadler & Buhlmann (2012), which is on missing data w/ Gaussian graphical model
# Try implementing their algorithm. Section 2.3.2
# Just remove the glasso part and we have a saturated mean and covariance matrix
em.cov <- function(dat, max.iter = 500, tol=1e-5, start=c("diag","pairwise","listwise","full"), debug=0,
                         ...){

  dat <- as.matrix(dat)
  
  # obtain starting values  
  if(is.list(start) && is.vector(start$mustart) && is.matrix(start$covstart)){
    mustart <- start$mustart
    covstart <- start$covstart
  } else if (is.character(start)) {
    start <- startvals.cov(dat, start)
    mustart <- start$mustart
    covstart <- start$covstart
  }

  # track change in means and covariances
  mu.est <- as.matrix(mustart)
  S.est.mat <- covstart
  S.est <- lav_matrix_vechr(S.est.mat)
  p.est <- c(mustart, S.est)

  # loop until convergence
  conv<-FALSE
  for(it in 1:max.iter){
    if(debug>1){
      print(paste0('iter ',it))
    }
    
    # Do one EM cycle
    EMres <- EMcyclecov(dat, mu.est, S.est.mat)
    mu.est <- EMres$mu
    S.est.mat <- EMres$S
    
    # Check S
    if(!is.symmetric.matrix(S.est.mat)){
      S.est.mat<-as.matrix(forceSymmetric(S.est.mat))
      warnings("Non-symmetric S found after E step. Could indicate identification or estimation problem.")
    }
    
    # force positive-definite S
    if(!is.positive.definite(S.est.mat)){
      S.est.mat<-as.matrix(nearPD(S.est.mat)$mat)
      warnings("Non-positive definite S found after E step. Could indicate identification or estimation problem.")
    }
    
    # Save estimates
    S.est <- lav_matrix_vechr(S.est.mat)
    p.old <- p.est
    p.est <- c(mu.est,S.est)
    
    # check for convergence
    if(debug>2){print(max(abs(p.old-p.est)))}
    if(max(abs(p.old-p.est)) < tol){
      conv<-TRUE
      break
    }
  }
  
  out<-list(p.est=p.est, mu=mu.est, S=S.est.mat, it=it, conv=conv)
  return(out)
  
}


#' Starting values for means and covariances
#' 
#' @param dat Data frame or matrix that contains the raw data.
#' @param start Starting value method (see details).
#' @details
#' Attempts to figure out a starting values for the means and covariances for use
#' with other functions that do the EM algorithm such as \code{\link{em.prec}} or
#' \code{\link{em.cov}}. Note that means are determined univariately using all
#' available cases. For covariances, several options are available:
#' 
#' - "diag" Use all available complete values to compute the variances of each variable and construct a diagonal covariance matrix.
#' - "pairwise" Pairwise (co)variances will be used to construct the starting covariance matrix.
#' - "listwise" Listwise deletion will be used and only those with complete data will contribute to the starting covariance matrix.
#' - "full" Cheat and use \code{lavaan} to obtain direct maximum likelihood estimates of covariances. This defeats the purpose to some extent, but not that \code{lavaan} may be quite slow compared to this implementation.
#' 
#' @return
#' A list consisting of:
#' \itemize{
#'   \item{\code{mustart} - starting values for means.}
#'   \item{\code{covstart} - starting values for covariances.}
#' }
#' @importFrom stats cov
#' @importFrom lavaan lavCor
#' @importFrom matrixcalc is.positive.definite
#' @importFrom Matrix nearPD
#' @export
#' @examples
#' \donttest{
#'   library(psych)
#'   data(bfi)
#'   startvals.cov(bfi[,1:25])
#' }
startvals.cov <- function(dat, start=c("diag","pairwise","listwise","full")){

  start <- match.arg(start)
  
  dat <- as.matrix(dat)
  
  # obtain starting values
  mustart <- colMeans(dat,na.rm=T)
  if(start=="listwise"){
    covstart <- cov(dat, use="complete.obs")    
  } else if (start=="pairwise"){
    covstart <- cov(dat,use="pairwise.complete.obs")
  } else if (start=="full"){
    covstart <- lavCor(as.data.frame(dat), missing="ml", output="cov")
  } else if (start=="diag"){
    covstart <- diag(diag(cov(dat,use="pairwise.complete.obs")))
  }
  
  if(any(is.na(covstart))){
    covstart[is.na(covstart)]<-0
    warning("Some NAs in starting vals for covariance matrix. Replacing with zeros.")
  }
  
  # force positive-definite covstart
  if(!is.positive.definite(covstart)){
    covstart<-as.matrix(nearPD(covstart)$mat)
  }
  
  out <- list(mustart = mustart,
              covstart = covstart)
  
  return(out)
  
}