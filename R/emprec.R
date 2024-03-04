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

#' EM algorithm for multivariate normal, precision matrix parameterization
#' 
#' @param dat Data frame or matrix that contains the raw data.
#' @param max.iter Max number of EM cycles.
#' @param tol Tolerance for change in parameter estimates across EM Cycles. If
#'   all changes are less than \code{tol}, the algorithm terminates.
#' @param start Starting value method (see details).
#' @param glassoversion Character indicating whether to do regularization (lasso),
#'   and if so, using which package. "glasso" uses the \code{\link[glasso]{glasso}}
#'   function and uses the E-step covariance matrix for starting values,
#'   "glassoFast" uses \code{\link{glassoFast}{glassoFast}} which also
#'   penalizes the diagonal of the precision matrix by default (glasso does not),
#'   "glassonostart" also uses \code{\link[glasso]{glasso}} but no "warm" starting
#'   values.
#' @param debug (Experimental) set an integer value > 1 for some information as the algorithm runs.
#' @param ... Arguments passed down to any of the glasso functions.
#' @details
#' This function computes all means and the precision matrix (inverse of covariance
#' matrix) among a set of variables using the Expectation-Maximization (EM)
#' algorithm to handle missing values, and assuming multivariate normality. The
#' EM code was originally developed based on Stadler and Buhlmann (2012) and for
#' use with the graphical lasso (i.e., glasso). This version allows the possibility
#' of using a lasso by specifying something other than "none" for \code{glassoversion}.
#' However, it can also be used without regularization to just estimate the precision matrix.
#' 
#' For starting values for the EM algorithm itself (not at the M-step), the
#' function accepts either a list that has \code{mu} and \code{S} slots
#' corresponding to the starting mean and covariance matrix. This is useful if
#' the user would like to use custom starting values. Otherwise, a character
#' corresponding to any of the options available in the \code{\link{startvals.cov}}
#' function will be used to take a guess at starting values.
#' 
#' @references Städler, N., & Bühlmann, P. (2012). Missing values: sparse inverse covariance estimation and an
#'   extension to sparse regression. Statistics and Computing, 22, 219–235. \doi{10.1007/s11222-010-9219-7}
#' 
#' @return 
#' A list with the following:
#' \itemize{
#'  \item{\code{p.est}: all parameter estimates as a vector (means followed by unique elements of precision matrix).}
#'  \item{\code{mu}: estimated means.}
#'  \item{\code{S}: estimated covariance matrix.}
#'  \item{\code{K}: estimated precision matrix.}
#'  \item{\code{it}: number of EM cycles completed.}
#'  \item{\code{conv}: boolean value indicating convergence (TRUE) or not (FALSE).}
#' }
#' 
#' @importFrom stats cov
#' @importFrom lavaan lav_matrix_vechr lav_matrix_vechr_reverse
#' @importFrom matrixcalc is.positive.definite is.symmetric.matrix
#' @importFrom Matrix nearPD forceSymmetric
#' @importFrom glasso glasso
#' @importFrom glassoFast glassoFast
#' @importFrom Rcpp evalCpp
#' @useDynLib EMgaussian
#' @export
#' @examples
#' \donttest{
#'   library(psych)
#'   data(bfi)
#'   test <- em.prec(bfi[,1:25])
#' }
##########################################################################################
# Städler, N., & Bühlmann, P. (2012). Missing values: sparse inverse covariance estimation and an extension to sparse regression.Statistics and Computing,22, 219–235.  doi:10.1007/s11222-010-9219-7
# Try implementing their algorithm. Section 2.3.2
# Note that setting the glasso method to "none" results in no regularization
em.prec <- function(dat, max.iter = 500, tol=1e-5, start=c("diag","pairwise","listwise","full"),
                    glassoversion = c("none","glassoFast","glasso","glassonostart"),
                    debug=0,
                     ...){
  
  dat <- as.matrix(dat)
  N <- nrow(dat)
  glassoversion <- match.arg(glassoversion)
  
  # obtain starting values  
  if(is.list(start) && is.vector(start$mustart) && is.matrix(start$covstart)){
    mustart <- start$mustart
    covstart <- start$covstart
  } else if (is.character(start)) {
    start <- startvals.cov(dat, start)
    mustart <- start$mustart
    covstart <- start$covstart
  }

  Kstart.mat<-solve(covstart)
  Kstart <- lav_matrix_vechr(Kstart.mat)
  
  mu.est <- as.matrix(mustart)
  
  K.est <- Kstart
  K.est.mat <- lav_matrix_vechr_reverse(K.est)

  pkstart<-c(mustart,Kstart)
  p.est<-pkstart

  # loop until convergence
  conv<-FALSE
  for(it in 1:max.iter){
    if(debug>1){
      print(paste0('iter ',it))
    }

    # Do one EM cycle
    EMres <- EMcycleprec(dat, mu.est, K.est.mat)
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
    
    if(glassoversion=="glasso"){
      gres <- glasso(S.est.mat, nobs=N, w.init=S.est.mat, wi.init=solve(S.est.mat), start="warm",...)
    } else if (glassoversion=="glassoFast"){
      gres <- glassoFast(S.est.mat, w.init=S.est.mat, wi.init=solve(S.est.mat), start="warm",...)
    } else if (glassoversion=="glassonostart"){
      gres <- glasso(S.est.mat, nobs=N, ...)
    } else if (glassoversion=="none"){
      gres <- list(wi = solve(S.est.mat))
    }
    
    # Save estimates
    K.est.mat <- K.est.mat <- gres$wi #solve(S.est.mat)
    K.est <- lav_matrix_vechr(K.est.mat)
    p.old <- p.est
    p.est <- c(mu.est,K.est)
    
    # check for convergence
    if(debug>2){print(max(abs(p.old-p.est)))}
    if(max(abs(p.old-p.est)) < tol){
      conv<-TRUE
      break
    }
  }
  
  out<-list(p.est=p.est, mu=mu.est, S=S.est.mat, K=K.est.mat, it=it, conv=conv)
  return(out)
  
}