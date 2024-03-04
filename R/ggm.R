######################################################################
## Functions for estimating regularized network models with
## missing data.
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
##
## Notice of other contributions:
##
## Some lines regarding computing pairwise counts are similar
## to those in bootnet (https://cran.r-project.org/package=bootnet) e.g., V. 1.5.1
## Those lines are acknowledged below (search for bootnet).
##
## Some code regarding creating the grid for the tuning parameter is from qgraph
## (https://cran.r-project.org/package=qgraph) or glasso
## (https://cran.r-project.org/package=glasso)
## Refer to the rhogrid function.
##
## At the time this package was shared publicly, all other code was released under GPL-2:
##  qgraph is Copyright (C) 2010 - 2019 Sacha Epskamp
##  bootnet is authored by Sacha Epskamp and Eiko I. Fried
##  glasso is authored by Jerome Friedman, Trevor Hastie and Rob Tibshirani
##

##########################################################################################

#' Regularized GGM under missing data with EBIC or k-fold cross validation
#' 
#' @param dat A data frame or matrix of the raw data.
#' @param max.iter Maxumum number of EM cycles for the EM algorithm, passed eventually to \code{\link{em.prec}}.
#' @param est.tol Tolerance for change in parameter estimates across EM Cycles. If
#'   all changes are less than \code{tol}, the algorithm terminates.
#' @param start Starting value method (see details of \code{\link{em.prec}}). Note that "none" will not do
#'   any regularization, making \code{rho} moot.
#' @param glassoversion Character indicating which function to use at the M-step for regularization. See 
#'   \code{\link{em.prec}}) for more details.
#' @param rho Vector of tuning parameter values. Defaults to just a single value of 0 (i.e., also no regularization).
#' @param rhoselect Method of selecting the tuning parameter. "ebic" will select the best value in \code{rho} based
#'   on that yielding the best EBIC. "kfold" will do k-fold cross-validation.
#' @param N Sample size to use in any EBIC calculations. Defaults to average pairwise sample size.
#' @param gam Value of gamma in EBIC calculations. Typical choice, and the default, is .5.
#' @param zero.tol Tolerance in EBIC calculations for declaring edges to be zero. Anything in absolute value
#'   below \code{zero.tol} will be considered zero and will not count as a parameter in EBIC calculations.
#' @param k If using k-fold cross validation, an integer specifying the number of folds. Defaults to 5.
#' @param seed Random number seed passed to function that does k-fold cross validation. Use if you want folds to
#'   be more or less replicable.
#' @param debug (Experimental) Pick an integer above 0 for some messages that may aide in debugging.
#' @param convfail (Experimental) What to do if optimization fails. If TRUE, will fill in the partial correlation matrix with 0's.
#' @param ... Other arguments passed down to \code{\link{em.prec}} and to functions that do the glasso.
#' @details
#' This function will estimate a regularized Gaussian graphical model (GGM) using either EBIC or k-fold
#' cross validation to select the tuning parameter. It was written as a wrapper function to the
#' \code{\link{em.prec}} function, which is an implementation of the EM algorithm from
#' Städler & Bühlmann (2012). These authors also used k-fold cross validation, which is implemented here.
#' 
#' For the tuning parameter (\code{rho}), typically a grid of of values is evaluated and the one
#' that results in the best EBIC or cross validation metric is selected and returned as the final model.
#' See \code{\link{rhogrid}} and examples for some ways to pick these candidate tuning parameter values.
#' 
#' This function is intended to be compatible (to my knowledge) with \code{\link[bootnet]{estimateNetwork}}
#' in the \code{bootnet} package in that one can input this as a custom function and then have all of the
#' benefits of plotting, centrality indices, and so on, that \code{bootnet} provides.
#' 
#' Most methods in examples were studied by Falk and Starr (under review). In particular use of this
#' function in conjunction with \code{\link[bootnet]{estimateNetwork}} worked well for both "ebic"
#' and "kfold" for model selection, though the original article used the \code{\link[glassoFast]{glassoFast}}
#' package for estimation, which by default penalizes the diagonal of the precision matrix. It seems
#' slightly more common to not penalize the diagonal. Therefore, the below use \code{\link[glasso]{glasso}};
#' this approach performed well but occasionally would get stuck while trying to find an optimal solution.
#' In addition, the two-stage approach studied by Falk and Starr (under review) also performed well, though
#' not as good as the present function; that approach is available in the bootnet package. Note also that
#' \code{\link[cglasso]{cglasso}} has an implementation of Städler & Bühlmann (2012), but we found in
#' simulations with a high proportion of missing data that our implementation was less likely to 
#' encounter estimation problems.
#' 
#' @references
#' 
#' Evidence that EBIC and k-fold works well with this package (cite if you use the EMglasso R package):
#' Falk, C. F., & Starr, J. (2023, July 19). Regularized cross-sectional network modeling with missing data:
#'   A comparison of methods. Preprint: \doi{10.31219/osf.io/dk6zv}
#'
#' Original publication on use of EM algorithm and k-fold cross-validation for glasso:
#' Städler, N., & Bühlmann, P. (2012). Missing values: sparse inverse covariance estimation and an
#'   extension to sparse regression. Statistics and Computing,22, 219–235. \doi{10.1007/s11222-010-9219-7}
#'
#' If you use this package with bootnet:
#' Epskamp,S., Borsboom, D., &  Fried, E. I. (2018). Estimating psychological networks and their accuracy:
#'   A tutorial paper. Behavior research methods, 50(1), 195–212.  \doi{10.3758/s13428-017-0862-1}
#'
#' If you also use this package with qgraph:
#' Epskamp, S., Cramer, A., Waldorp, L., Schmittmann, V. D., & Borsboom, D. (2012). qgraph: Network visualizations of relationships in psychometric data. Journal of Statistical Software, 48 (1), 1-18.
#'	
#' Foundational article on glasso:
#' Friedman, J. H., Hastie, T., & Tibshirani, R. (2008). Sparse inverse covariance estimation with the graphical lasso.
#'   Biostatistics, 9 (3), 432-441.
#' 
#' Foundational article on use of EBIC with glasso:
#' Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria for Gaussian graphical models.
#'
#' If glasso is also the estimation method:
#' Friedman, J. H., Hastie, T., & Tibshirani, R. (2014). glasso: Graphical lasso estimation of Gaussian graphical models.
#'   Retrieved from https://CRAN.R-project.org/package=glasso
#' 
#' If glassoFast is also the estimation method: 
#' Sustik M.A., Calderhead B. (2012). GLASSOFAST: An efficient GLASSO implementation. UTCS Technical Report TR-12-29:1-3.
#' 
#' @return A list with the following elements:
#' 
#' \itemize{
#'  \item{\code{results}: Contains the "best" or selected model, with elements mostly following that of \code{\link{em.prec}}'s output.
#'    Additional elements include \code{rho} and \code{crit}, which are the grid of tuning parameter values and value of the
#'    criterion. These are added here because if used with \code{\link[bootnet]{bootnet}}, it may not save the other slots below.}
#'  \item{\code{rho}: Grid for \code{rho}.}    
#'  \item{\code{crit}: Vector of the same length as the grid for \code{rho} with the criterion (EBIC or cross-validation).
#'    Smaller is better.}
#'  \item{\code{graph}: Estimated partial correlation matrix; \code{\link[bootnet]{bootnet}} appears to expect this in order
#'    to do plots, compute centrality indices and so on.}
#' }
#' 
#' @importFrom stats cov2cor
#' @export
#' @examples
#' \donttest{
#'   library(psych)
#'   library(bootnet)
#'   data(bfi)
#'
#'   # Regularized estimation with just a couple of tuning parameter values
#'   
#'   # EBIC
#'   rho <- seq(.01,.5,length.out = 50)
#'   ebic1 <- EMggm(bfi[,1:25], rho = rho, glassoversion = "glasso",
#'                  rhoselect = "ebic")
#'   
#'   # k-fold
#'   kfold1 <- EMggm(bfi[,1:25], rho = rho, glassoversion = "glasso",
#'                   rhoselect = "kfold")   
#'
#'   plot(rho, ebic1$crit) # values of EBIC along grid
#'   plot(rho, kfold1$crit) # values of kfold along grid
#'   
#'   # partial correlation matrix
#'   # ebic1$graph 
#'   # kfold1$graph
#'
#'   # Integration with bootnet package
#'   ebic2 <- estimateNetwork(bfi[,1:25], fun = EMggm, rho = rho,
#'                            glassoversion = "glasso", rhoselect = "ebic")
#'   kfold2 <- estimateNetwork(bfi[,1:25], fun = EMggm, rho = rho,
#'                            glassoversion = "glasso", rhoselect = "kfold")   
#'   
#'   # ebic2 and kfold2 now do just about anything one could normally do with an
#'   # object returned from estimateNetwork e.g., plotting
#'   plot(ebic2)
#'   plot(kfold2)
#'   
#'   # e.g., centrality indices
#'   library(qgraph)
#'   centralityPlot(ebic2)
#'   # and so on
#'
#'   # Other ways to pick grid for tuning parameter... 100 values using method
#'   # that qgraph basically uses; requires estimate of covariance matrix. Here
#'   # the covariance matrix is obtained directly from the data.
#'   rho <- rhogrid(100, method="qgraph", dat = bfi[,1:25])
#'  
#'   ebic3 <- estimateNetwork(bfi[,1:25], fun = EMggm, rho=rho,
#'                            glassoversion = "glasso", rhoselect = "ebic")
#'   kfold3 <- estimateNetwork(bfi[,1:25], fun = EMggm, rho=rho,
#'                            glassoversion = "glasso", rhoselect = "kfold")  
#'   
#'   # look at tuning parameter values vs criterion
#'   plot(ebic3$result$rho, ebic3$result$crit)
#'   plot(kfold3$result$rho, kfold3$result$crit)
#'   
#'   # EBIC with bootnet; does listwise deletion by default
#'   ebic.listwise <- estimateNetwork(bfi[,1:25], default="EBICglasso")
#'   
#'   # EBIC with bootnet; two-stage estimation
#'   # Note the constant added to bfi tricks lavaan into thinking the data are
#'   # not categorical
#'   ebic.ts <- estimateNetwork(bfi[,1:25] + 1e-10, default="EBICglasso",
#'                              corMethod="cor_auto", missing="fiml")
#'   
#' }
EMggm<-function(dat,
                max.iter = 500, est.tol = 1e-7, start=c("diag","pairwise","listwise","full"), glassoversion=c("glasso","glassoFast","glassonostart","none"), # em.prec options
                rho = 0, rhoselect = c("ebic","kfold"), # how to select model
                N=NULL, gam = .5, zero.tol = 1e-10, # stuff for EBIC
                k = 5, seed=NULL, # for kfold
                debug=0,
                convfail=FALSE, ...){
  
  rhoselect<-match.arg(rhoselect)
  
  if(rhoselect=="ebic"){
    
    # how to obtain "N"?
    if(is.null(N)){
      # next 3 lines copied and modified from bootnet, which does pairwise_average by default but accidentally counts the univariate (diagonal) elements
      xmat <- as.matrix(!is.na(dat))
      misMatrix <- t(xmat) %*% xmat
      N<-mean(misMatrix[lower.tri(misMatrix)])
    }
    
    # estimate models
    mods <- list()
    for(m in 1:length(rho)){
      if(debug>0){
        print(paste0('model', m))
      }
      mods[[m]] <- em.prec(dat,  max.iter = max.iter, tol=est.tol, start=start, debug=debug, glassoversion=glassoversion, rho=rho[m], ...)
    }
    
    # Calculate EBIC for all models
    crit<-vector("numeric")
    for(m in 1:length(rho)){
      if(convfail & !mods[[m]]$conv){
        val <- NA
      } else {
        val <- EBICggm(mods[[m]]$p.est, dat, N, gam=gam, tol=zero.tol)
      }
      crit<-c(crit,val)
    }
    
    # select the best model, and return it in the "graph" slot along w/ results
    best.idx<-which.min(crit)
    
    best.mod<-mods[[best.idx]]
  } else if (rhoselect == "kfold"){
    
    # do k-fold cross validation to select rho
    cv <- fiml.ggm.cv(dat, rho, max.iter, est.tol, start, V=k, seed=seed, debug=debug, glassoversion=glassoversion, ...)
    crit <- cv$rho.result

    # then, estimate model w/ best rho
    best.mod <- em.prec(dat, max.iter = max.iter, tol=est.tol, start=start, debug=debug, glassoversion=glassoversion, rho=cv$best.rho, ...)

    # troubleshoot if it did not converge
    if(!best.mod$conv){
      for(j in 1:(sum(!is.na(crit))-1)){
        warning("Best option for tuning parameter according to cross-validation did not converge; trying on next best tuning parameter value.")
        nextrho <- rho[order(crit)[j+1]]
        best.mod <- em.prec(dat, max.iter = max.iter, tol=est.tol, start=start, debug=debug, glassoversion=glassoversion, rho=nextrho, ...)
        if(best.mod$conv) break
      }
    }
  }
  
  out<-list()
  out$results <- best.mod
  out$results$rho <- rho
  out$results$crit <- crit
  out$rho <- rho
  out$crit <- crit
  
  pcor.net <- -cov2cor(best.mod$K)
  diag(pcor.net)<-0
  colnames(pcor.net)<-rownames(pcor.net)<-colnames(dat)
  pcor.net<-as.matrix(forceSymmetric(pcor.net))
  out$graph <- pcor.net
  
  # if convergence failed for the given model, fill it with zero's
  # Although it appears that's what bootnet/qgraph does with two-stage estimation when lavCor fails,
  # this may not actually be the case.
  if(convfail & !best.mod$conv){
    out$graph[]<-0
  }

  return(out)
}


# k-fold cross-validation
#' @importFrom caret createFolds
fiml.ggm.cv <- function(dat, rho,  max.iter = 500, est.tol = 1e-7, start=c("diag","pairwise","listwise","full"), V=10, seed=NULL, debug=0, glassoversion, ...){
  if(!is.null(seed)) {set.seed(seed)}
  N <- nrow(dat)
  
  # create folds
  folds <- createFolds(1:N, k=V)
  
  rho.result<-vector("numeric")
  
  # loop over folds
  for(r in 1:length(rho)){
    
    if(debug>0){
      print(paste0('model', r))
    }
    
    test.result <- 0
    for(i in 1:V){
      train <- dat[-folds[[i]],]
      test <- dat[folds[[i]],]
      
      train.result<-try(em.prec(train, max.iter = max.iter, tol = est.tol, start = start, debug=debug, glassoversion=glassoversion, rho = rho[[r]], ...))
      if(!inherits(train.result,"try-error") && train.result$conv){
        test.result <- test.result + nllggm.wrap(train.result$p.est, test)
      } else {
        test.result <- NA
      }
    }
    
    # save test.result
    rho.result<-c(rho.result,test.result)
  }
  
  out<-list()
  out$best.rho <- rho[which.min(rho.result)]
  out$rho.result <- rho.result
  out$rho <- rho
  
  return(out)
}

#' @importFrom lavaan lav_matrix_vechr_reverse
#' @useDynLib EMgaussian
nllggm.wrap<-function(p, dat){
  
  ni <- ncol(dat)
  np <- length(p)
  mu <- p[1:ni]
  K <- lav_matrix_vechr_reverse(p[(ni+1):np])
  
  nll <- nllprec(as.matrix(dat), mu, K)
  return(nll)
  
}

# EBIC
#' @importFrom lavaan lav_matrix_vechr_reverse
EBICggm <- function(p, dat, N=NULL, gam=.5, tol=1e-32){
  
  # how to obtain "N"?
  if(is.null(N)){
    # next 3 lines copied and modified from bootnet, which does pairwise_average by default but accidentally counts the univariate (diagonal) elements
    xmat <- as.matrix(!is.na(dat))
    misMatrix <- t(xmat) %*% xmat
    N<-mean(misMatrix[lower.tri(misMatrix)])
  }
  
  # -2L
  neg2l <- 2*nllggm.wrap(p, dat)
  
  # number of "nodes"
  P <- ncol(dat)
  
  # how do we know how many non-zero edges there are?
  
  # actually look at off-diagonal elements of K
  K <- lav_matrix_vechr_reverse(p[(P+1):length(p)])
  E<-sum(abs(K[lower.tri(K)])>tol)
  
  out <- neg2l + E*log(N) + 4*gam*E*log(P)
  
  return(out)
  
}

#' Create sequence of possible tuning parameter values
#' 
#' @param n.rho Integer corresponding to the number of tuning parameter values.
#' @param method Character corresponding to the method to create tuning parameter values ("qgraph" or "glassopath"); see Details.
#' @param rho.min.ratio Numeric value that mimics \code{\link[qgraph]{EBICglasso}} behavior for tuning parameter.
#'   i.e., "Ratio of lowest (tuning parameter) compared to maximal (tuning parameter)".
#' @param dat The raw data, if \code{S} is not provided used to estimate \code{S}. Not required if \code{S} is provided.
#' @param S Covariance matrix for the data. If provided, supersedes \code{dat}.
#' @param ... Other arguments passed down to \code{\link{em.prec}}.
#' @details
#' For regularized estimation of the Gaussian graphical model, a sequence or grid of possible tuning
#' parameter values is often tried, with the tuning parameter that optimizes some criterion (EBIC, k-fold
#' cross validation) chosen. This is an attempt to automate some of the sequence creation. Code
#' is borrowed from qgraph and glasso, acknowledged in references below. Both require some estimate
#' of the covariance matrix in order to do regularization; if not provided, \code{\link{em.prec}} with
#' default options is attempted.
#' 
#' For "qgraph" the max value is determined by the maximum absolute value of the
#' difference between the covariance matrix and an identity matrix. The min is \code{rho.min.ratio}
#' times the max value. A sequence that is equally spaced on a log scale is then constructed between
#' these two values.
#' 
#' For "glassopath", the max value is the max absolute value of the covariance matrix. The min is
#' the max divided by the number of desired tuning parameter values. A sequence that is equally spaced
#' between these two values is then constructed.
#' 
#' @references
#' qgraph:
#' Epskamp, S., Cramer, A., Waldorp, L., Schmittmann, V. D., & Borsboom, D. (2012). qgraph: Network visualizations of relationships in psychometric data. Journal of Statistical Software, 48 (1), 1-18.
#'	
#' glasso (i.e., glassopath option):
#' Friedman, J. H., Hastie, T., & Tibshirani, R. (2014). glasso: Graphical lasso estimation of Gaussian graphical models.
#'   Retrieved from https://CRAN.R-project.org/package=glasso
#' 
#' @return A vector of possible tuning parameter values.
#' 
#' @export
#' @examples
#' \donttest{
#' library(psych)
#' data(bfi)
#'
#' # pick 50 values using the approach qgraph uses; give data as input 
#' rho <- rhogrid(50, method="qgraph", dat = bfi[,1:25])
#' 
#' emresult <- em.cov(bfi[,1:25])
#' S<-emresult$S
#' 
#' # pick 50 values using the approach glasso uses; give S as input
#' rho2 <- rhogrid(50, method="glassopath", S = S)
#' 
#' }
# Helper function to pick some values for rho
rhogrid <- function(n.rho, method=c("qgraph","glassopath"), rho.min.ratio = .01,
                    dat = NULL, S = NULL, ...){
  
  method <- match.arg(method)
  
  if(is.null(S) & is.null(dat)){
    stop("Either provide raw data (dat) or some estimate of the covariance matrix (S)")
  }
  
  if(is.null(S)){
    # estimate saturated covariance matrix using EM algorithm
    sat <- em.cov(dat, ...)
    if(!sat$conv){
      warning("Estimation of covariance matrix may not have converged.")
    }
    S <- sat$S
  }
  
  if (method == "qgraph") {
    # This is basically from qgraph, EBICglassoCore function
    rho.max = max(max(S - diag(nrow(S))), -min(S - diag(nrow(S))))
    rho.min = rho.min.ratio * rho.max
    rho = exp(seq(log(rho.min), log(rho.max), length = n.rho))
  } else if (method =="glassopath"){
    # This is from the glasso package, glassopath function
    rho = seq(max(abs(S))/n.rho, max(abs(S)), length = n.rho)
  }

  return(rho)

}