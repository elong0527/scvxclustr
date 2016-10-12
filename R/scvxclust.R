#' @useDynLib scvxclustr
#' @importFrom Rcpp sourceCpp
#' @import cvxclustr
#' @importFrom stats dist
#' @importFrom stats rnorm
#' @import gglasso
#' @import mclust
#' @import Matrix
#' @import igraph
NULL

#' Sparse Convex Clustering Path
#'
#' The function estimates the sparse convex clustering path via AMA or ADMM.
#' Required inputs include a data matrix \code{X} (rows are samples; columns are features), a vector of weights
#' \code{w}, regularization parameters \code{Gamma1}, \code{Gamma2} and the adaptive weight \code{Gamma2_weight}. 
#'
#' @param X The data matrix to be clustered. The rows are the samples, and the columns are the features.
#' @param w A vector of nonnegative weights. The ith entry \code{w[i]} denotes the weight used between the ith pair of centroids. The weights are in dictionary order.
#' @param Gamma1 A regularization parameter controls cluster size .
#' @param Gamma2 A regularization parameter controls the number of informative features .
#' @param Gamma2_weight The weight  to adaptively penalize the features.
#' @param nu A positive penalty parameter for quadratic deviation term.
#' @param tol_abs The convergence tolerance (absolute).
#' @param tol_rel The convergence tolerance (relative).
#' @param max_iter The maximum number of iterations.
#' @param type An integer indicating the norm used: 2 = 2-norm. (Only L2 norm are supported for now)
#' @param verbose report convergence information
#' @param method method to fit the sparse convex clustering ("ama" or "admm"). Default is ama
#' @param init initial vlaue of the method
#' @return \code{U} A list of centroid matrices.
#' @return \code{V} A list of centroid difference matrices.
#' @return \code{Lambda} A list of Lagrange multiplier matrices.
#' @return \code{iters} number of iterations.
#' @return \code{eva}   the absolute difference of U between two most recent iteration.
#' @return \code{method} fitted method ("ama" or "admm")
#'  
#' @export
scvxclust <- function(X,w,Gamma1, Gamma2, Gamma2_weight, nu=1,tol_abs=1e-3,tol_rel=1e-4,max_iter=1e4,type=2, verbose = F, method = "ama", init = NULL){

  call <- match.call()
  nGamma <- length(Gamma1)

  list_A <- vector(mode="list",length=nGamma)
  list_V <- vector(mode="list",length=nGamma)
  list_Lambda <- vector(mode="list",length=nGamma)
  iter_vec <- integer(nGamma)

  gamma2 <- Gamma2

  if(! (type %in% 2) ){
    stop("Only type equals 2 are sportted")
  }
  if(! (method %in% c("ama", "admm"))){
    stop("Only method in `ama` or `admm` are sportted")
  }


  ###--------------------------
  ### AMA Part
  ###--------------------------
  if( method == "ama"){
    for (ig in 1:nGamma) {

      edge_info <- compactify_edges(w, n, method = "ama")
      nK <- length(which(w > 0))
      ix <- edge_info$ix

      gamma1 <- Gamma1[ig]
      fit <- ama_eigen( t(X), w[w > 0],ix-1, gamma1, gamma2, Gamma2_weight, nu, tol_abs, max_iter, type)

      #------------------------
      # Core algorithm end
      #------------------------
      iter_vec[ig] <- fit$iter
      list_A[[ig]] <- t(fit$A)
      list_V[[ig]] <- fit$V
      list_Lambda[[ig]] <- fit$Lambda
    }
    cvxclust_obj <- list(U=list_A, V = list_V, Lambda=list_Lambda,nGamma=ig,iters=iter_vec,call=call, method = method)
    class(cvxclust_obj) <- "cvxclustobject"
    return(cvxclust_obj)
  }

  ###--------------------------
  ### ADMM part
  ###--------------------------
  if( method == "admm"){
    if( is.null(init) ){
      n <- nrow(X)
      p <- ncol(X)
      nK <- n*(n-1)/2
      A <- X
      Lambda <- matrix(0,p,nK)
      V <- matrix(0,p,nK)
    }else{
      Lambda <- init$lambda
      V <- init$V
      A <- init$A
    }
    # Group Lasso function
    foo <- function(y, N, Gamma2, Gamma2.weight){
      fit0 <- gglasso(x = N, y = y, group = rep(1,n), lambda = Gamma2 * Gamma2.weight )
      fit0$beta
    }

    for (ig in 1:nGamma) {
      gamma1 <- Gamma1[ig]
      fit <- admm_eigen(X, w, gamma1, gamma2, Gamma2_weight, nu, tol_abs, max_iter, type, foo)

      #------------------------
      # Core algorithm end
      #------------------------
      iter_vec[ig] <- fit$iter
      list_A[[ig]] <- fit$A
      list_V[[ig]] <- fit$V
      list_Lambda[[ig]] <- fit$Lambda
    }
    cvxclust_obj <- list(U=list_A, V = list_V, Lambda=list_Lambda,nGamma=ig,iters=iter_vec,call=call, eva = fit$eva, method = method)
    class(cvxclust_obj) <- "cvxclustobject"
    return(cvxclust_obj)
  }
}
