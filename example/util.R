###-----------Simulation Functions------------------###

#' Simulate data (4group)
#' 
#' Simulation 1: n=80, p = 50, 200, 500, 1000 (informative features 1:50), true K=4. Section 7.1, Scenario I in SWF12 paper.
#' @param n sample size
#' @param p dimension
#' @param k number of clusters of samples
#' @param mu mean value for generating mean vectors
#' @param sigma standard deviation 
#' @param seed for random data generation
#' 
#' @examples
#' DATA1 = simu_4(n = 80, p = 50, k = 4, mu = 0.4, sigma = 1)
#' image(t(DATA1),main="X")
#' kmeans(DATA1,4)$cluster 
simu_4 = function(n, true_p, p, k, mu, sigma, seed = NULL){
  
  if(! is.null(seed) ){ set.seed(seed)  }
  n_p <- true_p / 2
  n_k <- n / k
  clust.ind <- rep(1:k, each = n_k)
  clust.mat <- rbind( c(rep( mu, n_p), rep(-mu, n_p)), 
                      c(rep(-mu, n_p), rep(-mu, n_p)),
                      c(rep(-mu, n_p), rep( mu, n_p)),
                      c(rep( mu, n_p), rep( mu, n_p))
  )
  
  X = matrix(0,n,p)
  for(i in 1:n){
    mu_mean <- c( clust.mat[clust.ind[i],], rep(0, p - true_p) )
    X[i,] <- rnorm(p, mu_mean,rep(sigma,p))
  }
  
  list(X = X, label = clust.ind, features = c(rep(TRUE, true_p), rep(FALSE, p - true_p)))
  
}

### Prediction results for convex clustering
fit_cvxclust <- function(X, gamma1, k_w, phi, method, data_valide, type = 2, dist.type = "euclidean", verbose = TRUE ){
  fit <- list()
  n <- ncol(X)
  p <- nrow(X)
  
  #   dist.type <- c("euclidean", "maximum", "manhattan", "canberra", "binary" , "minkowski")
  w <- dist_weight(X / sqrt(p),phi, dist.type = dist.type, p = 2 )
  if(verbose)  hist(w, breaks = 20, main = "Histogram of weight")
  w <- knn_weights(w,k = k_w,n)
  
  nu <- AMA_step_size(w,n)/2  # Largest nu we can have
  
  
  fit <- cvxclust(X = X, w = w, gamma = gamma1, method = method, nu = nu, type = type)
  
  label_est <- list()
  features_est <- list()
  rand_est <- vector()
  for(i in 1:length(gamma1)){
    A <- create_adjacency( round(fit$V[[i]],2), w, n, method = method)
    label_est[[i]] <- find_clusters(A)
    features_est[[i]] <- (! apply(fit$U[[i]], 1, sd) == 0)
    rand_est[i] <- adjustedRandIndex( label_est[[i]]$cluster, data$label)
  }
  
  fit$cluster  <- lapply(label_est, function(x) x$cluster)
  fit$size     <- lapply(label_est, function(x) x$size)
  fit$rand     <- rand_est
  fit$gamma    <- gamma1
  fit$features <- rep(TRUE, p)
  fit$X        <- X
  
  fit$predict <- list()
  for(iter in 1:length(gamma1) ){
    cluster0 = fit$cluster[[iter]]
    features0 = fit$features
    pred = predict_validate(X, cluster = cluster0, features = features0, data_valide = data_valide)
    fit$predict[[iter]] <- c(gamma1 = gamma1[iter], gamma2 = NA, pred)
  }
  fit$predict <- do.call(rbind, fit$predict)
  
  fit
}

### Prediction results for sparse convex clustering
fit_sparse <- function(X, gamma1, gamma2, Gamma2_weight, k_w, phi, method = "ama", data_valide, type = 2, dist.type = "euclidean", verbose = TRUE ){
  
  n <- nrow(X)
  p <- ncol(X)
  
  #   dist.type <- c("euclidean", "maximum", "manhattan", "canberra", "binary" , "minkowski")
  w <- dist_weight( t(X) / sqrt(p),phi, dist.type = dist.type, p = 2 )
  if(verbose)  hist(w, breaks = 20, main = "Histogram of weight")
  w <- knn_weights(w,k = k_w,n)
  w <- w / sum(w) * sqrt(p) * sqrt(n)   # Standardize Lambda1 weight
  
  nu <- AMA_step_size(w,n)/2  # Largest nu we can have
  
  
  fit <- scvxclustr::scvxclust(X = X, w = w, Gamma1 = gamma1, Gamma2 = gamma2, Gamma2_weight = Gamma2.weight, nu = nu, method = method, tol_abs = 1e-7)
  
  label_est <- list()
  features_est <- list()
  rand_est <- vector()
  for(i in 1:length(gamma1)){
    A <- create_adjacency( round(fit$V[[i]],2), w, n, method = method)
    label_est[[i]] <- find_clusters(A)
    features_est[[i]] <- (! apply(fit$U[[i]], 2, sd) == 0)
    rand_est[i] <- adjustedRandIndex( label_est[[i]]$cluster, data$label)
  }
  
  fit$cluster  <- lapply(label_est, function(x) x$cluster)
  fit$size     <- lapply(label_est, function(x) x$size)
  fit$rand     <- rand_est
  fit$gamma    <- gamma1
  fit$features <- features_est
  fit$X        <- X
  
  
  cluster0 = fit$cluster[[1]]
  features0 = fit$features[[1]]
  pred = predict_validate( t(X), cluster = cluster0, features = features0, data_valide)
  fit$predict = c(gamma1 = gamma1, gamma2 = gamma2, pred)
  
  fit
}