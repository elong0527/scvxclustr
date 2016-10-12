#' Weight function for w
#'
#' wrap up of the stat::dist function for the weight of data
#' @param X data matrix p by n
#' @param phi scale
#' @param dist.type c("euclidean", "maximum", "manhattan", "canberra", "binary" , "minkowski")
#' @param p power of distance
#'
#' @export
dist_weight <- function(X, phi, dist.type, p){
  dist_X <- as.numeric(dist(t(X), method = dist.type))
  exp(- phi * dist_X^p  )
}
