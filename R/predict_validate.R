#' Predict validation data
#'
#' report Rand index and variable stability in validate dataset
#'
#' @param X centralized data (p by n)
#' @param cluster cluster assignment for each subject
#' @param features used features
#' @param data_valide validate dataset list
#'
#' @export
predict_validate <- function(X, cluster, features, data_valide){

  size <- table(cluster)
  n.cluster <- length(size)

  Rand.valide <- vector()
  for(i.rand in 1:length(data_valide)){

    data_valide0 <- data_valide[[i.rand]]
    X_valide <- t(scale(data_valide0$X[, features],center=TRUE,scale=FALSE))
    p1  <- nrow(X_valide)
    n1  <- ncol(X_valide)

    center <- matrix(0, n.cluster, p1)
    dist <- matrix(0, n.cluster, n1)
    for(i.cluster in 1:n.cluster){
      center0 <- apply( as.matrix( X[features, cluster == i.cluster]), 1, mean)
      dist0 <- apply(X_valide, 2, function(x) sqrt( sum( (x - center0)^2 ) ) )
      center[i.cluster, ] <- center0
      dist[i.cluster, ] <- dist0
    }
    cluster.valide <- apply(dist, 2, which.min)
    Rand.valide[i.rand] <- adjustedRandIndex(cluster.valide, data_valide0$label)
  }
  Rand.valide <- mean(Rand.valide)
  #   features.stab <- var.stab(which(features), which(data$features), p = p1)

  # MCC
  true_p <- table(data_valide[[1]]$features)
  true_p <- rep(true_p, each = 2)
  features0 <- factor(features + 0, levels = c(0, 1))
  tab = table(features0, data_valide0$features)
  tmp = as.numeric(tab) / true_p
  names(tmp) <- c("trueF","typeI","typeII","trueT")

  tp  <- tab[2, 2]
  fp  <- tab[1, 2]
  fn  <- tab[2, 1]
  tn  <- tab[1, 1]
  mcc = MCC(tp, tn, fp, fn)

  c(Rand = Rand.valide, tmp, mcc = mcc)
}

MCC <- function(tp, tn, fp, fn) {
  # compute Matthews correlation coefficient
  prod <- log2(tp + fp) + log2(tp + fn) + log2(tn + fp) + log2(tn + fn)
  return((tp * tn - fp * fn) / sqrt(2 ^ prod))
}
