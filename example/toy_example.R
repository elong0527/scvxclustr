
set.seed(123)

# Load scvxclustr package
library(devtools)
library(scvxclustr)

# load necessary functions
source("util.R")


n = 60          # Sample size
true_p = 20     # Number of true features
p = 80        # Number of total features (p = 150, 500)
k = 4           # Number of true cluster
mu = 1.2        # mean of normal distribution (mu = 0.6, 0.9)
sigma = 1       # sd of normal distribution
method = "ama"  # Fitted method (method = "ama", "admm")

# Simiulate 4 cluster Gaussian data
data <- simu_4(n = n, true_p = true_p, p = p, k = k, mu = mu, sigma = sigma )

#standardize n by p data matrix
X <- scale(data$X,center=TRUE,scale=FALSE)

# Validation data
data_valide <- list()
for(i in 1:5){
  data_valide[[i]] <- simu_4(n = n, true_p = true_p, p = p, k = k, mu = mu, sigma = sigma )
}

# Adaptive Weight (if possible)
g1 <-6
g2 <- 0
Gamma2.weight <- c(rep(0.5, true_p), rep(1,p - true_p) )
k_w <- 5    # Number of nearest neighbors
phi <- 0.5  # scale of the kernel
verbose <- TRUE # show more information
w <- dist_weight( t(X) / sqrt(p),phi, dist.type = "euclidean", p = 2 )
w <- knn_weights(w,k = k_w,n)
nu <- AMA_step_size(w,n) /2


## Validate the cvxclust and scvxclust is the same when g2 = 0
# Fit a convex clustering model
fit1 <- cvxclust(X = t(X), w = w, gamma = g1, method = "ama", nu = nu, max_iter = 10000, tol = 1e-5)

# Fit a sparce convex clsutering model
fit2 <- scvxclust(X = X, w = w, Gamma1 = g1, Gamma2 = g2, Gamma2_weight = Gamma2.weight, method = method, nu = nu, max_iter = 10000, tol_abs = 1e-5)

fit1$iters
fit2$iters
diff_U <- as.numeric(fit1$U[[1]] - t(fit2$U[[1]]) )
summary( diff_U )
plot(diff_U)

## Validate the sparse convex clustring model create a correct clustering structure under the tody exmaple.
g1 <- 9
g2 <- 10
fit_predict <- fit_sparse(X = X, gamma1 = g1, gamma2 = g2, Gamma2.weight, k_w, phi, method = method, data_valide = data_valide, verbose = F)
fit_predict$predict
table(data$label, fit_predict$cluster[[1]])
