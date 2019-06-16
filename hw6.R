library(MASS)
set.seed(521)

X <<- NULL
centroids <<- NULL
assignments <<- NULL

X1 <- mvrnorm(50, c(+2.5, +2.5), matrix(c(+0.8, -0.6, -0.6, +0.8), 2, 2))
X2 <- mvrnorm(50, c(-2.5, +2.5), matrix(c(+0.8, +0.6, +0.6, +0.8), 2, 2))
X3 <- mvrnorm(50, c(-2.5, -2.5), matrix(c(+0.8, -0.6, -0.6, +0.8), 2, 2))
X4 <- mvrnorm(50, c(+2.5, -2.5), matrix(c(+0.8, +0.6, +0.6, +0.8), 2, 2))
X5 <- mvrnorm(100, c(0, 0), matrix(c(+1.6, 0, 0, +1.6), 2, 2))

X <- rbind(X1, X2, X3, X4, X5)

K <- 5
N <- length(X[,1])


centroids <- X[sample(1:N, K),]

k_iteration <- 2

for(i in 1:k_iteration){
  
  distances <- as.matrix(dist(rbind(centroids, X), method = "euclidean"))
  distances <- distances[1:nrow(centroids), (nrow(centroids) + 1):(nrow(centroids) + nrow(X))]
  assignments <- sapply(1:ncol(distances), function(i) {which.min(distances[,i])})
  
  centroids <- t(sapply(1:K, function(k){colMeans(rbind(X[assignments == k,]))}))
}

H <- matrix(sapply(assignments, function(component){ (1:K) == component }), N, K, byrow = TRUE)

iteration <- 100
i <- 0
D <- 2

gaussian_helper <- function(X, cov){X %*% solve(cov) %*% cbind(X)}

gaussian_density_func <- function(x_row, k){
  component_priors[k]*(det(covariances[,,k])**-0.5) * exp(-0.5*gaussian_helper(x_row - centroids[k,], covariances[,,k]))
}
multiplier_helper <- function(side, middle){ side %*% middle %*% t(side) }

while(i<iteration){
  

  

  covariances <- sapply(X = 1:K, FUN = function(k) {   (t(X) - matrix(centroids[k,], D, N)) %*% diag(H[,k]) %*% t(t(X) - matrix(centroids[k,], D, N))/ sum(H[,k]) }, simplify = "array")

  component_priors <- colMeans(H)
  
  H <- t(sapply(1:N, function(n){
    row <- sapply(1:K, function(k){gaussian_density_func(X[n,], k)})
    return(row / sum(row))
  }))
  
  centroids <- (t(H) %*% X ) / matrix(colSums(H), K, D)
  i <- i + 1
}



print(centroids)

assignments <- max.col(H)

colors <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#b15928", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99")



plot(X[,1], X[,2], type = "p", pch = 19, col = colors[assignments], las = 1,
     xlim = c(-6, 6), ylim = c(-6, 6),
     xlab = "x1", ylab = "x2")


x1_interval <- seq(from = -6, to = +6, by = 0.12)
x2_interval <- seq(from = -6, to = +6, by = 0.12)
x1_grid <- matrix(x1_interval, nrow = length(x1_interval), ncol = length(x1_interval), byrow = FALSE)
x2_grid <- matrix(x2_interval, nrow = length(x2_interval), ncol = length(x2_interval), byrow = TRUE)

f <-  function(k){ function(x1, x2){ gaussian_density_func(c(x1,x2), k)}}


learned_component_priors <- component_priors
learned_component_covariances <- covariances
learned_centroids <- centroids

#for plotting
class_means <- matrix(c(+2.5, +2.5,
                        -2.5, +2.5,
                        -2.5, -2.5,
                        +2.5, -2.5,
                        +0.0, +0.0), 2, 5)

real_component_covariances <- array(c(+0.8, -0.6, -0.6, +0.8,
                                      +0.8, +0.6, +0.6, +0.8,
                                      +0.8, -0.6, -0.6, +0.8,
                                      +0.8, +0.6, +0.6, +0.8,
                                      +1.6, +0.0, +0.0, +1.6), c(2, 2, 5))

sample_sizes <- c(rep(50, 4), 100)

for(i in 1:2){
  density_values <- lapply(1:K, function(k){matrix(mapply(f(k), x1_grid, x2_grid), nrow(x2_grid), ncol(x2_grid))})
  
  for(k in 1:K){
    contour(x1_interval, x2_interval, density_values[[k]], lty = i, levels = c(0.05), add = TRUE, lwd = 2, drawlabels = FALSE)
  }

  component_priors <- sample_sizes / sum(sample_sizes)
  covariances <- real_component_covariances
  centroids <- t(class_means)
  
  
}

