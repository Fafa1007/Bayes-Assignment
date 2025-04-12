### Q3.b ###

library(mvtnorm)

gibbs_sampler <- function(Y, X, n1 # Number of standard observations
                          , S = 5000 # Total number of Gibbs iterations
                          , burn = 1000 # Number of burn-in samples to discard
                          , a = 1, b = 1, T0 = diag(ncol(X)) # Prior covariance matrix
                          ) {
  
  n <- nrow(X)
  k <- ncol(X)
  n2 <- n - n1
  
  # Split data 
  # First n1 observations treated as regular
  # Remaning n2 = n - n1 are considered potential outliers
  X1 <- X[1:n1, ]
  Y1 <- Y[1:n1]
  X2 <- X[(n1 + 1):n, ]
  Y2 <- Y[(n1 + 1):n]
  
  # Storage
  beta_samples <- matrix(0, S, k)
  tau1_samples <- numeric(S)
  tau2_samples <- numeric(S)
  
  # Initialization
  beta <- rep(0, k)
  tau1 <- 1
  tau2 <- 0.5
  
  for (s in 1:S) {
    ## 1. Sample beta 
    # Beta|tau1,tau2,X,Y
    Sigma_beta_inv <- solve(T0) + tau1 * crossprod(X1) + tau2 * crossprod(X2)
    Sigma_beta <- solve(Sigma_beta_inv)
    mu_beta <- Sigma_beta %*% (tau1 * t(X1) %*% Y1 + tau2 * t(X2) %*% Y2)
    beta <- as.vector(mvtnorm::rmvnorm(1, mean = mu_beta, sigma = Sigma_beta))
    
    ## 2. Sample tau1
    #tau1|Beta,X,Y
    res1 <- Y1 - X1 %*% beta
    shape1 <- a + n1 / 2
    rate1 <- b + 0.5 * sum(res1^2)
    tau1 <- rgamma(1, shape = shape1, rate = rate1)
    
    ## 3. Sample tau2 (truncated gamma: tau2 < tau1)
    #tau2|B,tau1,X,Y
    res2 <- Y2 - X2 %*% beta
    shape2 <- a + n2 / 2
    rate2 <- b + 0.5 * sum(res2^2)
    
    repeat {
      tau2_candidate <- rgamma(1, shape = shape2, rate = rate2)
      if (tau2_candidate < tau1) {
        tau2 <- tau2_candidate
        break
      }
    }
    
    # Store samples
    beta_samples[s, ] <- beta
    tau1_samples[s] <- tau1
    tau2_samples[s] <- tau2
  }
  
  # Return posterior draws (discard burn-in)
  list(
    beta = beta_samples[(burn+1):S,],
    tau1 = tau1_samples[(burn+1):S],
    tau2 = tau2_samples[(burn+1):S]
  )
}

#Packagesrequired
require(MASS)
require(cubature)

#Lets simulate some data
set.seed(2025)
n=150 # Number of data points
n1 = n/2
sd1 = 1.2
sd2 = 1.5

X1.c=data.frame(matrix(rnorm(n = 5 * n1),ncol=5))
X2.c=data.frame(matrix(rnorm(n = 5* (n-n1)),ncol=5))
X.c=rbind(X1.c, X2.c)
colnames(X.c)= c("X1","X2","X3","X4", "X5")
X=as.matrix(cbind(1,X.c)) #Designmatrix

e1=matrix(rnorm(n = n1 , sd = sd1^2),ncol=1) #Errors
e2=matrix(rnorm(n= n-n1, sd = sd2^2),ncol=1) #Errors
beta.true=matrix(c(1,0,10,0,2,-3),ncol=1)
Y1 = as.matrix(cbind(1,X1.c))%*%beta.true+e1
Y2 = as.matrix(cbind(1,X2.c))%*%beta.true+e2
Y = rbind(Y1,Y2)

gibbs_sampler(Y, X, n1)
