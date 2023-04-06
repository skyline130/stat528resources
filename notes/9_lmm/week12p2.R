
rm(list = ls())

#####################
## Generate data
#####################

set.seed(13)

# intercept
m <- 2

# random intercepts
a <- -3:3

# sample sizes
n <- length(a)
T <- rpois(n, lambda = 20)

# response vector
y <- unlist(lapply(1:n , function(i){
  rnorm(T[i], m + a[i], 1)
}))



#########################################
## initial values and hypperparameters
#########################################

## initial values
mu <- 0
aest <- rep(0,n)
sigma_a <- sigma_mu <- sigma_e <- 0.5

## hyperparameters
alpha_a <- alpha_e <- 2
beta_a <- beta_e <- 100

## initialize matrix that will store MCMC output
parms <- c(sigma_a, sigma_e, mu, aest)
mat <- matrix(0, ncol = length(parms), nrow = 5e4)
mat[1, ] <- parms

## get indices of responses across the groups 
## (for organization)
index <- lapply(1:n, function(j){
  ( (1 + c(0, cumsum(T)[-n]))[j]):cumsum(T)[j]
})


########################
## Run Gibbs sampler
########################

library(MCMCpack)
system.time({for(j in 2:5e4){
  # generate sigma_a
  mat[j, 1] <- rinvgamma(1, alpha_a + n/2, 
                         beta_a + 0.5 * crossprod(mat[j-1, 4:10]))
  
  # generate sigma_e
  mat[j, 2] <- rinvgamma(1, alpha_e + 0.5 * sum(T), 
                         beta_e + 0.5 * crossprod(y - mat[j-1, 3] - 
                                                    rep(mat[j-1,4:10], times = T)))
  
  # generate intercept
  mat[j, 3]  <- rnorm(1, 
                      mean = sum(y -  rep(mat[j-1,4:10], times = T))/
                        (sum(T)/sigma_e + 1/sigma_mu), 
                      sd = 1/sqrt(sum(T)/sigma_e + 1/sigma_mu))
  
  # generate all a terms in parallel
  mat[j, 4:10] <- rnorm(n, 
                        mean = unlist(lapply(index, function(k) sum(y[k] - mat[j, 3])))/
                          (T + sigma_e/sigma_a), 
                        sd = 1/sqrt(T/sigma_e + 1/sigma_a))
}})



#################
## Diagnostics
#################

# We now check some diagnostics. The first of these are acf plots 
# which check how much dependence is present among realizations in 
# our Gibbs sampler. The top left panel corresponds to the global 
# intercept. The other three plots correspond to the first three 
# subject specific intercept terms. These plots indicate that not 
# much dependency exists in our Gibbs sampler.

par(mfrow = c(2,2))
for(j in 3:6) acf(mat[, j])

# The next four plots correspond to the last four subject specific 
# intercept terms. These plots indicate that not much dependency 
# exists in our Gibbs sampler.

par(mfrow = c(2,2))
for(j in 7:10) acf(mat[, j])

# Next we check mixing issues of our Gibbs sampler using traceplots. 
# A chain that mixes well has a traceplot that is bouncy and appears 
# to not have a trend. The top left panel corresponds to the global 
# intercept. The other three plots correspond to the first three 
# subject specific intercept terms. These plots indicate that our 
# Gibbs sampler is mixing well. 

# takes a while to load
#par(mfrow = c(2,2))
#for(j in 3:6){ 
#  plot(mat[, j])
#  lines(mat[, j])
#}

# takes a while to load
#par(mfrow = c(2,2))
#for(j in 7:10){ 
#  plot(mat[, j])
#  lines(mat[, j])
#}

# We now check the posterior distributions for the global intercept 
# and the subject-specific intercept terms. Note that our setup 
# specified $\mu = 2$ and $a = -3,-2,\ldots, 2,3$, and our posterior 
# mode for the global intercept is roughly $\hat\mu = 0$ and our 
# posterior modes estimate $\hat a$ closer to $-1,0,\ldots,4,5$.

par(mfrow = c(2,2))
plot(density(mat[, 3]), main = "posterior density of global intercept")
plot(density(mat[, 4]), main = "posterior density of a_1")
plot(density(mat[, 5]), main = "posterior density of a_2")
plot(density(mat[, 6]), main = "posterior density of a_3")

par(mfrow = c(2,2))
plot(density(mat[, 7]), main = "posterior density of a_4")
plot(density(mat[, 8]), main = "posterior density of a_5")
plot(density(mat[, 9]), main = "posterior density of a_6")
plot(density(mat[, 10]), main = "posterior density of a_7")


