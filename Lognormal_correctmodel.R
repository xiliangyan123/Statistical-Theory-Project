##################################
# Fitting A Correct Model 
# Lognormal Data / Lognormal Model 
##################################
## only estimating the mean 

## defining parameters
library(ggplot2)
set.seed(53)
n <- c(600, 400, 200)

## Generate data based on different sample sizes
X1 <- rlnorm(n[1], meanlog=1, sdlog=1) 
X2 <- rlnorm(n[2], meanlog=1, sdlog=1)
X3 <- rlnorm(n[3], meanlog=1, sdlog=1)

## Trying to calculate VaR and CTE
## Full Bayesian Analysis - True Values
VaR.hat <- qlnorm(0.95, meanlog=1, sdlog=1)
d <- function(x){exp(-(log(x)-1)^2/2)/(sqrt(2*pi))}
CTE.hat <- integrate(d, lower=VaR.hat, upper=Inf)$value/0.05

## Number of iterations and vectors to hold data. 
N <- 10000
keepers <- matrix(0, N, 3)
keepers2 <- matrix(0, N, 3)
keepers3 <- matrix(0, N, 3)
#mean

## result found on wikipedia.(Lognormal Distribution)
## want to use this to get the meanlog = 1 and sdlog = 1 for VaR.star
location <- log(1/sqrt(1+1))
shape <- sqrt(log(1+1))

for (i in 1:N) {
  Y <- rlnorm(n[1], meanlog=1, sdlog=1)
  mu <- mean(Y)
  Z <- rlnorm(n[1], location, shape)
  VaR.star <- qlnorm(0.95, meanlog = mean(Z), sdlog = sqrt(var(Z)))
  CTE.star <- integrate(d, lower=VaR.star, upper=Inf)$value/0.05
  keepers[i,] <- c(mu,VaR.star,CTE.star)
}

for (j in 1:N) {
  Y2 <- rlnorm(n[2], meanlog=1, sdlog=1)
  mu2 <- mean(Y2)
  Z2 <- rlnorm(n[2], location, shape)
  VaR.star2 <- qlnorm(0.95, meanlog = mean(Z2), sdlog = sqrt(var(Z2)))
  CTE.star2 <- integrate(d, lower=VaR.star2, upper=Inf)$value/0.05
  keepers2[j,] <- c(mu2,VaR.star2,CTE.star2)
}

for (k in 1:N) {
  Y3 <- rlnorm(n[3], meanlog=1, sdlog=1)
  mu3 <- mean(Y3)
  Z3 <- rlnorm(n[3], location, shape)
  VaR.star3 <- qlnorm(0.95, meanlog = mean(Z3), sdlog = sqrt(var(Z3)))
  CTE.star3 <- integrate(d, lower=VaR.star3, upper=Inf)$value/0.05
  keepers3[k,] <- c(mu3,VaR.star3,CTE.star3)
}

#Plots of the posterior means
plot(density(keepers[,1], from = 4, to = 6), 
     main="Mean Estimates", xlab = "Mu", ylim = c(0, 2), 
     lwd = 5, col = "Green")
par(new=T)
plot(density(keepers2[,1], from = 4, to = 6), 
     main = "", xlab = "", ylim = c(0,2), 
     lty = 3, lwd = 5, col = "Red")
par(new=T)
plot(density(keepers3[,1], from = 4, to = 6), 
     main = "", xlab = "", ylim = c(0,2), 
     lty = 2, lwd = 5, col = "Gold")
par(new=T)
points(4.48, 0, pch = "X")

#Plots of the VaR for different sample sizes
plot(density(keepers[,2], from = 0, to = 100), 
     main="Var Estimates", xlab = "VaR", ylim = c(0,0.2), 
     lwd = 5, col = "Green")
par(new=T)
plot(density(keepers2[,2], from = 0, to = 100), 
     main = "", xlab = "", ylim = c(0,0.2), 
     lty = 3, lwd = 5, col = "Red")
par(new=T)
plot(density(keepers3[,2], from = 0, to = 100), 
     main = "", xlab = "", ylim = c(0,0.2), 
     lty = 2, lwd = 5, col = "Gold")
par(new=T)
points(VaR.hat, 0, pch = "X")

#Plots of the CTE for different sample sizes. 
plot(density(keepers[,3], from = 0, to = 100), 
     main="CTE Estimates", xlab = "Mu", ylim = c(0, 0.10), 
     lwd = 5, col = "Green")
par(new=T)
plot(density(keepers2[,3], from = 0, to = 100), 
     main = "", xlab = "", ylim = c(0,0.1), 
     lty = 3, lwd = 5, col = "Red")
par(new=T)
plot(density(keepers3[,3], from = 0, to = 100), 
     main = "", xlab = "", ylim = c(0,0.1), 
     lty = 2, lwd = 5, col = "Gold")
par(new=T)
points(CTE.hat, 0, pch = "X")

# ## Plots for sample size of 600
# plot(1:N, keepers[,1], type="l", main="Trace plot for Credibility Estimator", ylab = "Mu")
# hist(keepers[,1], main="Mean Estimate with n = 600", xlab = "Mu", prob = TRUE) 
# points(mean(X1), 0, pch = "X") #predictive mean
## Histogram Plots
# plot(1:N, keepers[,2], type="l", main="Trace plot for VaR", ylab = "VaR")
# hist(keepers[,2], main="VaR estimate with n = 600", xlab = "VaR", prob = TRUE)
# points(VaR.hat, 0, pch = "X") #predicted VaR estimate
# 
# plot(1:N, keepers[,3], type="l", main="Trace plot for CTE", ylab = "CTE")
# hist(keepers[,3], main="CTE estimate with n = 600", xlab = "CTE", prob = TRUE)
# points(CTE.hat, 0, pch = "X") #predicted CTE estimate
# 
# ## Plots for sample size of 400
# plot(1:N, keepers2[,1], type="l", main="Trace plot for Credibility Estimator", ylab = "Mu")
# hist(keepers2[,1], main="Mean Estimate with n = 400", xlab = "Mu", prob = TRUE)
# points(mean(X2), 0, pch = "X") #predictive mean
# 
# plot(1:N, keepers2[,2], type="l", main="Trace plot for VaR", ylab = "VaR")
# hist(keepers2[,2], main="VaR estimate with n = 400", xlab = "VaR", prob = TRUE)
# points(VaR.hat, 0, pch = "X") #predicted VaR estimate
# 
# plot(1:N, keepers2[,3], type="l", main="Trace plot for CTE", ylab = "CTE")
# hist(keepers2[,3], main="CTE estimate with n = 400", xlab = "CTE", prob = TRUE)
# points(CTE.hat, 0, pch = "X") #predicted CTE estimate
# points(mean(CTE.star), 0, pch = "O")
# 
# ## Plots for sample size of 200
# plot(1:N, keepers3[,1], type="l", main="Trace plot for Credibility Estimator", ylab = "Mu")
# hist(keepers3[,1], main="Mean Estimate with n = 200", xlab = "Mu", prob = TRUE)
# points(mean(X3), 0, pch = "X") #predictive mean
# 
# plot(1:1000, keepers3[c(1:1000),2], type="l", main="Trace plot for VaR", ylab = "VaR")
# hist(keepers3[c(1:1000),2], main="VaR estimate with n = 200", xlab = "VaR", prob = TRUE)
# points(VaR.hat, 0, pch = "X") #predicted VaR estimate
# 
# plot(1:N, keepers3[,3], type="l", main="Trace plot for CTE", ylab = "CTE")
# hist(keepers3[,3], main="CTE estimate with n = 200", xlab = "CTE", prob = TRUE)
# points(CTE.hat, 0, pch = "X") #predicted CTE estimate
