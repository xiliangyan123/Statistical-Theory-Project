##################################
# Fitting A Correct Model 
# Lognormal Data / Lognormal Model 
##################################
## only estimating the mean 

## defining parameters
set.seed(53)
n <- c(600, 400, 200)

## Generate data based on different sample sizes
X1 <- rlnorm(n[1], meanlog=1, sdlog=1) 
X2 <- rlnorm(n[2], meanlog=1, sdlog=1)
X3 <- rlnorm(n[3], meanlog=1, sdlog=1)

# grid1 <- seq(0,29.95/3,.05/3)
# grid2 <- seq(0,19.95/2,.05/2)
# grid3 <- seq(0,9.95,.05)

#X should be concentrating around 4.48. 
# Not sure if we need this, since we have it below.(In histogram)
# dpost1 <- plot(grid1, dlnorm(grid1, 1, 1), type = "l", xlab = "Mean Parameter Estimate", ylab = "f(x)",
#                main = "Parameter Estimation under Bayesian Framework")
# lines(density(X1), col = "Red", lty = 2, lwd = 5)
# points(4.48, 0, pch = "X")
# dpost2 <- plot(grid2, dlnorm(grid2, 1, 1), type = "l", xlab = "Mean Parameter Estimate", ylab = "f(x)",
#                main = "Parameter Estimation under Bayesian Framework")
# lines(density(X2), col = "Green", lty = 3, lwd = 5)
# points(4.48, 0, pch = "X")
# dpost3 <- plot(grid3, dlnorm(grid3, 1, 1), type = "l", xlab = "Mean Parameter Estimate", ylab = "f(x)",
#                main = "Parameter Estimation under Bayesian Framework")
# lines(density(X3), col = "Blue", lty = 4, lwd = 5)
# points(4.48, 0, pch = "X")

## Trying to calculate VaR and CTE
VaR.hat <- qlnorm(0.95, meanlog=1, sdlog=1)
d <- function(x){exp(-(log(x)-1)^2/2)/(sqrt(2*pi))}
CTE.hat <- integrate(d, lower=VaR.hat, upper=Inf)$value/0.05

N <- 10000
keepers <- matrix(0, N, 3)
keepers2 <- matrix(0, N, 3)
keepers3 <- matrix(0, N, 3)
#mean

#we know arithmetic mean and variance. 
#result found on wikipedia.(Lognormal Distribution)
#want to use this to get the meanlog = 1 and sdlog = 1 for VaR.star
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

## Plots for sample size of 600
plot(1:N, keepers[,1], type="l", main="Trace plot for Credibility Estimator")
hist(keepers[,1], main="Plot of the predicted mean")
points(mean(X1), 0, pch = "X") #predictive mean

plot(1:N, keepers[,2], type="l", main="Trace plot for VaR")
hist(keepers[,2], main="Plot of the VaR estimate")
points(VaR.hat, 0, pch = "X") #predicted VaR estimate

plot(1:N, keepers[,3], type="l", main="Trace plot for CTE")
hist(keepers[,3], main="Plot of the CTE estimate")
points(CTE.hat, 0, pch = "X") #predicted CTE estimate

## Plots for sample size of 400
plot(1:N, keepers2[,1], type="l", main="Trace plot for Credibility Estimator")
hist(keepers2[,1], main="Plot of the predicted mean")
points(mean(X2), 0, pch = "X") #predictive mean

plot(1:N, keepers2[,2], type="l", main="Trace plot for VaR")
hist(keepers2[,2], main="Plot of the VaR estimate")
points(VaR.hat, 0, pch = "X") #predicted VaR estimate

plot(1:N, keepers2[,3], type="l", main="Trace plot for CTE")
hist(keepers2[,3], main="Plot of the CTE estimate")
points(CTE.hat, 0, pch = "X") #predicted CTE estimate

## Plots for sample size of 200
plot(1:N, keepers3[,1], type="l", main="Trace plot for Credibility Estimator")
hist(keepers3[,1], main="Plot of the predicted mean")
points(mean(X3), 0, pch = "X") #predictive mean

plot(1:N, keepers3[,2], type="l", main="Trace plot for VaR")
hist(keepers3[,2], main="Plot of the VaR estimate")
points(VaR.hat, 0, pch = "X") #predicted VaR estimate

plot(1:N, keepers3[,3], type="l", main="Trace plot for CTE")
hist(keepers3[,3], main="Plot of the CTE estimate")
points(CTE.hat, 0, pch = "X") #predicted CTE estimate

