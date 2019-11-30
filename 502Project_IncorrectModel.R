##-----------------------------------------------------------------------------------
## TITLE:     Model correctlt specified, Bayesian versus credibility estimation,
##            and Gibbs posterior
##
## AUTHORS:    L. Hong and R. Martin
##-----------------------------------------------------------------------------------

set.seed(53)
alpha <- 3 
lambda <- 0.2
n <- c(600, 400, 200)

X <- rlnorm(n[1], meanlog=1, sdlog=1) #generate random number log normal distribution
X2 <- rlnorm(n[2], meanlog=1, sdlog=1)
X3 <- rlnorm(n[3], meanlog=1, sdlog=1)

#Define true parameter value
VaR.hat <- qlnorm(0.95, meanlog=1, sdlog=1)
d <- function(x){exp(-(log(x)-1)^2/2)/sqrt(2*pi)} # Log normal distribution with var = 1 and mean = 1. 
CTE.hat <- integrate(d, lower=VaR.hat, upper=Inf)$value/0.05 #lower tail expectation. 

#alpha is shape and lambda is 1/rate 
#simulating Gibbs

#---------------------------------------------------
# plots of the posterior for different sample sizes
#---------------------------------------------------
#Simulation for one trial with n=600

dpost1 <- function(x){dnorm(x, mean=mean(X), sd=sqrt(n[1]/((n[1]/(var(X)) + 1))))}
curve(dpost1, from=4, to=5, xlab="parameter", ylab="posterior density", lty=1)
points(mean(X), pch = "X") #True value of 

N = 10000
#Record
keepers <- matrix(0, N, 3)
sigma.hat2 <- n[1]/((n[1]/(var(X)) + 1))

Y <- dnorm(1000, mean=mean(X), sd=sqrt(sigma.hat2))

for (i in 1:N) {
  Y <- rnorm(n[1], mean=mean(X), sd=sqrt(sigma.hat2))
  mu <- mean(Y)
  VaR.star <- qnorm(0.95, mean = mean(Y), sd = sqrt(var(Y)))
  CTE.star <- integrate(d, lower=VaR.star, upper=Inf)$value/0.05
  keepers[i,] <- c(mu,VaR.star,CTE.star)
}

plot(1:N, keepers[,1], type="l", main="Trace plot for Credibility Estimator")
hist(keepers[,1], main="Plot of the predicted mean")
points(mean(X), 0, pch = "X") #predictive mean

plot(1:N, keepers[,2], type="l", main="Trace plot for VaR")
hist(keepers[,2], main="Plot of the VaR estimate")
points(VaR.hat, 0, pch = "X") #predicted VaR estimate

plot(1:N, keepers[,3], type="l", main="Trace plot for CTE")
hist(keepers[,3], main="Plot of the CTE estimate")
points(CTE.hat, 0, pch = "X") #predicted CTE estimate