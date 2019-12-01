##-----------------------------------------------------------------------------------
## TITLE:     Model correctlt specified, Bayesian versus credibility estimation,
##            and Gibbs posterior
##
## AUTHORS:    L. Hong and R. Martin
##-----------------------------------------------------------------------------------

## defining our sample size
set.seed(100)
n <- c(600, 400, 200)

#generate random number log normal distribution, for different n. 
X1 <- rlnorm(n[1], meanlog=1, sdlog=1) 
X2 <- rlnorm(n[2], meanlog=1, sdlog=1)
X3 <- rlnorm(n[3], meanlog=1, sdlog=1)

#Define true parameter values
True.Mu <- exp(1 + 0.50)
VaR.True <- qlnorm(0.95, meanlog=1, sdlog=1)
d <- function(x){exp(-(log(x)-1)^2/2)/sqrt(2*pi)}
CTE.True <- integrate(d, lower=VaR.True, upper=Inf)$value/0.05

## Number of iterations
N = 10000

## Matrices to hold our parameter values
keepersA <- matrix(0, N, 3)
keepersB <- matrix(0, N, 3)
keepersC <- matrix(0, N, 3)

## variances of the gibbs posterior. 
sigma2.n1 <- n[1]/((n[1]/(var(X1)) + 1))
sigma2.n2 <- n[2]/((n[2]/(var(X2)) + 1))
sigma2.n3 <- n[3]/((n[3]/(var(X3)) + 1))

## Generating parameter values for n = 600. 
for (m in 1:N) {
  Y1 <- rnorm(n[1], mean=mean(X1), sd=sqrt(sigma2.n1))
  mu1 <- mean(Y1)
  VaR.est1 <- qnorm(0.95, mean = mean(Y1), sd = sqrt(var(Y1)))
  CTE.est1 <- integrate(d, lower=VaR.est1, upper=Inf)$value/0.05
  keepersA[m,] <- c(mu1,VaR.est1,CTE.est1)
}

## Generating parameter values for n = 400. 
for (q in 1:N) {
  Y2 <- rnorm(n[2], mean=mean(X2), sd=sqrt(sigma2.n2))
  mu2 <- mean(Y2)
  VaR.est2 <- qnorm(0.95, mean = mean(Y2), sd = sqrt(var(Y2)))
  CTE.est2 <- integrate(d, lower=VaR.est2, upper=Inf)$value/0.05
  keepersB[q,] <- c(mu2,VaR.est2,CTE.est2)
}

## Generating parameter values for n = 200. 
for (t in 1:N) {
  Y3 <- rnorm(n[3], mean=mean(X3), sd=sqrt(sigma2.n3))
  mu3 <- mean(Y3)
  VaR.est3 <- qnorm(0.95, mean = mean(Y3), sd = sqrt(var(Y3)))
  CTE.est3 <- integrate(d, lower=VaR.est3, upper=Inf)$value/0.05
  keepersC[t,] <- c(mu3,VaR.est3,CTE.est3)
}

#Plots of the posterior means
plot(density(keepersA[,1], from = 4, to = 6), 
     main="Mean Estimates", xlab = "Mu", ylim = c(0, 2), 
     lwd = 5, col = "Green")
par(new=T)
plot(density(keepersB[,1], from = 4, to = 6), 
     main = "", xlab = "", ylim = c(0,2), 
     lty = 3, lwd = 5, col = "Red")
par(new=T)
plot(density(keepersC[,1], from = 4, to = 6), 
     main = "", xlab = "", ylim = c(0,2), 
     lty = 2, lwd = 5, col = "Brown")
par(new=T)
points(True.Mu, 0, pch = "X")
legend("topright", title = "Sample Sizes", 
       legend =c("n=600", "n=400", "n=200"), 
       col = c("Green", "Red", "Brown"),
       bg="lightblue", lwd=5, lty=1:3)

#Plots of the VaR for different sample sizes
plot(density(keepersA[,2], from = 10, to = 20), 
     main="Var Estimates", xlab = "VaR", ylim = c(0,2.0), 
     lwd = 5, col = "Green")
par(new=T)
plot(density(keepersB[,2], from = 10, to = 20), 
     main = "", xlab = "", ylim = c(0,2.0), 
     lty = 3, lwd = 5, col = "Red")
par(new=T)
plot(density(keepersC[,2], from = 10, to = 20), 
     main = "", xlab = "", ylim = c(0,2.0), 
     lty = 2, lwd = 5, col = "Brown")
par(new=T)
points(VaR.True, 0, pch = "X")
legend("topright", title = "Sample Sizes", 
       legend =c("n=600", "n=400", "n=200"), 
       col = c("Green", "Red", "Brown"),
       bg="lightblue", lwd=5, lty=1:3)

#Plots of the CTE for different sample sizes. 
plot(density(keepersA[,3], from = 18, to = 28), 
     main="CTE Estimates", xlab = "CTE", ylim = c(0, 1.0), 
     lwd = 5, col = "Green")
par(new=T)
plot(density(keepersB[,3], from = 18, to = 28), 
     main = "", xlab = "", ylim = c(0,1.0), 
     lty = 3, lwd = 5, col = "Red")
par(new=T)
plot(density(keepersC[,3], from = 18, to = 28), 
     main = "", xlab = "", ylim = c(0,1.0), 
     lty = 2, lwd = 5, col = "Brown")
par(new=T)
points(CTE.True, 0, pch = "X")
legend("topright", title = "Sample Sizes", 
       legend =c("n=600", "n=400", "n=200"), 
       col = c("Green", "Red", "Brown"),
       bg="lightblue", lwd=5, lty=1:3)