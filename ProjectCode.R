##-----------------------------------------------------------------------------------
## TITLE:     Model misspecification, Bayesian versus credibility estimation,
##            and Gibbs posterior
##
## AUTHORS:    L. Hong and R. Martin
##-----------------------------------------------------------------------------------

#Information to generate datasets. 
#Generating Lognormal Data while fitting a Gamma Model

set.seed(53) 
alpha <- 3 
lambda <- 0.2 
n <- c(600, 400, 200) 
N <- 500 
h <- 0.2 
X <- rlnorm(n[1], meanlog=1, sdlog=1) 

shape <- n*alpha+1 
rate1 <- lambda+sum(X) 
rate2 <- lambda+sum(X[1:400])
rate3 <- lambda+sum(X[1:200])
mean(X) #posterior mean of the log normal distribution. 

#---------------------------------------------------
# plots of the posterior for different sample sizes
#---------------------------------------------------

## Plot of the kull-back liebler minimizer = 0.669
dpost1 <- function(x){dgamma(x, shape=shape[1], rate=rate1)}
curve(dpost1, from=0.5, to=0.8, xlab="parameter", ylab="posterior density", lty=1)
dpost2 <- function(x){dgamma(x, shape=shape[2], rate=rate2)}
curve(dpost2, from=0.5, to=0.8, xlab="", ylab="", lty=2, add=TRUE)
dpost3 <- function(x){dgamma(x, shape=shape[3], rate=rate3)}
curve(dpost3, from=0.5, to=0.8, xlab="", ylab="", lty=3, add=TRUE)
points(0.669 , 0, pch = "X")

#------------------------------------
# Posteriors of the predictive mean
#------------------------------------

Y1 <- rgamma(N, shape=shape[1], rate=rate1)
Y2 <- rgamma(N, shape=shape[2], rate=rate2)
Y3 <- rgamma(N, shape=shape[3], rate=rate3)

mu1 <- alpha/Y1
mu2 <- alpha/Y2
mu3 <- alpha/Y3

f1.hat <- function(x){
  y <- mean(dnorm((x-mu1)/h))/h
  return(y)
}

f2.hat <- function(x){
  y <- mean(dnorm((x-mu2)/h))/h
  return(y)
}

f3.hat <- function(x){
  y <- mean(dnorm((x-mu3)/h))/h
  return(y)
}

x <- seq(3.5, 5.5, 0.01)
y1 <- sapply(x, f1.hat)
y2 <- sapply(x, f2.hat)
y3 <- sapply(x, f3.hat)

plot(x, y1, type="l", xlim=c(3.5, 5.5), ylim=c(0, 2), lty=1, xlab="mean", ylab="posterior density", 
     main = "Plot of the posterior mean", col = "Green", lwd = 5)
par(new=T)
plot(x, y2, type="l", xlim=c(3.5, 5.5), ylim=c(0, 2), lty=2, xlab="", ylab="", col = "Red", lwd = 5)
par(new=T)
plot(x, y3, type="l", xlim=c(3.5, 5.5), ylim=c(0, 2), lty=3, xlab="", ylab="", col = "Brown", lwd = 5)
par(new=T)
points(4.482 , 0, pch = "X")
legend("topleft", title = "Sample Sizes", 
       legend =c("n=600", "n=400", "n=200"), 
       col = c("Green", "Red", "Brown"),
       bg="lightblue", lwd=5, lty=1:3)
#------------------
# posterior of VaR
#------------------

VaR.star <- qlnorm(0.95, meanlog=1, sdlog=1)
d <- function(x){exp(-(log(x)-1)^2/2)/(sqrt(2*pi))}
CTE.star <- integrate(d, lower=VaR.star, upper=Inf)$value/0.05

g <- function(x){
  y <- qgamma(0.95, shape=3, rate=x)
  return(y)
}
VaR1 <- sapply(Y1, g)
VaR2 <- sapply(Y2, g)
VaR3 <- sapply(Y3, g)

g1.hat <- function(x){ 
  y <- mean(dnorm((x-VaR1)/h))/h
  return(y)
}

g2.hat <- function(x){ 
  y <- mean(dnorm((x-VaR2)/h))/h
  return(y)
}

g3.hat <- function(x){ 
  y <- mean(dnorm((x-VaR3)/h))/h
  return(y)
}

x <- seq(7, 12, 0.01)
y1 <- sapply(x, g1.hat)
y2 <- sapply(x, g2.hat)
y3 <- sapply(x, g3.hat)

plot(x, y1, type="l", xlim=c(7, 12), ylim=c(0, 1.5), lty=1, xlab="VaR", ylab="posterior density", 
     main = "Plot of the VaR estimate", col = "Green", lwd = 5)
par(new=T)
plot(x, y2, type="l", xlim=c(7, 12), ylim=c(0, 1.5), lty=2, xlab="", ylab="", col = "Red", lwd = 5)
par(new=T)
plot(x, y3, type="l", xlim=c(7, 12), ylim=c(0, 1.5), lty=3, xlab="", ylab="", col = "Brown", lwd = 5)
par(new=T)
points(9.411, 0, pch = "X")
legend("topleft", title = "Sample Sizes", 
       legend =c("n=600", "n=400", "n=200"), 
       col = c("Green", "Red", "Brown"),
       bg="lightblue", lwd=5, lty=1:3)

#-----------------------------------
# posterior of CTE
#-----------------------------------

b <- function(theta){
  integrand<-function(x){x*dgamma(x, shape=3, rate=theta)}
  y<-integrate(integrand, lower=qgamma(0.95, 3, theta), upper=Inf)$value/0.05  
  return(y)
}

CTE1 <- sapply(Y1, b)
CTE2 <- sapply(Y2, b)
CTE3 <- sapply(Y3, b)

h1.hat<-function(x){
  y<-mean(dnorm((x-CTE1)/h))/h
  return(y)
}

h2.hat<-function(x){
  y<-mean(dnorm((x-CTE2)/h))/h
  return(y)
}

h3.hat<-function(x){
  y<-mean(dnorm((x-CTE3)/h))/h
  return(y)
}

x <- seq(9, 15, 0.01)
y1 <- sapply(x, h1.hat)
y2 <- sapply(x, h2.hat)
y3 <- sapply(x, h3.hat)

plot(x, y1, type="l", xlim=c(9, 15),  ylim=c(0, 1.5), lty=1, xlab="CTE", ylab="posterior density", 
     main = "Plot of the CTE estimate", col = "Green", lwd = 5)
par(new=T)
plot(x, y2, type="l", xlim=c(9, 15), ylim=c(0, 1.5), lty=2, xlab="", ylab="", col = "Red", lwd = 5)
par(new=T)
plot(x, y3, type="l", xlim=c(9, 15), ylim=c(0, 1.5), lty=3, xlab="", ylab="", col = "Brown", lwd = 5)
par(new=T)
points(11.363, 0, pch = "X")
legend("topleft", title = "Sample Sizes", 
       legend =c("n=600", "n=400", "n=200"), 
       col = c("Green", "Red", "Brown"),
       bg="lightblue", lwd=5, lty=1:3)
