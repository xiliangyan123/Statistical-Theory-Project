##-----------------------------------------------------------------------------------
## TITLE:     Model correctlt specified, Bayesian versus credibility estimation,
##            and Gibbs posterio
##
## AUTHORS:    L. Hong and R. Martin
##-----------------------------------------------------------------------------------

set.seed(55)
alpha <- 3 #shape
lambda <- 0.2 #rate parameter
n <- c(600, 400, 200)
N <- 500
h <- 0.2 #what is h correnspond for?

#found this on page 16, middle paragraph 
#They call 'h' the hazard rate. Note that this is different from the rate parameter. 
#if h = 0.20, then posterior will have gamma distribution with :
  #shape parameter (3n+1) and 
  #rate parameters 1/5 + sum(x[1:n]^-1)

#X <- rlnorm(n[1], meanlog=1, sdlog=1) #generate random number log normal distribution
X <- rgamma(n[1], shape=3, rate=lambda)
shape <- n*alpha+1
rate1 <- lambda+sum(X) #RATE for n = 600
rate2 <- lambda+sum(X[1:400]) #rate for n=400
rate3 <- lambda+sum(X[1:200]) #rate for n=200
mean(X) #credibility mean

#alpha is shape and lamda is 1/rate 

#---------------------------------------------------
# plots of the posterior for different sample sizes
#---------------------------------------------------

dpost1 <- function(x){dgamma(x, shape=shape[1], rate=rate1)}
curve(dpost1, from=0.15, to=0.25, xlab="parameter", ylab="posterior density", lty=1)
dpost2 <- function(x){dgamma(x, shape=shape[2], rate=rate2)}
curve(dpost2, from=0.15, to=0.25, xlab="", ylab="", lty=2, add=TRUE)
dpost3 <- function(x){dgamma(x, shape=shape[3], rate=rate3)}
curve(dpost3, from=0.15, to=0.25, xlab="", ylab="", lty=3, add=TRUE)

dpost4 <- function(x){dlnorm(x, meanlog = 1, sdlog = 1)}
curve(dpost4, from=0.15, to=0.25, xlab="parameter", ylab="posterior density", lty=1)
dpost5 <- function(x){dlnorm(x, meanlog = 2, sdlog = 2)}
curve(dpost5, from=0.15, to=0.25, xlab="parameter", ylab="posterior density", lty=2, add=TRUE)
dpost6 <- function(x){dlnorm(x, meanlog = 3, sdlog = 3)}
curve(dpost6, from=0.15, to=0.25, xlab="parameter", ylab="posterior density", lty=3, add=TRUE)
#points(0.669 , 0, pch = "X") #need to calculate K(P*, P) by calculus

#------------------------------------
# Posteriors of the predictive mean
#------------------------------------


#produce random # from gamma w/ different alpha and betha
Y1 <- rgamma(N, shape=shape[1], rate=rate1) 
Y2 <- rgamma(N, shape=shape[2], rate=rate2)
Y3 <- rgamma(N, shape=shape[3], rate=rate3)


mu1 <- alpha/Y1
mu2 <- alpha/Y2
mu3 <- alpha/Y3

#get function to produce density of standard normal
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

x <- seq(13.5, 17.5, 0.01) #there is 201 # generated from the sequence from 3.5 - 5.5 by 0.01
y1 <- sapply(x, f1.hat)
y2 <- sapply(x, f2.hat)
y3 <- sapply(x, f3.hat)

plot(x, y1, type="l", xlim=c(13.5, 17), ylim=c(0, 2), lty=1, xlab="mean", ylab="posterior density")
par(new=T)
plot(x, y2, type="l", xlim=c(13.5, 17), ylim=c(0, 2), lty=2, xlab="", ylab="")
par(new=T)
plot(x, y3, type="l", xlim=c(13.5, 17), ylim=c(0, 2), lty=3, xlab="", ylab="")
par(new=T)
points(15.21233 , 0, pch = "X") #credibility mean


#------------------
# posterior of VaR
#------------------

#VaR.star <- qlnorm(0.95, meanlog=1, sdlog=1)
VaR.star <- qgamma(0.95, shape=3, rate=lambda)
#d <- function(x){exp(-(log(x)-1)^2/2)/(sqrt(2*pi))}
d <- function(x){(lambda^3/gamma(3))*x^(3)*exp(-lambda*x)}
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


x <- seq(29.5, 34.5, 0.01)
y1 <- sapply(x, g1.hat)
y2 <- sapply(x, g2.hat)
y3 <- sapply(x, g3.hat)

plot(x, y1, type="l", xlim=c(29.5, 34.5),  ylim=c(0, 1.5), lty=1, xlab="VaR", ylab="posterior density")
par(new=T)
plot(x, y2, type="l", xlim=c(29.5, 34.5), ylim=c(0, 1.5), lty=2, xlab="", ylab="")
par(new=T)
plot(x, y3, type="l", xlim=c(29.5, 34.5), ylim=c(0, 1.5), lty=3, xlab="", ylab="")
par(new=T)
points(31.47897, 0, pch = "X") #this point is true VaR

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


x <- seq(35.5, 41.5, 0.01)
y1 <- sapply(x, h1.hat)
y2 <- sapply(x, h2.hat)
y3 <- sapply(x, h3.hat)

plot(x, y1, type="l", xlim=c(35.5, 41.5),  ylim=c(0, 1.5), lty=1, xlab="CTE", ylab="posterior density")
par(new=T)
plot(x, y2, type="l", xlim=c(35.5, 41.5), ylim=c(0, 1.5), lty=2, xlab="", ylab="")
par(new=T)
plot(x, y3, type="l", xlim=c(35.5, 41.5), ylim=c(0, 1.5), lty=3, xlab="", ylab="")
par(new=T)
points(38.0088, 0, pch = "X") #This point is true CTE



#simulating Gibbs
#Work in progress
set.seed(54)
X6 <- rlnorm(n[1], meanlog=1, sdlog=1)
X4 <- rlnorm(n[2], meanlog=1, sdlog=1)
X2 <- rlnorm(n[3], meanlog=1, sdlog=1)
mu6 <- mean(X6)
mu4 <- mean(X4)
mu2 <- mean(X2)
w6 <- 1/(2*var(X6))
w4 <- 1/(2*var(X4))
w2<- 1/(2*var(X2))
w <- cbind(w6,w4,w2)
sigma <- 1/(2*n*w+1)


#---------------------------------------------------
# plots of the posterior for different sample sizes
#---------------------------------------------------
?dnorm
dpost1 <- function(x){dnorm(x, mean=mu6, sd=sigma[1])}
curve(dpost1, from=3, to=6, xlab="parameter", ylab="posterior density", lty=1)
dpost2 <- function(x){dnorm(x, mean=mu4, sd=sigma[2])}
curve(dpost2, from=3, to=6, xlab="", ylab="", lty=2, add=TRUE)
dpost3 <- function(x){dnorm(x, mean=mu2, sd=sigma[3])}
curve(dpost3, from=3, to=6, xlab="", ylab="", lty=3, add=TRUE)
points(4.48 , 0, pch = "X")

