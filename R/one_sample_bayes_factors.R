##*********************************************
##*
##* @file: one_sample_bayes_factors.R
##*
##* One sample inference bayes factors
##*
##* Author:
##* Israel Almodovar-Rivera PhD
##* Department of Mathematical Sciences
##* University of Puerto Rico at Mayaguez
##* israel.almodovar@upr.edu
##* Copyright July 2025
##*********************************************


##**********************
##* One sample mean inference
##***********************

coth <- function(x){(exp(2*x)+1)/(exp(2*x)-1)}

one.sample.t.test.bf <- function(x,mu0,p0 = 0.5){
  ## remove NA
  x <- x[!is.na(x)]
  xbar <- mean(x)
  n <- length(x)
  sx <- sd(x)
  tstat <- (xbar-mu0)/(sx/sqrt(n))
  nu <- n-1
  ## intrinsic bayes factor
  bfi <- sqrt(2*n) * (1+tstat^2/nu)^(-n/2) * tstat^2/nu *1/(1-exp(-tstat^2/nu))
  ## robust bayes factor
  bfr <-sqrt(2/(n+1)) * ((n-2)/(n-1)) * tstat^2 * (1+tstat^2/nu)^(-n/2)*(1-(1+2*tstat^2/(n^2-1))^(-(n-2)/2))^(-1)
  ## jeffreys
  bfj <- sqrt(pi*nu/2) * (1+tstat^2/nu)^(-(nu-1)/2)
  
  
  ## posterior distributions of H_0
  
  
  p1 <- 1-p0
  
  pH0i.data <- (1+p1/p0 * 1/bfi)^(-1)
  pH0r.data <- (1+p1/p0 * 1/bfr)^(-1)
  pH0j.data <- (1+p1/p0 * 1/bfj)^(-1)
  
  bfs <- data.frame(statistic=tstat,BF01 = c(bfi,bfr,bfj),log2BF01 = 2*log(c(bfi,bfr,bfj)),posterior = c(pH0i.data,pH0r.data,pH0i.data))
  
  names(bfs) <- c("t","\\( B_{01} \\)","\\( 2 \\log (B_{01}) \\)",paste("\\(P (H_0: \\mu =", mu0,"|data) \\)"))
  row.names(bfs) <- c("Intrinsic","Robust","Jeffreys")
  
  bfs 
}

one.sample.t.test.bf.sm <- function(xbar,sx,n,mu0,p0 = 0.5){
  tstat <- (xbar-mu0)/(sx/sqrt(n))
  nu <- n-1
  ## intrinsic bayes factor
  bfi <- sqrt(2*n) * (1+tstat^2/nu)^(-n/2) * tstat^2/nu *1/(1-exp(-tstat^2/nu))
  ## robust bayes factor
  bfr <-sqrt(2/(n+1)) * ((n-2)/(n-1)) * tstat^2 * (1+tstat^2/nu)^(-n/2)*(1-(1+2*tstat^2/(n^2-1))^(-(n-2)/2))^(-1)
  ## jeffreys
  bfj <- sqrt(pi*nu/2) * (1+tstat^2/nu)^(-(nu-1)/2)
  
  
  ## posterior distributions of H_0
  
  
  p1 <- 1-p0
  
  pH0i.data <- (1+p1/p0 * 1/bfi)^(-1)
  pH0r.data <- (1+p1/p0 * 1/bfr)^(-1)
  pH0j.data <- (1+p1/p0 * 1/bfj)^(-1)
  
  bfs <- data.frame(statistic=tstat,BF01 = c(bfi,bfr,bfj),log2BF01 = 2*log(c(bfi,bfr,bfj)),posterior = c(pH0i.data,pH0r.data,pH0i.data))
  
  names(bfs) <- c("t","\\( B_{01} \\)","\\( 2 \\log (B_{01}) \\)",paste("\\(P (H_0: \\mu =", mu0,"|data) \\)"))
  row.names(bfs) <- c("Intrinsic","Robust","Jeffreys")
  
  bfs 
}

