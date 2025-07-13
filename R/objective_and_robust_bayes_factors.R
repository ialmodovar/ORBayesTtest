##*********************************************
##*
##* @file: objective_and_robust_bayes_factors.R
##*
##* Objective and Robust Bayes factors
##*
##* Author:
##* Israel Almodovar-Rivera PhD
##* Department of Mathematical Sciences
##* University of Puerto Rico at Mayaguez
##* israel.almodovar@upr.edu
##* Copyright July 2025
##*********************************************

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
  p1 <- 1-p0
  
  pH0i.data <- (1+p1/p0 * 1/bfi)^(-1)
  pH0r.data <- (1+p1/p0 * 1/bfr)^(-1)
  pH0j.data <- (1+p1/p0 * 1/bfj)^(-1)
  
  bfs <- data.frame(statistic=tstat,BF01 = c(bfi,bfr,bfj),log2BF01 = 2*log(c(bfi,bfr,bfj)),posterior = c(pH0i.data,pH0r.data,pH0i.data))
  
  names(bfs) <- c("t","\\( B_{01} \\)","\\( 2 \\log (B_{01}) \\)",paste("\\(P (H_0: \\mu =", mu0,"|data) \\)"))
  row.names(bfs) <- c("Intrinsic","Robust","Jeffreys")
  
  bfs 
}

two.sample.t.test.bf <- function(x,y,paired=FALSE,mu0,p0 = 0.5){
  
  if(paired){
    one.sample.t.test.bf(x-y,mu0 = mu0,p0 = p0)
  } else{
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    s21 <- var(x)
    s22 <- var(y)
    xbar1 <- mean(x)
    xbar2 <- mean(y)
    n1 <- length(x)
    n2 <- length(y)
    n <- n1+n2
    nu <- n-2
    nd <- (1/n1+1/n2)
    
    Sp <- ((n1-1)*s21+(n2-1)*s22)/(nu)
    tstat <- (xbar1-xbar2)/sqrt(Sp*nd)
    
    bfj <-  sqrt(pi*1/nd/2) * (1+tstat^2/(nu))^(-(n-1)/2)
    bfi <- tstat^2* n/(nu)*sqrt(nd)* (1+tstat^2/nu)^(-(n-1)/2)*(coth(nd*tstat^2/(nu))+1)
    ##----
    d <- nd/4
    b <- max(abs(1/n1),abs(-1/n2))^(-2)/4
    bfr <- sqrt(8/(b+d)*nd) * (n-3) * (tstat^2/(4*nu)) * (1 + tstat^2/nu)^(-(n-1)/2) *(1 - (1 + tstat^2*nd/(2*nu*(b+d)))^(-(n-3)/2))^(-1)
    
    ne0 <- (max(abs(1/n1),abs(-1/n2)))^(-2)* 1/nd
    bft <- sqrt(ne0) * (1+tstat^2/nu)^(-n/2)
    
    bfs <- sqrt(n) * (1+tstat^2/nu)^(-n/2)
    
    
    ## posterior distributions of H_0
    p1 <- 1-p0
    
    pH0i.data <- (1+p1/p0 * 1/bfi)^(-1)
    pH0r.data <- (1+p1/p0 * 1/bfr)^(-1)
    pH0j.data <- (1+p1/p0 * 1/bfj)^(-1)
    pH0t.data <- (1+p1/p0 * 1/bft)^(-1)
    pH0s.data <- (1+p1/p0 * 1/bfs)^(-1)
    
    bf.all <- data.frame(statistic=tstat,
                         BF01 = c(bfi,bfr,bfj,bft,bfs),
                         logBF01 = 2*log(c(bfi,bfr,bfj,bft,bfs)),
                         posterior = c(pH0i.data,pH0r.data,pH0j.data,pH0t.data,pH0s.data))
    
    names(bf.all) <- c("t","\\( B_{01} \\)","\\(  2\\log (B_{01}) \\)",paste("\\(P (H_0: \\mu_1 -\\mu_2 = ", mu0,"|data) \\)"))
    row.names(bf.all) <- c("Intrinsic","Robust","Jeffreys","TESS","Schwarz")
    
    bf.all 
  }
}

ORBayesianTtest <- function(x,y = NULL, mu0 = 0, paired=FALSE,pH0 = 0.5){
  if(is.null(y)){
    res <- one.sample.t.test.bf(x = x, mu0 = mu0, p0 = pH0)
  } else if(!is.null(y)){
    res <- two.sample.t.test.bf(x = x, y = y,paired=paired,mu0 = mu0,p0 = pH0)
  }
  res
}

##*****************
##* Summary
##******************

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


two.sample.t.test.bf.sm <- function(xbar1,s21,n1,xbar2,s22,n2,mu0,p0 = 0.5)
{
  
  n <- n1+n2
  nu <- n-2
  nd <- (1/n1+1/n2)
  
  Sp <- ((n1-1)*s21+(n2-1)*s22)/(nu)
  tstat <- (xbar1-xbar2-mu0)/sqrt(Sp*nd)
  
  bfj <-  sqrt(pi*1/nd/2) * (1+tstat^2/(nu))^(-(n-1)/2)
  bfi <- tstat^2* n/(nu)*sqrt(nd)* (1+tstat^2/nu)^(-(n-1)/2)*(coth(nd*tstat^2/(nu))+1)
  ##----
  d <- nd/4
  b <- max(abs(1/n1),abs(-1/n2))^(-2)/4
  bfr <- sqrt(8/(b+d)*nd) * (n-3) * (tstat^2/(4*nu)) * (1 + tstat^2/nu)^(-(n-1)/2) *(1 - (1 + tstat^2*nd/(2*nu*(b+d)))^(-(n-3)/2))^(-1)
  
  ne0 <- (max(abs(1/n1),abs(-1/n2)))^(-2)* 1/nd
  bft <- sqrt(ne0) * (1+tstat^2/nu)^(-n/2)
  
  bfs <- sqrt(n) * (1+tstat^2/nu)^(-n/2)
  
  p1 <- 1-p0
  
  pH0i.data <- (1+p1/p0 * 1/bfi)^(-1)
  pH0r.data <- (1+p1/p0 * 1/bfr)^(-1)
  pH0j.data <- (1+p1/p0 * 1/bfj)^(-1)
  pH0t.data <- (1+p1/p0 * 1/bft)^(-1)
  pH0s.data <- (1+p1/p0 * 1/bfs)^(-1)

  bf.all <- data.frame(statistic=tstat,
                       BF01 = c(bfi,bfr,bft,bfj,bfs),
                       logBF01 = 2*log(c(bfi,bfr,bft,bfj,bfs)),
                       posterior = c(pH0i.data,pH0r.data,pH0t.data,pH0j.data,pH0s.data))
  
  names(bf.all) <- c("t","\\( B_{01} \\)","\\(  2\\log (B_{01}) \\)",paste("\\(P (H_0: \\mu_1 -\\mu_2 = ", mu0,"|data) \\)"))
  row.names(bf.all) <- c("Intrinsic","Robust","TESS","Jeffreys","Schwarz")
  
  bf.all 
}

ORBayesianTtest.sm <- function(xbar1,s1,n1,xbar2 = NULL,s2 = NULL,n2 = NULL, mu0 = 0, pH0 = 0.5){
  
  if(is.null(xbar2) & is.null(s2) & is.null(n2) ){
    res <- one.sample.t.test.bf.sm(xbar = xbar1, sx = s1, n = n1, mu0 = mu0, p0 = pH0)
  } else if(!is.null(xbar2) & !is.null(s2) & !is.null(n2)){
    res <- two.sample.t.test.bf.sm(xbar1 = xbar1, s21  = s1^2, n1 = n1,xbar2 = xbar2,s22 = s2^2,n2 = n2,mu0 = mu0,p0 = pH0)
  }
  res
}