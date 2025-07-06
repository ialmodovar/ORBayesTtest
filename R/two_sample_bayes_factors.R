##*********************************************
##*
##* @file: two_sample_bayes_factors.R
##*
##* Two sample inference bayes factors
##*
##* Author:
##* Israel Almodovar-Rivera PhD
##* Department of Mathematical Sciences
##* University of Puerto Rico at Mayaguez
##* israel.almodovar@upr.edu
##* Copyright July 2025
##*********************************************



##****************************
##* Two sample means comparison
##****************************


two.sample.t.test.bf <- function(x,y,paired=FALSE,mu0,p0 = 0.5,lambda =0, sigma2d = 1/3){
  
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
    
    ##
    
    postv <- 1+nd * sigma2d
    nc <- sqrt(nd/postv) * lambda
    bfg <- (dt(tstat,nu)/(dt(tstat/sqrt(postv),nu,nc)/sqrt(postv)))
    
    
    ## posterior distributions of H_0
    p1 <- 1-p0
    
    pH0i.data <- (1+p1/p0 * 1/bfi)^(-1)
    pH0r.data <- (1+p1/p0 * 1/bfr)^(-1)
    pH0j.data <- (1+p1/p0 * 1/bfj)^(-1)
    pH0t.data <- (1+p1/p0 * 1/bft)^(-1)
    pH0s.data <- (1+p1/p0 * 1/bfs)^(-1)
    pH0g.data <- (1+p1/p0 * 1/bfg)^(-1)
    
    
    bf.all <- data.frame(statistic=tstat,
                         BF01 = c(bfi,bfr,bfj,bft,bfs,bfg),
                         logBF01 = 2*log(c(bfi,bfr,bfj,bft,bfs,bfg)),
                         posterior = c(pH0i.data,pH0r.data,pH0j.data,pH0t.data,pH0s.data,pH0g.data))
    
    names(bf.all) <- c("t","\\( B_{01} \\)","\\(  2\\log (B_{01}) \\)",paste("\\(P (H_0: \\mu_1 -\\mu_2 = ", mu0,"|data) \\)"))
    row.names(bf.all) <- c("Intrinsic","Robust","Jeffreys","TESS","Schwarz","Gonen")
    
    bf.all 
  }
}

two.sample.t.test.bf.sm <- function(xbar1,s21,n1,xbar2,s22,n2,mu0,p0 = 0.5,lambda =0, sigma2d = 1/3)
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
  
  ##
  
  postv <- 1+nd * sigma2d
  nc <- sqrt(nd/postv) * lambda
  bfg <- (dt(tstat,nu)/(dt(tstat/sqrt(postv),nu,nc)/sqrt(postv)))
  
  
  ## posterior distributions of H_0
  p1 <- 1-p0
  
  pH0i.data <- (1+p1/p0 * 1/bfi)^(-1)
  pH0r.data <- (1+p1/p0 * 1/bfr)^(-1)
  pH0j.data <- (1+p1/p0 * 1/bfj)^(-1)
  pH0t.data <- (1+p1/p0 * 1/bft)^(-1)
  pH0s.data <- (1+p1/p0 * 1/bfs)^(-1)
  pH0g.data <- (1+p1/p0 * 1/bfg)^(-1)
  
  
  bf.all <- data.frame(statistic=tstat,
                       BF01 = c(bfi,bfr,bft,bfj,bfs,bfg),
                       logBF01 = 2*log(c(bfi,bfr,bft,bfj,bfs,bfg)),
                       posterior = c(pH0i.data,pH0r.data,pH0t.data,pH0j.data,pH0s.data,pH0g.data))
  
  names(bf.all) <- c("t","\\( B_{01} \\)","\\(  2\\log (B_{01}) \\)",paste("\\(P (H_0: \\mu_1 -\\mu_2 = ", mu0,"|data) \\)"))
  row.names(bf.all) <- c("Intrinsic","Robust","TESS","Jeffreys","Schwarz","Gonen")
  
  bf.all 
}
