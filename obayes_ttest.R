##*************************************************
##*
##* @file: obayes_ttest.R
##*
##*
##* Perform Jeffrey's, Intrinsic and 
##* Robust based on Berger's Priors Bayes Factors,
##* for the one sample hypothesis test for a 
##* population mean.
##*
##* Author:
##* 
##* Israel A. Almodovar-Rivera, PhD
##* University of Puerto Rico at Mayaguez
##* Department of Mathematical Sciences
##* @email: israel.almodovar@upr.edu
##*
##****************************************************

## Jeffrey's Bayes Factors

obayes.t.test <- function(x,y = NULL,mu0 = 0, method="robust", pH0 = 0.5, paired = FALSE, ...){

  method <- match.arg(method,choices = c("robust","intrinsic","Jeffreys","TESS"))
  if (!missing(mu0) && (length(mu0) != 1 || is.na(mu0))){
    stop("'mu0' must be a single number")
  }
  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    if (paired) 
      xok <- yok <- complete.cases(x, y)
    else {
      yok <- !is.na(y)
      xok <- !is.na(x)
    }
    y <- y[ok]
  } else {
    dname <- deparse(substitute(x))
    if (paired) 
      stop("'y' is missing for paired test")
    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]
  if (paired) {
    x <- x - y
    y <- NULL
  }
   
  n1 <- length(x)
  xbar <- mean(x)
  Sx <- var(x)
  
  p0 <- pH0
  p1 <- 1-p0
  
  if(is.null(y)) {
    if(n1 < 2) stop("not enough 'x' observations")
    nu <- n1-1
    stderr <- sqrt(Sx/n1)
    if(stderr < 10 *.Machine$double.eps * abs(xbar))
      stop("data are essentially constant")
    tstat <- (xbar-mu0)/stderr
    method <- if(paired) "Paired Bayesian t-test" else "One Sample Bayesian t-test"
    estimate <-
      setNames(xbar, if(paired)"mean of the differences" else "mean of x")
    
    if(method=="Jeffrey"){
      bf <- sqrt(pi*nu/2) * (1+tstat^2/nu)^(-(nu-1)/2)
    } else if(method=="intrinsic"){
      bf <- sqrt(2*n1) * (1+tstat^2/nu)^(-n1/2) * tstat^2/nu *1/(1-exp(-tstat^2/nu))
    } else if (method=="robust"){
      bf <- sqrt(2/(n1+1)) * ((n1-2)/(n1-1)) * tstat^2 * (1+tstat^2/nu)^(-n1/2)*(1-(1+2*tstat^2/(n1^2-1))^(-(n1-2)/2))^(-1)
    }
    
    pH0.data <- (1+p1/p0 * 1/bf)^(-1)
    
  } else {
    n2 <- length(y)
    if(n1 < 1 || (!var.equal && n1 < 2))
      stop("not enough 'x' observations")
    if(n2 < 1 || (!var.equal && n2 < 2))
      stop("not enough 'y' observations")
    if(var.equal && n1+n2 < 3) stop("not enough observations")
    ybar <- mean(y)
    Sy <- var(y)
    method <- "Two Sample Bayesian t-test"#paste("Two Sample", met ,"t-test")
    estimate <- c(xbar,ybar)
    names(estimate) <- c("mean of x","mean of y")

    nu <- n1+n2-2
    v <- 0
    if(n1 > 1) v <- v + (n2-1)*S1
    if(n1 > 1) v <- v + (n2-1)*S2
    v <- v/nu
    stderr <- sqrt(v*(1/n1+1/n2))
     
    if(stderr < 10 *.Machine$double.eps * max(abs(xbar), abs(ybar)))
      stop("data are essentially constant")
    
    tstat <- (xbar - ybar - mu0)/stderr
    
  } 
  
  names(tstat) <- "t"
  names(nu) <- "nu"
  names(mu0) <- if(paired || !is.null(y)) "difference in means" else "mean"
  ##**********
  ##* output
  ##*********
  
  rval <- list(statistic = tstat, parameter = nu, posterior = pH0.data,
               cestimate = estimate, null.value = mu0,
               method = method, data.name = dname)
   return(rval)
  
}

Bayesian.t.test <- function(dat1,dat2 = NULL, alternative = c("two.sided", "less", "greater"),betastar = 0,prior = c(0.5,0.5),paired = FALSE,...) 
{
  alternative <- match.arg(alternative)
  if (!missing(betastar) && (length(betastar) != 1 || is.na(betastar))) 
    stop("'betastar' must be a single number")
  if (!is.null(dat2)) {
    dname <- paste(deparse(substitute(dat1)), "and", deparse(substitute(dat2)))
    if (paired) 
      dat1ok <- dat2ok <- complete.cases(dat1, dat2)
    else {
      dat2ok <- !is.na(dat2)
      dat1ok <- !is.na(dat1)
    }
    dat2 <- dat2[dat2ok]
  }
  else {
    dname <- deparse(substitute(dat1))
    if (paired) 
      stop("'data 2' is missing for paired test")
    dat1ok <- !is.na(dat1)
    dat2ok <- NULL
  }
  dat1 <- dat1[dat1ok]
  if (paired) {
    dat1 <- dat1 - dat2
    dat2 <- NULL
  }
  n1 <- length(dat1)
  ybar1 <- mean(dat1)
  S1 <- var(dat1)
  if (is.null(dat2)) {
    if (n1 < 2) 
      stop("not enough 'data 1' observations")
    df <- n1 - 1
    stderr <- sqrt(S1/n1)
    if (stderr < 10 * .Machine$double.eps * abs(ybar1)) 
      stop("data are essentially constant")
    tvalue <- (ybar1 - betastar)/stderr
    method <- ifelse(paired, "Paired Robust Bayesian t-test", "Robust Bayesian One Sample t-test")
    
    bfs <- sqrt(n1)*(1+tvalue^2/(n1-1))^(-n1/2)
    bfr <- sqrt(n1)^(-1) * bfs * sqrt(2/(n1+1))*(n1-2)/(n1-1) * tvalue^2 * (1-(1+2*tvalue^2/(n1^2-1))^(-(n1-2)/2))^(-1)
    pH0 <- bfr/(1+prior[2]/prior[1] *bfr)
    names(ybar1) <- ifelse(paired, "mean of the differences", "mean of Data 1")
    names(tvalue) <- "t"
    names(betastar) <- if (paired || !is.null(dat2)) 
      "difference in means"
    else "mean"
    btt <- list(statistic = tvalue, posterior = pH0 , null.value = betastar, alternative = alternative, method = method, data.name = dname)
    class(btt) <- "bhtest"
    return(btt)
  }
  else {
    n2 <- length(dat2)
    if (n1 < 1 || (n1 < 2)) 
      stop("not enough 'data 1' observations")
    if (n2 < 1 || (n2 < 2)) 
      stop("not enough 'data 2' observations")
    if (n1 + n2 < 3) 
      stop("not enough observations")
    ybar2 <- mean(dat2)
    S2 <- var(dat2)
    method <- paste("Robust Bayesian", "Two Sample t-test")
    n <- n1+n2
    ndelta <- n1*n2 /n
    diff <- ybar1-ybar2
    Spool <- ((n1-1)*S1+(n2-1)*S2)/(n-2)
    d <- ndelta
    b <- (max(1/n1,abs(-1/n2)))^(-2)
    
    if (sqrt(Spool/ndelta) < 10 * .Machine$double.eps * max(abs(ybar1), abs(ybar2))) 
      stop("data are essentially constant")
    tvalue <- sqrt(ndelta)* abs(diff-betastar)/ sqrt(Spool)
  }
  if (alternative == "less") {
    const <- S1+S2+n1*ybar1^2+n2*ybar2 -((n1*ybar1+n2*ybar2)^2 - (n1*ybar1-n2*ybar2)^2)/n
    a <- (n1*ybar1 + n2*ybar2)/n
    b <- (n1 * ybar1 - n2*ybar2)/n
    integrad <- function(x) { exp(-log(2*sqrt(n)) +lgamma((n-1)/2)-(n-1)/2*log(pi)-(n-1)/2*log(const+b^2*n-2*a*((n2-n1)+b*n)*x + 4*ndelta*x^2))}
    B01 <- integrate(integrad,lower = 0,upper = Inf)$value /integrate(integrad,lower = -Inf,upper = 0)$value
    pH0 <- B01/(1+B01)        
  }
  else if (alternative == "greater") {
    const <- S1+S2+n1*ybar1^2+n2*ybar2 -((n1*ybar1+n2*ybar2)^2 - (n1*ybar1-n2*ybar2)^2)/n
    a <- (n1*ybar1 + n2*ybar2)/n
    b <- (n1 * ybar1 - n2*ybar2)/n
    integrad <- function(x) { exp(-log(2*sqrt(n)) +lgamma((n-1)/2)-(n-1)/2*log(pi)-(n-1)/2*log(const+b^2*n-2*a*((n2-n1)+b*n)*x + 4*ndelta*x^2))}
    B01 <- integrate(integrad,lower = -Inf,upper = 0)$value /integrate(integrad,lower = 0,upper = Inf)$value
    pH0 <- B01/(1+prior[2]/prior[1] * B01)        
  }
  else {
    bfs  <- (1+tvalue^2/(n-2))^(-(n-1)/2)
    bfr <- sqrt(8 * ndelta/(d+b)) * exp(lgamma((n-1)/2)-lgamma((n-3)/2) ) * tvalue^2/(n-2) * bfs * (1- (1 + (2*ndelta*tvalue^2/((n-2)*(d+b))) )^(-(n-3)/2) )^(-1)
    pH0 <- bfr/(1+prior[2]/prior[1] *bfr)
  }
  names(tvalue) <- "t"
  names(betastar) <- if (paired || !is.null(dat2)) 
    "difference in means"
  else "mean"
  btt <- list(statistic = tvalue, posterior = pH0 ,null.value = betastar, alternative = alternative, method = method, data.name = dname)
  class(btt) <- "bhtest"
  return(btt)
}
