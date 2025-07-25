\name{ORBayesianTtest}
\alias{ORBayesianTtest}
\title{
A Objective and Robust Bayesian \eqn{t}-test
}
\description{
Compute a Bayes Factor based on the intrinsic, Jeffrey's and Berger's robust prior. For a population mean
}
\usage{
ORBayesianTtest(x, y = NULL, mu0 = 0, type = "Robust", paired=FALSE, pH0 = 0.5)
}
\arguments{
  \item{x}{a (non-empty) numeric vector of data values.}
  \item{y}{an optional (non-empty) numeric vector of data values.}
  \item{mu0}{a number indicating the true value of the mean.}
  \item{type}{Method of computing the Bayes Factor. Default is approach is Robust which based on the Berger's Robust prior, other choices are "Intrinsic" based on the intrinsic priors, "TESS", the BIC based on the  effective sample size, and "Jeffreys" based on Jeffreys prior.}
  \item{paired}{a logical indicating whether you want a paired samples.}
  \item{pH0}{Prior probability of the null hypothesis, i.e.,\eqn{P(H_0)}. Default is 0.5. }
}
\value{
A list containig the \eqn{t}, the Bayes Factor of \eqn{H_0} versus \eqn{H_1} and the posterior probability of \eqn{H_0}, i.e., \eqn{P(H_0|data)}.
}
\references{
Almodóvar-Rivera, I. A., & Pericchi-Guerra, L. R. (2024). An objective and robust bayes factor for the hypothesis test one sample and two population means. Entropy, 26(1), 88.
}
\author{
Israel A. Almodóvar-Rivera, Ph.D.
}
\examples{
\dontrun{
set.seed(787)
x <- rnorm(n=787)
## one sample
ORBayesianTtest(x=x,mu0 = 0,pH0 = 0.5)
## two-samples
set.seed(787)
x <- rnorm(n=787)
y <- rnorm(n=939)
ORBayesianTtest(x=x,y=y,mu0 = 0,pH0 = 0.5)
}
}