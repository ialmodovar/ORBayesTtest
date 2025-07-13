\name{ORBayesianTtest.sm}
\alias{ORBayesianTtest.sm}
\title{
A Objective and Robust Bayesian \eqn{t}-test
}
\description{
Compute a Bayes Factor based on the intrinsic, Jeffrey's and Berger's robust prior using summaries statistics.
}
\usage{
ORBayesianTtest.sm(xbar1, s1, n1, xbar2 = NULL, s2 = NULL, n2 = NULL, mu0 = 0, type="Robust", pH0 = 0.5)
}
\arguments{
  \item{xbar1}{sample mean of group 1}
  \item{s1}{sample standard deviation of group 1}
  \item{n1}{sample size of group 1}
  \item{xbar2}{sample mean of group 2}
  \item{s2}{sample standard deviation of group 2}
  \item{n2}{sample size of group 2}
  \item{mu0}{a number indicating the true value of the mean.}
  \item{type}{Method of computing the Bayes Factor. Default is approach is Robust which based on the Berger's Robust prior, other choices are "Intrinsic" based on the intrinsic priors, "TESS", the BIC based on the  effective sample size, and "Jeffreys" based on Jeffreys prior.}
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
## one sample
ORBayesianTtest.sm(xbar1=1,s1=0.5,n1=50)
## two sample
ORBayesianTtest.sm(xbar1=1,s1=0.5,n1=50,xbar2=9,s2=10,n2=5)
}
}