% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_theta_ig.r
\name{get_theta_ig}
\alias{get_theta_ig}
\title{Find Scale Parameter for Inverse Gamma Hyperprior}
\usage{
get_theta_ig(alpha = 0.01, method = "integrate", Z, c = 3,
  eps = .Machine$double.eps, Kinv, equals = FALSE, a = 1,
  type = "marginalt")
}
\arguments{
\item{alpha}{denotes the 1-\eqn{\alpha} level.}

\item{method}{with \code{integrate} as default.
Currently no further method implemented.}

\item{Z}{the design matrix.}

\item{c}{denotes the expected range of the function.}

\item{eps}{denotes the error tolerance of the result, default is \code{.Machine$double.eps}.}

\item{Kinv}{the generalised inverse of K.}

\item{equals}{saying whether \code{a}=\code{b}. The default is FALSE due to the fact that a is a shape parameter.}

\item{a}{is the shape parameter of the inverse gamma distribution, default is 1.}

\item{type}{is either numerical integration (\code{integrate}) of to obtain the marginal distribution of \eqn{z_p'\beta}
or the theoretical marginal t-distribution (\code{marginalt}). \code{marginalt} is the default value.}
}
\value{
an object of class \code{list} with values from \code{\link{uniroot}}.
}
\description{
This function implements a optimisation routine that computes the scale parameter \code{b}
of the inverse gamma prior for \eqn{\tau^2} when \eqn{a=b=\epsilon} with \eqn{\epsilon} small
for a given design matrix and prior precision matrix
such that approximately \eqn{P(|f(x_{k}|\le c,k=1,\ldots,p)\ge 1-\alpha}
When \code{a} unequal to \code{a} the shape parameter \code{a} has to be specified.
}
\details{
Currently, the implementation only works properly for the cases \code{a} unequal \code{b}.
}
\examples{
set.seed(123)
library(MASS)
# prior precision matrix (second order differences) 
# of a spline of degree l=3 and with m=20 inner knots
# yielding dim(K)=m+l-1=22
K <- t(diff(diag(22), differences=2))\%*\%diff(diag(22), differences=2)
# generalised inverse of K
Kinv <- ginv(K)
# covariate x
x <- runif(1)
Z <- matrix(DesignM(x)$Z_B,nrow=1)
theta <- get_theta_ig(alpha = 0.01, method = "integrate", Z = Z, 
                      c = 3, eps = .Machine$double.eps, Kinv = Kinv, 
					 equals = FALSE, a = 1, type="marginalt")$root

}
\references{
Nadja Klein and Thomas Kneib (2015). Scale-Dependent Priors for Variance Parameters in Structured Additive Distributional Regression. 
\emph{Working Paper}.

Stefan Lang and Andreas Brezger (2004). Bayesian P-Splines. 
\emph{Journal of Computational and Graphical Statistics}, \bold{13}, 183-212.
}
\author{
Nadja Klein
}
