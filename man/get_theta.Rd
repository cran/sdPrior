% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_theta.r
\name{get_theta}
\alias{get_theta}
\title{Find Scale Parameter for (Scale Dependent) Hyperprior}
\usage{
get_theta(alpha = 0.01, method = "integrate", Z, c = 3,
  eps = .Machine$double.eps, Kinv)
}
\arguments{
\item{alpha}{denotes the 1-\eqn{\alpha} level.}

\item{method}{either \code{integrate} or \code{trapezoid} with \code{integrate} as default.
\code{trapezoid} is a self-implemented version of the trapezoid rule.}

\item{Z}{the design matrix.}

\item{c}{denotes the expected range of the function.}

\item{eps}{denotes the error tolerance of the result, default is \code{.Machine$double.eps}.}

\item{Kinv}{the generalised inverse of K.}
}
\value{
an object of class \code{list} with values from \code{\link{uniroot}}.
}
\description{
This function implements a optimisation routine that computes the scale parameter \eqn{\theta}
of the scale dependent hyperprior for a given design matrix and prior precision matrix
such that approximately \eqn{P(|f(x_{k}|\le c,k=1,\ldots,p)\ge 1-\alpha}
}
\examples{
\dontrun{

set.seed(91179)
library(BayesX)
library(MASS)
# prior precision matrix to zambia data set
K <- read.gra(system.file("examples/zambia.gra", package="sdPrior"))
# generalised inverse of K
Kinv <- ginv(K)

# read data
dat <- read.table(system.file("examples/zambia_height92.raw", package="sdPrior"), header = TRUE)

# design matrix for spatial component
Z <- t(sapply(dat$district, FUN=function(x){1*(x==rownames(K))}))

# get scale parameter
theta <- get_theta(alpha = 0.01, method = "integrate", Z = Z, 
                            c = 3, eps = .Machine$double.eps, Kinv = Kinv)$root
} 

}
\references{
Nadja Klein and Thomas Kneib (2015). Scale-Dependent Priors for Variance Parameters in Structured Additive Distributional Regression. 
\emph{Working Paper}.
}
\author{
Nadja Klein
}
