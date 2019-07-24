\name{getSigma}
\alias{getSigma}
\alias{getRho}

\title{
Covariance matrix and correlation estimation for repeated measurements
}

\description{
\code{getRho} is an implementation of the proposed moment estimator (see reference) for estimating the unknown correlation in repeated measurements. Based on the output from \code{getRho}, the \code{getSigma} function returns an estimated covariance matrix.
}

\usage{
getSigma(YY, xx, id, weights=NULL, correction=TRUE)
getRho(YY, xx, id, weights, correction=TRUE, return.all=FALSE)
}


\arguments{
  \item{YY}{Vector or matrix of response values. If gene expression matrix, genes in columns and subjects in rows.}
  \item{xx}{Vector of covariate values.}
  \item{id}{Vector of subject IDs (for repeated measurements).}
  \item{weights}{Vector of \emph{a priori} weights. If \code{NULL}, equal weights of 1 will be assigned.}
  \item{correction}{Logical value. If \code{TRUE}, the correlation estimate is corrected by the Olkin and Pratt's method for downward bias in small samples.}
  \item{return.all}{Logical value. If \code{TRUE}, individual correlation estimates for each column of the \code{reponse} matrix are returned. If \code{FALSE}, average correlation across columns of the \code{reponse} matrix is returned. In \code{getSigma()}, the average correlation is used.}
}

\value{
The \code{getSigma} function returns a variance-covariance matrix based on the correlation and weights.

The \code{getRho} function returns the individual or average correlation estimate(s).
}

\references{
  \cite{Zhang et al. (2019) Highly efficient hypothesis testing methods for regression-type tests with correlated observations and heterogeneous variance structure. BMC Bioinformatics, 20:185.}
}

\author{
Yun Zhang, Xing Qiu
}

\examples{
library(PBtest)
data(sim1)

## attach list elements to the global environment
list2env(sim1,globalenv())

rho.hat <- getRho(YY=simdat, xx=Grp, id=Subj, weights=ww)
rho.all <- getRho(YY=simdat, xx=Grp, id=Subj, weights=ww, return.all=TRUE)
## Given the true correlation is 0.2, the mean square error is
mse <- mean((rho.all-0.2)^2)

## Sigma matrix
Sigma <- getSigma(YY=simdat, xx=Grp, id=Subj, weights=ww)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
