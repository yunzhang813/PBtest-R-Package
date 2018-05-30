\name{PBtest}
\alias{PBtest}
\alias{PBmap}

\title{
  PB-transformed test
}
\description{
  \code{PBmap} calculates the proposed transformations (see reference), i.e. B-map and P-map.

  \code{PBtest} conducts the PB-transformed parametric or nonparametric test (see reference).
}

\usage{
PBmap(xx, Sigma)
PBtest(YY, xx, Sigma=NULL, test="t", id=NULL, weights=NULL, debug=FALSE)
}

\arguments{
  \item{YY}{Vector or matrix of response values. If gene expression matrix, genes in columns and subjects in rows.}
  \item{xx}{Vector of covariate values.}
  \item{Sigma}{Vairance-covariance matrix.}
  \item{test}{Parametric or nonparametric test. Options: \code{"t"} (default), \code{"wilcox"}.}
  \item{id}{Vector of subject IDs (for repeated measurements).}
  \item{weights}{Vector of \emph{a priori} weights.}
  \item{debug}{Logical value. If \code{TRUE}, intermediate results are returned.}
}

\value{
A list of testing results:
  \item{p.value}{P-value(s) of two-sided test.}
  \item{statistic}{Test statistic(s).}
  \item{df}{Adjusted degrees of freedom by \code{KRapprox}.}

If \code{debug=TRUE}, the following intermediate results are also returned.
  \item{YY.tilde}{PB-transformed data.}
  \item{Sigma}{Variance-covariance matrix.}
  \item{sigma.sq}{A scalar value. Dcompose the covariance matrix into a scale parameter and a shape parameter, such that \code{Sigma=sigma.sq*SS}.}
  \item{SS.inv}{Inverse of \code{SS}.}
  \item{SS.tilde}{Shape parameter after centering, i.e. \code{SS - JJ}, where \code{JJ} is the matrix of all 1's.}
  \item{LL}{Diagnal matrix of non-zero eigen-values of \code{SS.tilde}.}
  \item{TT}{Matrix of corresponding eigen-vectors of \code{SS.tilde}.}
  \item{BB}{B-map.}
  \item{zz}{Transformed covariate.}
  \item{QQ}{Q matrix of the QR decomposition.}
  \item{RR}{R matrix of the QR decomposition.}
  \item{Rot}{Rotation matrix by construction.}
  \item{PP}{P-map.}
}

\note{
Depending on different levels of \emph{a priori} information, the following options are available.

1. If \code{Sigma} is known, set \code{Sigma=Sigma} and leave \code{id=NULL, weights=NULL}.

2. If \code{Sigma} is unknown, estimate using \code{getSigma}. In the arguments here, set \code{Sigma=NULL, id=id, weights=weights} or \code{NULL}. If \code{weights=NULL}, equal weights of 1 will be assigned.

3. If \code{Sigma=NULL, id=NULL, weights=NULL}, the identity matrix will be used for \code{Sigma}.
}

\references{
  \cite{Yun Zhang, Gautam Bandyopadhyay, and Xing Qiu. (2018) Highly efficient hypothesis testing methods for regression-type tests with correlated observations and heterogeneous variance structure. Submitted.}
}

\author{
Yun Zhang, Xing Qiu
}

\seealso{
  \code{\link[PBtest]{getSigma}}, \code{\link[PBtest]{KRapprox}}.
}

\examples{
library(PBtest)
data(sim1)

## attach list elements to the global environment
list2env(sim1,globalenv())

## the PB-transformed t-test, with Sigma estimation
out.PB <- PBtest(YY=simdat, xx=Grp, test="t", id=Subj, weights=ww)
## p-values
pvec.PB <- out.PB$p.value
## degrees of freedom
df.PB <- out.PB$df

## estimated correlation
rho.hat.PB <- getRho(YY=simdat, xx=Grp, id=Subj, weights=ww)

## Sigma matrix
Sigma <- getSigma(YY=simdat, xx=Grp, id=Subj, weights=ww)
## or
Sigma <- PBtest(YY=simdat, xx=Grp, id=Subj, weights=ww, debug=TRUE)$Sigma

## the estimated correlation is also the off-diagonal non-zero value of the following matrix
cor.mat <- diag(sqrt(ww)) \%*\% Sigma \%*\% diag(sqrt(ww))

## the transformation maps
mymap <- PBmap(xx=Grp, Sigma)
Bmap <- mymap$BB
Pmap <- mymap$PP

## the PB-transformed data
simdat.PB <- PBtest(YY=simdat, xx=Grp, id=Subj, weights=ww, debug=TRUE)$YY.tilde
dim(simdat.PB) #this should have one less row than the original data dim(simdat)

## PBtest for two-sample t-test
out1.PB <- PBtest(YY=simdat, xx=Grp)
pvec1.PB <- out1.PB$p.value
df1.PB <- out1.PB$df
## two-sample t-test
pvec.t <- genefilter::colttests(simdat, Grp)$p.value
summary(pvec1.PB-pvec.t) #numerially same
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
