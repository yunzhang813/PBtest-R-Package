\name{PBtest}
\alias{PBtest}
\alias{PBmap}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
  PB-transformed test
}
\description{
  \code{PBmap} calculates the proposed transformations, i.e. B-map and P-map.

  \code{PBtest} conducts PB-transformed parametric or nonparametric test.
}

\usage{
PBmap(Delta, Sigma=NULL)
PBtest(XX, Delta, Sigma=NULL, test="fast.ttest", debug=FALSE, ...)
}

\arguments{
  \item{XX}{Numeric data matrix. Observations in rows and repetitions (e.g. genes) in columns.}
  \item{Sigma}{Covariance matrix. If \code{NULL}, identity matrix will be use.}
  \item{Delta}{Numeric vector of group indication, indicated by 1 and -1.}
  \item{test}{Parametric or nonparametric test. Options: \code{"fast.ttest"} (default), \code{"ttest"} and \code{"wilcox"}.}
  \item{debug}{If \code{TRUE}, it returns a long list of intermediate calculation results.}
  \item{...}{further arguments to be passed to or from \code{t.test()} or \code{wilcox.test()}.}
}

\value{
\code{PBtest} returns the p-values of the PB-transformed test.

If \code{debug=TRUE}, \code{PBtest} and \code{PBmap} return a long list of values, including:
  \item{p.value}{P-values.}
  \item{sigma.sq}{A scalar value. Dcompose the covariance matrix into a scale parameter and a shape parameter, such that \code{Sigma=sigma.sq*Sigma0}.}
  \item{Sigma0.inv}{Inverse of \code{Sigma0}.}
  \item{mu.hat}{The weighted average.}
  \item{XX.tilde}{The centered data.}
  \item{Sigma0.tilde}{Shape parameter for the centered data.}
  \item{LL}{Diagnal matrix of non-zero eigen-values of \code{Sigma0.tilde}.}
  \item{TT}{Matrix of corresponding eigen-vectors of \code{Sigma0.tilde}.}
  \item{BB}{The B-map.}
  \item{YY}{B-transformed data.}
  \item{cc}{The coefficient vector.}
  \item{QQ}{Q matrix of the QR decomposition.}
  \item{RR}{R matrix of the QR decomposition.}
  \item{Rot}{Rotation matrix by construction.}
  \item{PP}{The P-map.}
  \item{YY.tilde}{PB-transformed data.}
}

\note{
In \code{PBtest}, three test options are available. By default, \code{"fast.ttest"} conducts speed efficient t-tests for large-scale data, which uses \code{rwottests()} from the \code{genefilter} package. However, this fast t-test does not provide one-sided test. One-sided test can be specified using \code{(test="ttest", alternative="greater")} or \code{(test="wilcox", alternative="greater")}.
}

\references{
  \cite{Zhang, Y., Topham, D.J., Falsey, A.R. and Qiu, X. (2017). Highly Efficient Two-Group Comparison Methods for Correlated Observations with Heterogeneous Variance Structure. \emph{Submited}.}
}

\author{
Yun Zhang, Xing Qiu
}

\seealso{
  \code{\link[genefilter]{rowttests}}, \code{\link{t.test}}, \code{\link{wilcox.test}}.
}

\examples{
## Generate some data.
n1 <- 5; n2 <- 5; n <- n1+n2
XX <- matrix(rnorm(n), nrow=n)
Delta <- rep(c(1,-1), c(n1,n2))
## A general coveriance matrix.
weights <- runif(n, 1, 6)
corr <- matrix(0.2, n, n)
diag(corr) <- 1
Sigma <- diag(weights) \%*\% corr \%*\% diag(weights)

## PBmap.
mymap <- PBmap(Delta, Sigma)
mymap$BB  #matrix for B-map
mymap$PP  #matrix for P-map

## By default, two-sided t-test is performed. 
PBtest(XX, Delta, Sigma)

## One-sided t-test.
PBtest(XX, Delta, Sigma, test="ttest", alternative="greater")

## Two-sided rank-based test.
PBtest(XX, Delta, Sigma, test="wilcox")
## One-sided rank-based test.
PBtest(XX, Delta, Sigma, test="wilcox", alternativ="less")

## If Sigma=NULL (i.e. identity matrix), the parametric test reduces to the two-sample Student's t-test.
PBtest(XX[,1], Delta) #just use the first column of data
t.test(XX[Delta==1,1], XX[Delta==-1,1], var.equal=TRUE)$p.value

## debug
PBtest(XX, Delta, Sigma, debug=TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
