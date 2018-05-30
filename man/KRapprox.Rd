\name{KRapprox}
\alias{KRapprox}

\title{
Kenward-Roger approximation to the degrees of freedom
}

\description{
This implementation is based on the \code{pbkrtest} package. See reference.
}

\usage{
KRapprox(xx, Sigma, L=c(0,1))
}

\arguments{
  \item{xx}{Vector of covariate values.}
  \item{Sigma}{Vairance-covariance matrix.}
  \item{L}{Linear contrast. Default: \code{L=c(0,1)}.}
}

\value{
A list with components:
  \item{df}{Adjusted degrees of freedom.}
  \item{scaling}{Scaling factor for the F-statistic. Here always equal to 1.}
}

\references{
  \cite{Halekoh, U., & Højsgaard, S. (2014). A kenward-roger approximation and parametric bootstrap methods for tests in linear mixed models–the R package pbkrtest. Journal of Statistical Software, 59(9), 1-30.}
  
  \cite{Kenward, M. G., & Roger, J. H. (1997). Small sample inference for fixed effects from restricted maximum likelihood. Biometrics, 983-997.}
}

\author{
Yun Zhang, Xing Qiu
}

\seealso{
For more details, please see also the \code{pbketest} package.
}

\examples{
## Consider two groups: A and B
n <- 20
xx <- rep(c("A","B"),n)

## e.g. degrees of freedom for two-sample t-test, i.e. covariance is the identity matrix
KRapprox(xx, diag(2*n)) #df = 2n-2

## e.g. degrees of freedom for paired t-test, i.e. block-diagnal covariance matrix
library(Matrix)
rho <- 0.5 #this may be any non-zero correlation coefficient
mat <- matrix(rho,2,2)
diag(mat) <- 1
Sigma <- as.matrix(bdiag(replicate(n,mat,simplify=FALSE)))
KRapprox(xx, Sigma) #df = n-1
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
