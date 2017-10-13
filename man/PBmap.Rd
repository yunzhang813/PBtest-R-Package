\name{PBmap}
\alias{PBmap}
\alias{PBtest}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
  Two-group comparison test with PB-transformation
}
\description{
  \code{PBmap()} calculates the proposed data transformations, i.e. B-map and P-map, which transform the original data and hypothesis testing problem to the equivalent form that can be tested by standard tests.

  \code{PBtest()} takes in the original data and conducts parametric t-test or nonparametric Wilcoxon (rank-based) test for the transformed data from \code{PBmap()}.
}

\usage{
PBmap(XX, Sigma, Delta)
PBtest(XX, Sigma, Delta, test="t", ...)
}

\arguments{
  \item{XX}{Data vector of observations combined from both groups.}
  \item{Sigma}{Covariance matrix for combined data.}
  \item{Delta}{Numeric vector for groups: one group indicated by 1 and the other group indicated by -1.}
  \item{test}{Parametric or nonparametric test. Options: \code{"t"} (default) and \code{"wilcox"}.}
  \item{...}{further arguments to be passed to or from test functions.}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
\code{PBmap()} returns a list of values:
  \item{sigma.sq}{A scalar value such that \code{sigma.sq*Sigma0=Sigma}.}
  \item{Sigma0}{Standardized structure matrix.}
  \item{mu.hat}{Estimated weighted average.}
  \item{XX.tilde}{Centered data.}
  \item{Sigma0.tilde}{Strucure matrix for \code{XX.tilde}.}
  \item{LL}{Diagnal matrix of non-zero eigen-values of \code{Sigma0.tilde}.}
  \item{TT}{Matrix of corresponding eigen-vectors of \code{Sigma0.tilde}.}
  \item{BB}{The B-map.}
  \item{YY}{Transformed data after B-map.}
  \item{cc}{Coefficient vector of the transformed alternative hypothesis after B-map.}
  \item{QQ}{Q matrix of the QR decomposition.}
  \item{RR}{R matrix of the QR decomposition.}
  \item{Rot}{Rotation matrix by construction.}
  \item{PP}{The P-map.}
  \item{YY.tilde}{Final transformed data.}
  
\code{PBtest()} returns the p-value for the overall test.
}

\references{
%% ~put references to the literature/web site here ~
}

\author{
Yun Zhang, Xing Qiu
}

\seealso{
  \link{t.test}, \link{wilcox.test}.
}

\examples{
## Generate some data.
n1 <- 5; n2 <- 5; n <- n1+n2
XX <- rnorm(n)
Delta <- rep(c(1,-1), c(n1,n2))
## A general coveriance matrix.
weights <- runif(n, 1, 6)
corr <- matrix(0.2, n, n)
diag(corr) <- 1
Sigma <- diag(weights) \%*\% corr \%*\% diag(weights)
## An identity matrix.
Sigma.I <- diag(n)


## Data transformation.
mydat <- PBmap(XX, Sigma, Delta)
mydat$YY.tilde  #the transformed data
mydat$BB  #matrix for B-map
mydat$PP  #matrix for P-map

## By default, two-sided t-test is performed. 
PBtest(XX, Sigma, Delta)

## One-sided t-test.
PBtest(XX, Sigma, Delta, alternative="greater")

## Two-sided rank-based test.
PBtest(XX, Sigma, Delta, test="wilcox")

## If Sigma.I (i.e. identity covariance matrix), compare with standard tests.
PBtest(XX, Sigma.I, Delta)
t.test(XX[Delta==1], XX[Delta==-1])$p.value #only numerical difference
PBtest(XX, Sigma.I, Delta, test="wilcox")
wilcox.test(XX[Delta==1], XX[Delta==-1])$p.value #difference in p-value

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
