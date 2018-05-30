\name{sim1}
\alias{sim1}
\docType{data}

\title{
A set of simulated data
}

\description{
Simulated data of 2000 genes (1000 null genes and 1000 differentially expressed genes) and 40 samples. Paired samples are collected from 20 subjects. These samples come from two groups. Subject IDs and weights are available.
}

\usage{data("sim1")}

\format{
  The object \code{sim1} is a list with the following elements:
  \describe{
    \item{\code{simdat}}{Data matrix with samples in rows and genes in columns. Column names indicate whether the gene is a null gene (e.g. \code{null.G1}) or a differentially expressed gene (e.g. \code{sig.G1}).}
    \item{\code{Grp}}{Vector of group labels.}
    \item{\code{ww}}{Vector of weights.}
    \item{\code{Subj}}{Vector of subject IDs, which are also the row names of \code{simdat}.}
  }
}

\references{
  \cite{Yun Zhang, Gautam Bandyopadhyay, and Xing Qiu. (2018) Highly efficient hypothesis testing methods for regression-type tests with correlated observations and heterogeneous variance structure. Submitted.}
}

\examples{
library(PBtest)
data(sim1)

## attach list elements to the global environment
list2env(sim1,globalenv())

## the PB-transformed t-test, with Sigma estimation
pvec.PB <- PBtest(YY=simdat, xx=Grp, test="t", id=Subj, weights=ww)$p.value
}

\keyword{datasets}
