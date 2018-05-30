\name{HER2}
\alias{HER2}
\docType{data}

\title{
Example data of HER2-positive breast cancer
}

\description{
We downloaded a set of RNA-seq gene expression data from The Cancer Genome Atlas (TCGA). The data are sequenced on the Illumina GA platform with tissues collected from breast cancer subjects.  In particular, we selected 28 samples from the tissue source site ``BH'', which are controlled for white female subjects with the HER2-positive biomarkers. After data preprocessing based on nonspecific filtering, a total number of 11,453 genes are kept for subsequent analyses.  Among these data are 12 pairs of matched tumor and normal samples, and 4 unmatched tumor samples. The sequencing depths of the selected samples range from 23.80 million reads to 76.08 million reads. In the proposed analysis, we inversely weigh the sequencing samples by their sequencing depths (per million reads).
}

\usage{data("HER2")}

\format{
  The object \code{HER2} is a list with the following elements:
  \describe{
    \item{\code{dat}}{Data matrix with samples in rows and genes in columns. Column names are gene symbols.}
    \item{\code{Grp}}{Vector of group labels.}
    \item{\code{seq.depth}}{Vector of sequencing depth, corresponding to the rows of \code{dat}.}
    \item{\code{Subj}}{Vector of subject IDs, corresponding to the rows of \code{dat}.}
  }
}

\references{
  The Cancer Genome Atlas (TCGA) \url{https://cancergenome.nih.gov/}.

  \cite{Yun Zhang, Gautam Bandyopadhyay, and Xing Qiu. (2018) Highly efficient hypothesis testing methods for regression-type tests with correlated observations and heterogeneous variance structure. Submitted.}
}

\examples{
library(PBtest)
data(HER2)

## attach list elements to the global environment
list2env(HER2,globalenv())

## the PB-transformed t-test, with Sigma estimation
pvec.PB <- PBtest(YY=dat, xx=Grp, test="t", id=Subj, weights=seq.depth)$p.value
}

\keyword{datasets}
