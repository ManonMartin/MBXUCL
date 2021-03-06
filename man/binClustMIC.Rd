% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binClustMIC.R
\name{binClustMIC}
\alias{binClustMIC}
\title{Unsupervised clustering on (GPL) binary presence/absence of peaks and associated MIC indexes}
\usage{
binClustMIC(Positions, Distance, nClust, Trcl, Dendr = TRUE)
}
\arguments{
\item{Positions}{A binary numeric matrix}

\item{Distance}{Choice of the distance measure adapted to binary objects ('Jaccard' or 'Ochiai')}

\item{nClust}{The number of groups to retrieve (donors, mixtures, ...).}

\item{Trcl}{The real groups' memberships of the samples, true class.}

\item{Dendr}{Logical argument (TRUE/FALSE) to obtain graphical dendrogram based on the Ward algorithm.}
}
\value{
A list of MIC quality indexes (Dunn, Davies-Bouldin, Rand and adjusted-Rand):
\describe{
  \item{\code{DunnW}}{Dunn index for Ward clustering}
  \item{\code{DBW}}{Davies-Bouldin index for Ward clustering}
  \item{\code{RandW}}{Rand index for Ward clustering}
  \item{\code{AdjRandW}}{Adjusted Rand index for Ward clustering}
   }
}
\description{
Unsupervised clustering on (GPL) binary presence/absence of peaks based on Ward algorithm.
Calculation of MIC statistical criteria of clustering quality:
Dunn, Davies-Bouldin, Rand and adjusted-Rand indexes.
}
\examples{
Pos = t(GPL[, seq(4,dim(GPL)[2], 2)])
Trcl = c(rep(1,8), rep(2,8), rep(3,8))
binClustMIC(Positions = Pos, Distance = 'Jaccard', nClust = 3, Trcl, Dendr = TRUE)

}
\author{
Baptiste Feraud
}
