\name{GPL}
\alias{GPL}
\docType{data}
\title{
example of an typical real individual Global Peak List
}
\description{
Matrix containing a real typical individual Global Peak List, with coordinates, binary presences and intensities
}
\usage{data("GPL")}
\format{
  A data frame with 1730 rows and 51 columns (ID, two coordinates, 24 columns of intensities and 24 columns of binary presence/absence).
  }

\source{
Nantes (E. Martineau), used during the UCL-ULg-Nantes repetition project
}

\examples{
Int = GPL[,seq(3,dim(GPL)[2]- 1, 2)]
Coord = GPL[,c(1,2)]
IntT <- pTreatGPL(Int, Coord, 4.63, 5.1)
IntT[,1:2]
}
\keyword{datasets}
