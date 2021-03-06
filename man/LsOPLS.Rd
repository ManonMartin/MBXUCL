% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LsOPLS.R
\name{LsOPLS}
\alias{LsOPLS}
\title{Light-sparse-OPLS method for biomarker discovery}
\usage{
LsOPLS(X, Y, Auto, No, cv, Np, Eta, Method)
}
\arguments{
\item{X}{A numeric matrix containing the spectra.}

\item{Y}{A numeric vector containing the target to predict.}

\item{Auto}{Argument (TRUE or FALSE). If TRUE, the first optimalized step of L-sOPLS is performed. If FALSE, No has to be chosen.}

\item{No}{The number (integer) of orthogonal components to extract (with a maximum of 9 orthogonal components). If Auto=TRUE, No must be such that No=NA.}

\item{cv}{The number of cross-validation segments (for LOO: cv = dim(X)[1] )}

\item{Np}{The number (integer) of predictive components in the L-sOPLS. Can be a sequence of values to test in the optimization step (ex: seq(1,5,by=1)).}

\item{Eta}{Parameter for the level of sparsity into the sPLS step. Can be a sequence of values to test in the optimization step (ex: seq(0.6,0.99,by=0.01)).}

\item{Method}{Parameter for the choice of the PLS algorithm ("pls2" or "simpls")}
}
\value{
A print of cross-validated optimal parameters, a grid of RMSEP values and the sPLS final step containing the list of final selected biomarkers:

The number of final orthogonal dimensions (No)

Cross-validation results (containing a MSPE matrix, the optimal Eta parameter $eta.opt and the optimal number of predictive dimensions $K.opt)

A RMSEP matrix

Sparse results with the final selected variables
}
\description{
Implementation of the Light-sparse-OPLS algorithm
Application of the sPLS method to a deflated matrix obtained by OPLS
(matrix without the Y-orthogonal components)
Cross-validated two-step optimization for the numbers of orthogonal and predictive components.
}
\examples{
data("HumanSerum")
X <- as.matrix(HumanSerumSpectra)
Ex1 <- LsOPLS(X, ClassHS, Auto = TRUE, No = NA, cv = dim(X)[1],
Np = seq(1,3,by=1), Eta = seq(0.8,0.99,by=0.01), Method = "simpls")
Ex2 <- LsOPLS(X, ClassHS, Auto = FALSE, No = 2, cv = dim(X)[1],
Np = seq(1,3,by=1), Eta = seq(0.8,0.99,by=0.01), Method = "pls2")
print(Ex1)
print(Ex2)

}
