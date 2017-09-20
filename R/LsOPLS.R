
## Second implementation of the Light-sparse-OPLS (L-sOPLS) algorithm
## Application of the sPLS method to a deflated matrix obtained by OPLS
## Matrix without the Y-orthogonal components.
## Cross-validated two-step optimization for the numbers of orthogonal and
## predictive components.


#' @export LsOPLS
#'
#' @title Light-sparse-OPLS method for biomarker discovery
#'
#' @description
#' Implementation of the Light-sparse-OPLS algorithm
#' Application of the sPLS method to a deflated matrix obtained by OPLS
#' (matrix without the Y-orthogonal components)
#' Cross-validated two-step optimization for the numbers of orthogonal and predictive components.
#'
#' @param X A numeric matrix containing the spectra.
#' @param Y A numeric vector containing the target to predict.
#' @param Auto Argument (TRUE or FALSE). If TRUE, the first optimalized step of L-sOPLS is performed. If FALSE, No has to be chosen.
#' @param No The number (integer) of orthogonal components to extract (with a maximum of 9 orthogonal components). If Auto=TRUE, No must be such that No=NA.
#' @param cv The number of cross-validation segments (for LOO: cv = dim(X)[1] )
#' @param Np The number (integer) of predictive components in the L-sOPLS. Can be a sequence of values to test in the optimization step (ex: seq(1,5,by=1)).
#' @param Eta Parameter for the level of sparsity into the sPLS step. Can be a sequence of values to test in the optimization step (ex: seq(0.6,0.99,by=0.01)).
#' @param Method Parameter for the choice of the PLS algorithm ("pls2" or "simpls")
#'
#' @return A print of cross-validated optimal parameters, a grid of RMSEP values and the sPLS final step containing the list of final selected biomarkers:
#' @return The number of final orthogonal dimensions (No)
#' @return Cross-validation results (containing a MSPE matrix, the optimal Eta parameter $eta.opt and the optimal number of predictive dimensions $K.opt)
#' @return A RMSEP matrix
#' @return Sparse results with the final selected variables
#'
#' @examples
#' data("HumanSerum")
#' X <- as.matrix(HumanSerumSpectra)
#' Ex1 <- LsOPLS(X, ClassHS, Auto = TRUE, No = NA, cv = dim(X)[1],
#' Np = seq(1,3,by=1), Eta = seq(0.8,0.99,by=0.01), Method = "simpls")
#' Ex2 <- LsOPLS(X, ClassHS, Auto = FALSE, No = 2, cv = dim(X)[1],
#' Np = seq(1,3,by=1), Eta = seq(0.8,0.99,by=0.01), Method = "pls2")
#' print(Ex1)
#' print(Ex2)
#'
#' @importFrom ropls opls
#' @importFrom spls cv.spls
#' @importFrom spls spls
#'



LsOPLS = function(X, Y, Auto, No, cv, Np, Eta, Method){

# checks

if (missing(X)){
warning("X is missing with no default value")
}

if (missing(Y)){
warning("Y is missing with no default value")
}

if (missing(cv)){
warning("cv is missing with no default value")
}

if (missing(Np)){
warning("Np is missing with no default value")
}

if (missing(Eta)){
warning("Eta is missing with no default value")
}

if (missing(Method)){
warning("Method is missing with no default value")
}


if (! is.numeric(X)) {
  stop(deparse(substitute(X)), " is not numeric.")
}

if (! is.matrix(X)) {
  stop(deparse(substitute(X)), " is not a matrix.")
}

if (! is.numeric(Y)) {
    stop(deparse(substitute(Y)), " is not numeric.")
}

if (!(No %in% c(1:9,NA))) {
    stop(deparse(substitute(No)), " is not appropriate.")
}

if (! is.numeric(cv)) {
    stop(deparse(substitute(cv)), " is not numeric.")
}

if (! is.numeric(Np)) {
    stop(deparse(substitute(Np)), " is not numeric.")
}

if (! is.numeric(Eta)) {
    stop(deparse(substitute(Eta)), " is not numeric.")
}

if (TRUE %in% (Eta < 0)) {
    stop(deparse(substitute(Eta)), " is not correct.")
}

if (TRUE %in% (Eta > 1)) {
    stop(deparse(substitute(Eta)), " is not correct.")
}

if (!(Auto %in% c('TRUE','FALSE'))) {
    stop(deparse(substitute(Auto)), " is not appropriate.")
}

if (!(Method %in% c('pls2','simpls'))) {
    stop(deparse(substitute(Method)), " is not appropriate.")
}



################

#Automatic OPLS:
if ((Auto == 'TRUE')&!(No %in% c(1:9))) {
oplscv1 <- ropls::opls(X, Y,  predI = 1, orthoI = No, crossvalI = cv, scaleC="center", plotL = FALSE)
No <- dim(oplscv1@modelDF)[1] - 3
#Deflated matrix:
start1 <- OPLSDA(as.matrix(X), Y, impT = FALSE, impG = FALSE, no = No, nb = dim(X)[2], out.path = ".")
Xopls1 <- start1$Xopls
#sPLS step:
CV1 <- spls::cv.spls(Xopls1, Y, fold = cv, K = Np, eta = Eta, kappa=0.5, select=Method, fit="simpls",
scale.x=FALSE, scale.y=FALSE, plot.it=TRUE)
#Optimal L-sOPLS:
Opt1 <- spls::spls(Xopls1, Y, K = CV1$K.opt, eta = CV1$eta.opt, kappa=0.5, select=Method, fit="simpls",
scale.x=FALSE, scale.y=FALSE, eps=1e-4, maxstep=100, trace=TRUE)

output <- list("Number of orthogonal dimensions" = No, "Cross-validation results" = CV1, "RMSEP matrix" = sqrt(CV1$mspemat), "Sparse results" = Opt1)
}


else if ((Auto == 'TRUE')&!(No == 'NA')) {
warning("With Auto = TRUE, No has to be optimally chosen by the first L-sOPLS optimization step, i.e. No = NA")

output <- list()
}



#Manual OPLS
if ((Auto == 'FALSE')&(No %in% c(1:9))) {
oplscv2 <- ropls::opls(X, Y,  predI = 1, orthoI = No, crossvalI = cv, scaleC="center", plotL = FALSE)
#Deflated matrix:
start2 <- OPLSDA(as.matrix(X), Y, impT = FALSE, impG = FALSE, no = No, nb = dim(X)[2], out.path = ".")
Xopls2 <- start2$Xopls
#sPLS step:
CV2 <- spls::cv.spls(Xopls2, Y, fold = cv, K = Np, eta = Eta, kappa=0.5, select=Method, fit="simpls",
scale.x=FALSE, scale.y=FALSE, plot.it=TRUE)
#Optimal L-sOPLS:
Opt2 <- spls::spls(Xopls2, Y, K = CV2$K.opt, eta = CV2$eta.opt, kappa=0.5, select=Method, fit="simpls",
scale.x=FALSE, scale.y=FALSE, eps=1e-4, maxstep=100, trace=TRUE)

output <- list("Number of orthogonal dimensions" = No, "Cross-validation results" = CV2, "RMSEP matrix" = sqrt(CV2$mspemat), "Sparse results" = Opt2)
}


else if ((Auto == 'FALSE')&!(No %in% c(1:9))) {
warning("No is needed for manual L-sOPLS and is missing with no default value")

output <- list("Number of orthogonal dimensions" = No)
}



return(output)

}










