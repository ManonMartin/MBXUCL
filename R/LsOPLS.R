
## Implementation of the Light-sparse-OPLS algorithm Application of the sPLS
## method to a deflated matrix obtained by OPLS Matrix without the Y-orthogonal
## components


#' @export LsOPLS
#' @title Light-sparse-OPLS method for biomarker discovery
#' @description
#' Implementation of the Light-sparse-OPLS algorithm
#' Application of the sPLS method to a deflated matrix obtained by OPLS
#' (matrix without the Y-orthogonal components)
#'
#' @param X A numeric matrix containing the spectra.
#' @param Y A numeric vector containing the target to predict.
#' @param No The number of orthogonal components to extract (may be optimalized previously).
#' @param Eta Parameter for the level of sparsity into the sPLS method
#' @param Method Parameter for the choice of the PLS algorithm ('pls2' or 'simpls')
#'
#' @return A print of the sPLS final step containing the list of final selected biomarkers.
#'
#' @examples
#' data('iris')
#' X <- as.matrix(iris[,1:4])
#' Y <- as.numeric(iris[,5])
#' LsOPLS(X, Y, 2, 0.9, 'simpls')
#'
#' @importFrom spls spls
#'



LsOPLS <- function(X, Y, No, Eta, Method) {

  # checks

  if (missing(X))  {
    warning("X is missing with no default value")
  }

  if (missing(Y)) {
    warning("Y is missing with no default value")
  }

  if (missing(No)) {
    warning("No is missing with no default value")
  }

  if (missing(Eta))  {
    warning("Eta is missing with no default value")
  }

  if (missing(Method)) {
    warning("Method is missing with no default value")
  }


  if (!is.numeric(X)) {
    stop(deparse(substitute(X)), " is not numeric.")
  }

  if (!is.matrix(X)) {
    stop(deparse(substitute(X)), " is not a matrix.")
  }

  if (!is.numeric(Y)) {
    stop(deparse(substitute(Y)), " is not numeric.")
  }

  if (!is.numeric(No)) {
    stop(deparse(substitute(No)), " is not numeric.")
  }

  if (!is.numeric(Eta)) {
    stop(deparse(substitute(Eta)), " is not numeric.")
  }

  if (Eta < 0) {
    stop(deparse(substitute(Eta)), " is not correct.")
  }

  if (Eta > 1) {
    stop(deparse(substitute(Eta)), " is not correct.")
  }

  if (!(Method %in% c("pls2", "simpls"))) {
    stop(deparse(substitute(Method)), " is not appropriate.")
  }


  #### Deflated matrix with the y-orthogonal component removed ####


  # Debut de l'algorithme de l'OPLS-DA avec le calcul des poids de X : 1.
  w <- (t(X) %*% Y) %*% solve(t(Y) %*% Y)  # X-weights w
  w <- as.vector(w)
  # 2.
  norm.w <- sqrt(w %*% w)  # norm of w (vector)
  w <- w %*% solve(norm.w)  # normalised w


  # creation de la boucle pour extraire les composantes orthogonales

  for (j in 1:No) {
    # 3.
    t <- (X %*% w) %*% solve(t(w) %*% w)  # X-scores
    # 4.
    c <- t(Y) %*% t %*% solve(t(t) %*% t)  # y-weights
    # 5.
    u <- Y %*% c %*% solve(t(c) %*% c)  # y-scores
    # 6.
    p <- (t(X) %*% t) %*% solve(t(t) %*% t)  # 1st PLS loading
    p <- as.vector(p)
    norm.p <- sqrt(p %*% p)
    # 7.
    wo <- p - (w %*% (t(w) %*% p) %*% solve(t(w) %*% w))  # Xortho-weights
    wo <- as.vector(wo)
    # 8.
    norm.wo <- sqrt(wo %*% wo)
    wo <- wo %*% solve(norm.wo)  # normalised wortho
    cv <- norm.wo %*% solve(norm.p)  # number of orthogonal components to keep
    # 9.
    to <- X %*% wo %*% solve(t(wo) %*% wo)  # Xortho-scores
    # 10.
    po <- t(X) %*% to %*% solve(t(to) %*% to)  # ortho. loading

    varortho <- to %*% t(po)
    var <- t %*% t(p)
    # 11.
    XDefl <- X - to %*% t(po)
  }

  #### sPLS ####

  object <- spls::spls(XDefl, Y, K = No, eta = Eta, kappa = 0.5, select = Method,
                      fit = "widekernelpls", scale.x = FALSE, scale.y = FALSE,
                      eps = 1e-04, maxstep = 100, trace = FALSE)

  print(object)

}
