#' @export SVDforPCA
#' @title Singular value decomposition for PCA analysis
#'
#' @description
#' PCA over a X matrix by singular value decomposition, the preprocessing involves the mean-centering of X.
#'
#' @param x   A data matrix on which will be based the analysis.
#' @param ncomp  Number of Principal Components.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{scores}}{Scores}
#'   \item{\code{loadings}}{Loadings}
#'   \item{\code{eigval}}{Eigenvalues}
#'   \item{\code{pcd}}{Singular values}
#'   \item{\code{pcu}}{Normalized scores}
#'   \item{\code{var}}{Explained variance}
#'   \item{\code{cumvar}}{Cumulated explained variance}
#'  \item{\code{original.dataset}}{Original dataset}
#'
#' }
#'
#' @examples
#' data('HumanSerum')
#' PCA.res = SVDforPCA(HumanSerumSpectra)


SVDforPCA <- function(x, ncomp = min(dim(x))) {


  original.dataset <- x
  # column centering of X
  x <- x - matrix(apply(x, 2, mean), nrow = dim(x)[1], ncol = dim(x)[2], byrow = TRUE)  # centring of x over the columns

  # SVD of X
  x.svd <- svd(x)  # Compute the singular-value decomposition
  x.scores <- x.svd$u %*% diag(x.svd$d)  # scores
  x.normscores <- x.svd$u  # normalised scores
  x.loadings <- x.svd$v  # loadings
  x.singularval <- x.svd$d  # singular values
  names(x.singularval) <- paste0("PC", 1:length(x.singularval))

  # X-Variance explained
  x.vars <- x.singularval^2/(nrow(x) - 1)
  x.eigval <- x.singularval^2
  names(x.eigval) <- paste0("PC", 1:length(x.eigval))

  x.totalvar <- sum(x.vars)
  x.relvars <- x.vars/x.totalvar

  x.variances <- 100 * x.relvars  # variance
  names(x.variances) <- paste0("PC", 1:length(x.variances))

  x.cumvariances <- cumsum(x.variances)  # cumulative variance
  names(x.cumvariances) <- paste0("PC", 1:length(x.cumvariances))

  # as matrix

  x.scores <- as.matrix(x.scores)
  dimnames(x.scores) <- list(rownames(x), paste0("PC", 1:ncol(x.scores)))

  x.normscores <- as.matrix(x.normscores)
  dimnames(x.normscores) <- list(rownames(x), paste0("PC", 1:ncol(x.normscores)))

  x.loadings <- as.matrix(x.loadings)

  dimnames(x.loadings) <- list(colnames(x), paste0("PC", 1:ncol(x.loadings)))


  # selection of the first n components

  x.scores <- x.scores[, 1:ncomp]
  x.normscores <- x.normscores[, 1:ncomp]
  x.loadings <- x.loadings[, 1:ncomp]
  x.singularval <- x.singularval[1:ncomp]

  x.variances <- x.variances[1:ncomp]
  x.cumvariances <- x.cumvariances[1:ncomp]



  res <- list(scores = x.scores, loadings = x.loadings, eigval = x.eigval, pcu = x.normscores,
              pcd = x.singularval, var = x.variances, cumvar = x.cumvariances,
              original.dataset = original.dataset)

  return(res)

}


