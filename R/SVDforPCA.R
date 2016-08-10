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
#'   \item{\code{pcs}}{Scores}
#'   \item{\code{pcu}}{Normalized scores}
#'   \item{\code{pcv}}{Loadings}
#'   \item{\code{pcd}}{Xingular values}
#'   \item{\code{var}}{Explained variance}
#'   \item{\code{cumvar}}{Cumulated explained variance}
#'   \item{\code{eigval}}{Eigenvalues}
#' }
#'
#' @examples

#' data("iris")
#' data = iris[,1:4]
#' class = as.numeric(iris[,5])
#' PCA.res = SVDforPCA(data)
#'
#' RowNames = row.names(data)
#'
#' Xax=1 # composante x
#' Yax=2 # composante y
#' Xlim=c(min(PCA.res$pcs[,Xax])*1.4, max(PCA.res$pcs[,Xax])*1.4)
#' Ylim=c(min(PCA.res$pcs[,Yax])*1.4, max(PCA.res$pcs[,Yax])*1.4)
#'
#' plot(PCA.res$pcs[,Xax],PCA.res$pcs[,Yax], col=class,
#'     xlab=paste0("PC",Xax," (", PCA.res$var[Xax] ,"%)"), xlim=Xlim,
#'     ylab=paste0("PC",Yax," (", PCA.res$var[Yax] ,"%)"), ylim=Ylim,
#'     main=paste0("PCA score plot"))
#' text(PCA.res$pcs[,Xax],PCA.res$pcs[,Yax],labels=RowNames, pos=c(2,3), cex = 0.7, col=class)

#' abline(h=0, v=0)
#'

SVDforPCA= function(x, ncomp=min(dim(x))) {



# column centering of X
x=x-matrix(apply(x,2,mean),nrow=dim(x)[1],ncol=dim(x)[2],byrow=T) # centring of x over the columns

# SVD of X
x.svd = svd(x) # Compute the singular-value decomposition
x.scores = x.svd$u %*% diag(x.svd$d) # scores
x.normscores = x.svd$u # normalised scores
x.loadings = x.svd$v # loadings
x.singularval = x.svd$d # singular values

# X-Variance explained
x.vars = x.singularval^2 / (nrow(x) - 1)
x.eigval = x.singularval^2
x.totalvar = sum(x.vars)
x.relvars = x.vars / x.totalvar

x.variances = 100 * round(x.relvars, digits = 3) # variance
x.cumvariances = cumsum(x.variances) # cumulative variance

# as matrix

x.scores = as.matrix(x.scores)
x.normscores = as.matrix(x.normscores)
x.loadings = as.matrix(x.loadings)

# selection of the first n components

x.scores = x.scores[,1:ncomp]
x.normscores = x.normscores[,1:ncomp]
x.loadings = x.loadings[,1:ncomp]
x.singularval = x.singularval[1:ncomp]

x.variances = x.variances[1:ncomp]
x.cumvariances = x.cumvariances[1:ncomp]



res=list(pcs=x.scores, pcu=x.normscores, pcv=x.loadings, eigval=x.eigval, pcd=x.singularval, var=x.variances, cumvar=x.cumvariances)

return(res)

}


