#' @export UnsupClustering
#' @title Unsupervised classification of spectra
#' @description Runs Ward and K-means clustering algorithms and computes quality indexes from these clustering
#'
#' @param x  Datamatrix of spectra
#' @param RowNames  Optional, give the rownames of x
#' @param nClust  Number of clusters
#' @param y  The true class of spectra (as numeric)
#' @param print  If \code{TRUE}, prints the clustering results
#' @param impG  If \code{TRUE}, will save graphics (pdf format)
#' @param out.path  Path to save the graphs
#'
#'@return A list with the following elements:
#' \describe{
#'   \item{\code{m}}{ Number of columns of x}
#'   \item{\code{DI.Ward}}{Dunn index for Ward clustering}
#'   \item{\code{DI.Kmeans}}{Dunn index for Kmeans clustering}
#'   \item{\code{DBI.Wald}}{Davies-Bouldin index for Ward clustering}
#'   \item{\code{DBI.Kmeans}}{Davies-Bouldin index for Kmeans clustering}
#'   \item{\code{RI.Ward}}{Rand index for Ward clustering}
#'   \item{\code{RI.Kmeans}}{Rand index for Kmeans clustering}
#'   \item{\code{ARI.Ward}}{Adjusted Rand index for Ward clustering}
#'   \item{\code{ARI.Kmeans}}{Adjusted Rand index for Kmeans clustering}
#'    }
#'
#' @examples
#' data("iris")
#' data = as.matrix((iris[,1:4]))
#' class = as.numeric(iris[,5])
#'
#'Res = UnsupClustering(x = data, nClust = 3, y = class)
#'
#'
#'@importFrom clValid dunn
#'@importFrom clusterSim index.DB
#'@importFrom phyclust RRand


UnsupClustering = function(x, RowNames=NULL, nClust=3, y, print = FALSE ,impG=FALSE, out.path=".") {

# checks

  if (! is.numeric(x)) {
    stop(deparse(substitute(x)), " is not numeric.")
  }


  if (!is.null(RowNames)) {
    RowNames = as.character(RowNames)
    if (!length(RowNames) == dim(x)[1]) {
      stop("length of ", deparse(substitute(RowNames))," (",length(RowNames),")" , " is not equal to the row number of x ","(",dim(x)[1],")")
    }
  }


  if (! is.numeric(nClust)) {
    stop(deparse(substitute(nClust)), " is not numeric.")
  }  else if (length(nClust)>1) {
    stop(deparse(substitute(nClust)), " has a length > 1.")
  }

  if (! is.numeric(y)) {
    stop(deparse(substitute(y)), " is not numeric.")
  }


  if (! is.logical(print)) {
    stop(deparse(substitute(print)), " is not logical.")
  }

  if (! is.logical(impG)) {
    stop(deparse(substitute(impG)), " is not logical.")
  }

  if (! is.character(out.path)) {
    stop(deparse(substitute(print)), " is not character.")
  } else if (length(out.path)>1) {
    stop(deparse(substitute(out.path)), " has a length > 1.")
  }


m = dim(x)[2]


# WARD - hierarchical clustering
###################################

EuclDist <- dist(x, method ="euclidean", diag=FALSE, upper=FALSE)

WardClassif <- hclust(d=EuclDist,method="ward.D2")

CutTree <-cutree(tree=WardClassif, k=nClust) # class assignment
# cbind(y,CutTree)



# KMEANS:
###################################

KMeans <- kmeans(x, centers=nClust, iter.max = 100, nstart = 20)
KMeans$cluster
KMeans$iter


# Analysis of spectra quality
#__________________________________________________________

################################
##  Homogeneite des groupes   ##
################################


# DUNN INDEX
################
#Index to be maximized, between O and +inf
DIWard = round(dunn(distance=EuclDist, clusters=CutTree),5)
DIKmeans = round(dunn(distance=EuclDist, clusters=KMeans$cluster),5)




# DAVIES-BOULDIN INDEX
#######################
# to be minimized

DBIWald <- round(index.DB(x=x, cl=CutTree, d=EuclDist, centrotypes="medoids")$DB,5) # ou medoids
DBIKmeans <- round(index.DB(x=x, cl=KMeans$cluster, d=EuclDist, centrotypes="medoids")$DB,5)




################################
##   Qualite du clustering    ##
################################

#ADJUSTED RAND INDEX
#######################
# between 0 and 1


#True class

if(0 %in% unique(y)) {
  nlev=length(unique(y))
  trcl1 <- as.numeric(factor(y, levels=unique(y), labels=1:nlev))
}else {trcl1 <- as.numeric(y)}



RIWard <- round(as.numeric(RRand(trcl1, CutTree, lab = NULL)[1]),5)
ARIWard <- round(as.numeric(RRand(trcl1, CutTree, lab = NULL)[2]),5)

RIKmeans <- round(as.numeric(RRand(trcl1, KMeans$cluster, lab = NULL)[1]),5)
ARIKmeans <- round(as.numeric(RRand(trcl1, KMeans$cluster, lab = NULL)[2]),5)



## RESULTS


res=list(m=m, DI.Ward = round(DIWard,3), DI.Kmeans = round(DIKmeans,3),
     DBI.Wald = round(DBIWald,3), DBI.Kmeans = round(DBIKmeans,3),
     RI.Ward = round(RIWard,3), RI.Kmeans = round(RIKmeans,3),
     ARI.Ward = round(ARIWard,3), ARI.Kmeans = round(ARIKmeans,3))


table=data.frame(res)





if (print==TRUE) {
  cat("\n Dunn Indexes: ", round(DIWard,3), "(Ward)",round(DIKmeans,3), "(Kmeans)")

  cat("\n Davies-Bouldind Indexes: ", round(DBIWald,3), "(Ward)",round(DBIKmeans,3), "(Kmeans)")

  cat("\n Rand Indexes: ",round(RIWard,3),  "(Ward)", round(RIKmeans,3),"(Kmeans)")

  cat("\n Adjusted Rand Indexes: ",round(ARIWard,3), "(Ward)",round(ARIKmeans,3),"(Kmeans)")
}


if (impG==TRUE) {


col.Ward=as.numeric(CutTree)
col.Kmeans=as.numeric(KMeans$cluster)


res.PCA=SVDforPCA(x)

pdf(file.path(out.path,"PCAClust.pdf"), width=9, height=3)

par(mfrow=c(1,3))


  Xax=1
  Yax=2
  Xlim=c(min(res.PCA$pcs[,Xax])*1.2, max(res.PCA$pcs[,Xax])*1.2)
  Ylim=c(min(res.PCA$pcs[,Yax])*1.2, max(res.PCA$pcs[,Yax])*1.2)

  plot(res.PCA$pcs[,Xax], res.PCA$pcs[,Yax],xlim=Xlim, ylim=Ylim, pch=16, col=trcl1,
       xlab=paste0("PC",Xax," (", res.PCA$var[Xax] ,"%)"),
       ylab=paste0("PC",Yax," (", res.PCA$var[Yax] ,"%)"), main="PCA score plot \n True Groups")

  if (!is.null(RowNames)) {
  text(res.PCA$pcs[,Xax],res.PCA$pcs[,Yax],labels=RowNames,
       pos=c(2,3), col=trcl1)
  }
  abline(h=0, v=0)


  Xax=1
  Yax=2
  plot(res.PCA$pcs[,Xax],res.PCA$pcs[,Yax],xlim=Xlim, ylim=Ylim, pch=16,  col=col.Ward,
       xlab=paste0("PC",Xax," (", res.PCA$var[Xax] ,"%)"),
       ylab=paste0("PC",Yax," (", res.PCA$var[Yax] ,"%)"), main="PCA score plot \n  Method: Ward")
    if (!is.null(RowNames)) {
      text(res.PCA$pcs[,Xax],res.PCA$pcs[,Yax],labels=RowNames,
           pos=c(2,3), col=col.Ward)
    }

  abline(h=0, v=0)


  Xax=1
  Yax=2
  plot(res.PCA$pcs[,Xax],res.PCA$pcs[,Yax],xlim=Xlim, ylim=Ylim, pch=16, col=col.Kmeans,
       xlab=paste0("PC",Xax," (", res.PCA$var[Xax] ,"%)"),
       ylab=paste0("PC",Yax," (", res.PCA$var[Yax] ,"%)"), main="PCA score plot \n  Method: K-means")

  if (!is.null(RowNames)) {
  text(res.PCA$pcs[,Xax],res.PCA$pcs[,Yax],labels=RowNames,
       pos=c(2,3), col=col.Kmeans)
  }
  abline(h=0, v=0)


dev.off()

# Dendrogram Ward Classif
pdf(file.path(out.path,"DendroWardClust.pdf"))

A2Rplot(WardClassif,
        k=nClust,
        boxes = FALSE)

dev.off()


}


return(table)
}

