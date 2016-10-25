### CLUSTERING ##############################################

## Ward et K-means unsupervised clustering
## MIC criteria: Dunn, Davies-Bouldin, Rand and Adjusted-Rand


#' @export ClustMIC
#' @title Unsupervised clustering on (GPL) intensities and associated MIC indexes
#'
#' @description
#' Unsupervised clustering on (GPL) intensities based on Ward and K-Means algorithms.
#' Calculation of MIC statistical criteria of clustering quality:
#' Dunn, Davies-Bouldin, Rand and adjusted-Rand indexes.
#'
#' @param Intensities A numeric matrix of intensities (can be a pTreatGPL object).
#' @param nClust The number of groups to retrieve (donors, mixtures, ...).
#' @param Trcl The real groups' memberships of the samples, true class.
#' @param Dendr Logical argument (TRUE/FALSE) to obtain graphical dendrogram based on the Ward algorithm.
#'
#' @return A list of MIC quality indexes (Dunn, Davies-Bouldin, Rand and adjusted-Rand):
#' \describe{
#'   \item{\code{DunnW}}{Dunn index for Ward clustering}
#'   \item{\code{DunnKM}}{Dunn index for Kmeans clustering}
#'   \item{\code{DBW}}{Davies-Bouldin index for Ward clustering}
#'   \item{\code{DBKM}}{Davies-Bouldin index for Kmeans clustering}
#'   \item{\code{RandW}}{Rand index for Ward clustering}
#'   \item{\code{RandKM}}{Rand index for Kmeans clustering}
#'   \item{\code{AdjRandW}}{Adjusted Rand index for Ward clustering}
#'   \item{\code{AdjRandKM}}{Adjusted Rand index for Kmeans clustering}
#'    }
#' @author Baptiste Feraud, Manon Martin
#'
#' @examples
#' data("HumanSerum")
#' ClustMIC(Intensities = HumanSerumSpectra, nClust = 4, Trcl = ClassHS, Dendr = TRUE)
#'
#'@importFrom proxy dist
#'@importFrom clValid dunn
#'@importFrom clusterSim index.DB
#'@importFrom phyclust RRand
#'


ClustMIC = function(Intensities, nClust, Trcl, Dendr=TRUE){

# checks

if (missing(Intensities)){
warning("Intensities is missing with no default value")
}

if (missing(nClust)){
warning("nClust is missing with no default value")
}


if (sum(Trcl<=0) >0) {
  warning("Automatic rewritting of Trcl since values below 1 are not permitted")
  Trcl = Trcl+(1-min(Trcl))

}

if (! is.numeric(Intensities)) {
  stop(deparse(substitute(Intensities)), " is not numeric.")
}

  if (! is.numeric(nClust)) {
    stop(deparse(substitute(nClust)), " is not numeric.")
  }  else if (length(nClust)>1) {
    stop(deparse(substitute(nClust)), " has a length > 1.")
  }

  if (missing(Trcl) | sum(Trcl%%1!=0)>0){
    warning("Trcl is missing with no default value or is not an integer")
  }


  names(Trcl) = 1:length(Trcl)

if (! is.logical(Dendr)) {
    stop(deparse(substitute(Dendr)), " is not logical.")
}


#### WARD
DistI <- proxy::dist(Intensities, method ="euclidean", diag=FALSE, upper=FALSE)

WardI <- stats::hclust(DistI, "ward.D2")
clWI <- stats::cutree(WardI, nClust)


#### K-MEANS
clKMI <- stats::kmeans(Intensities, nClust, iter.max = 20, nstart = 1,
         algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"))


#### INDICES MIC ##############################

################################
##  Homogeneite des groupes   ##
################################

#### DUNN
DunnW <- clValid::dunn(DistI, clWI)
DunnKM <- clValid::dunn(DistI, clKMI$cluster)


#### DAVIES-BOULDIN
DBW <- clusterSim::index.DB(Intensities, clWI, DistI, centrotypes="medoids")
DBKM <- clusterSim::index.DB(Intensities, clKMI$cluster, DistI, centrotypes="medoids")


################################
##   Qualite du clustering    ##
################################

#### RAND
RandW <- phyclust::RRand(Trcl, clWI, lab = NULL)
RandKM <- phyclust::RRand(Trcl, clKMI$cluster, lab = NULL)


#### DENDROGRAM ################################

if (Dendr == TRUE) {

graphics::plot(WardI,lwd=1,cex=.5)
stats::rect.hclust(WardI, nClust, border="red")

}

res=list(DunnW = DunnW, DunnKM = DunnKM,
     DBW = DBW$DB, DBKM = DBKM$DB,
     RandW = RandW[[1]], RandKM = RandKM[[1]],
     AdjRandW = RandW[[2]], AdjRandKM = RandKM[[2]])


table=data.frame(res)
return(table)

}

