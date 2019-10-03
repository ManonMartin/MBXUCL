### CLUSTERING ON BINARY POSITIONS ######################################

## Ward unsupervised clustering MIC criteria: Dunn, Davies-Bouldin, Rand and
## Adjusted-Rand


#' @export binClustMIC
#' @title Unsupervised clustering on (GPL) binary presence/absence of peaks and associated MIC indexes
#'
#' @description
#' Unsupervised clustering on (GPL) binary presence/absence of peaks based on Ward algorithm.
#' Calculation of MIC statistical criteria of clustering quality:
#' Dunn, Davies-Bouldin, Rand and adjusted-Rand indexes.
#'
#' @param Positions A binary numeric matrix
#' @param Distance Choice of the distance measure adapted to binary objects ('Jaccard' or 'Ochiai')
#' @param nClust The number of groups to retrieve (donors, mixtures, ...).
#' @param Trcl The real groups' memberships of the samples, true class.
#' @param Dendr Logical argument (TRUE/FALSE) to obtain graphical dendrogram based on the Ward algorithm.
#'
#' @return A list of MIC quality indexes (Dunn, Davies-Bouldin, Rand and adjusted-Rand):
#' \describe{
#'   \item{\code{DunnW}}{Dunn index for Ward clustering}
#'   \item{\code{DBW}}{Davies-Bouldin index for Ward clustering}
#'   \item{\code{RandW}}{Rand index for Ward clustering}
#'   \item{\code{AdjRandW}}{Adjusted Rand index for Ward clustering}
#'    }
#' @author Baptiste Feraud
#'
#' @examples
#' Pos = t(GPL[, seq(4,dim(GPL)[2], 2)])
#' Trcl = c(rep(1,8), rep(2,8), rep(3,8))
#' binClustMIC(Positions = Pos, Distance = 'Jaccard', nClust = 3, Trcl, Dendr = TRUE)
#'
#'@importFrom proxy dist
#'@importFrom clValid dunn
#'@importFrom clusterSim index.DB
#'@importFrom phyclust RRand


binClustMIC <- function(Positions, Distance, nClust, Trcl, Dendr = TRUE) {

  # checks

  if (missing(Positions))  {
    warning("Positions is missing with no default value")
  }

  if (missing(Distance)) {
    warning("Distance is missing with no default value")
  }

  if (missing(nClust)) {
    warning("nClust is missing with no default value")
  }

  if (missing(Trcl))  {
    warning("Trcl is missing with no default value")
  }


  if (!is.numeric(Positions)) {
    stop(deparse(substitute(Positions)), " is not numeric.")
  }

  if (!(all(as.matrix(Positions) %in% 0:1)) == TRUE) {
    stop(deparse(substitute(Positions)), " is not a binary matrix.")
  }

  if (!(Distance %in% c("Jaccard", "Ochiai"))) {
    stop(deparse(substitute(Distance)), " is an incorrect distance measure for binary object.")
  }

  if (!is.numeric(nClust)) {
    stop(deparse(substitute(nClust)), " is not numeric.")
  } else if (length(nClust) > 1) {
    stop(deparse(substitute(nClust)), " has a length > 1.")
  }

  if (!is.numeric(Trcl))  {
    stop(deparse(substitute(Trcl)), " is not numeric.")
  }

  if (!is.logical(Dendr)) {
    stop(deparse(substitute(Dendr)), " is not logical.")
  }


  #### WARD
  DistP <- proxy::dist(Positions, method = Distance, diag = FALSE, upper = FALSE)

  WardP <- stats::hclust(DistP, "ward.D2")
  clWP <- stats::cutree(WardP, nClust)



  #### INDICES MIC ##############################

  ################################ Homogeneite des groupes ##

  #### DUNN
  DunnW <- clValid::dunn(DistP, clWP)


  #### DAVIES-BOULDIN
  DBW <- clusterSim::index.DB(Positions, clWP, DistP, centrotypes = "medoids")


  ################################ Qualite du clustering ##

  #### RAND
  RandW <- phyclust::RRand(Trcl, clWP, lab = NULL)


  #### DENDROGRAM ################################

  if (Dendr == TRUE) {

    graphics::plot(WardP, lwd = 1, cex = 0.5)
    stats::rect.hclust(WardP, nClust, border = "red")

  }

  res <- list(DunnW = DunnW, DBW = DBW$DB, RandW = RandW[[1]], AdjRandW = RandW[[2]])


  table <- data.frame(res)
  return(table)

}


