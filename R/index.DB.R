#' Function from ClusterSim (version 0.49-2)
#' @param x data
#' @param cl vector of integers indicating the cluster to which each object is allocated
#' @param d optional distance matrix, used for calculations if centrotypes="medoids"
#' @param centrotypes "centroids" or "medoids"
#' @param p the power of the Minkowski distance between centroids or medoids of clusters: p=1 - Manhattan distance; p=2 - Euclidean distance
#' @param q the power of dispersion measure of a cluster: q=1 - the average distance of objects in the r-th cluster to the centroid or medoid of the r-th cluster; q=2 - the standard deviation of the distance of objects in the r-th cluster to the centroid or medoid of the r-th cluster
#' @NoRd


index.DB <- function(x, cl, d = NULL, centrotypes = "centroids", 
                      p = 2, q = 2) {
  
  ###### .medoid function
  .medoid<-function(x,d)
  {
    minj<-0
    minsumdist<-sum(d)
    if(is.null(dim(x)) && is.null(dim(d))){
      dim(x)<-c(1,length(x))
      x
    }
    else{
      if(is.null(dim(d))){
        dim(d)<-c(1,1)
      }
      if(is.null(dim(x))){
        dim(x)<-c(length(x),1)
      }
      for(j in 1:nrow(d)){
        if (sum(d[j,])<=minsumdist){
          #minj<-row.names(d)[j]
          minj<-j
          minsumdist<-sum(d[j,])
        }
      }
      resul<-as.matrix(x[minj,])  
      resul
    }
  }
  
  
  
  ######
  
  if (sum(c("centroids", "medoids") == centrotypes) == 0) 
    stop("Wrong centrotypes argument")
  if ("medoids" == centrotypes && is.null(d)) 
    stop("For argument centrotypes = 'medoids' d cannot be null")
  if (!is.null(d)) {
    if (!is.matrix(d)) {
      d <- as.matrix(d)
    }
    row.names(d) <- row.names(x)
  }
  if (is.null(dim(x))) {
    dim(x) <- c(length(x), 1)
  }
  x <- as.matrix(x)
  n <- length(cl)
  k <- max(cl)
  dAm <- d
  centers <- matrix(nrow = k, ncol = ncol(x))
  if (centrotypes == "centroids") {
    for (i in 1:k) {
      for (j in 1:ncol(x)) {
        centers[i, j] <- mean(x[cl == i, j])
      }
    }
  }
  else if (centrotypes == "medoids") {
    for (i in 1:k) {
      clAi <- dAm[cl == i, cl == i]
      if (is.null(clAi)) {
        centers[i, ] <- NULL
      }
      else {
        centers[i, ] <- .medoid(x[cl == i, ], dAm[cl == 
                                                    i, cl == i])
      }
    }
  }
  else {
    stop("wrong centrotypes argument")
  }
  S <- rep(0, k)
  for (i in 1:k) {
    ind <- (cl == i)
    if (sum(ind) > 1) {
      centerI <- centers[i, ]
      centerI <- rep(centerI, sum(ind))
      centerI <- matrix(centerI, nrow = sum(ind), ncol = ncol(x), 
                        byrow = TRUE)
      S[i] <- mean(sqrt(apply((x[ind, ] - centerI)^2, 1, 
                              sum))^q)^(1/q)
    }
    else S[i] <- 0
  }
  M <- as.matrix(dist(centers, method = "minkowski", p = p))
  R <- array(Inf, c(k, k))
  r = rep(0, k)
  for (i in 1:k) {
    for (j in 1:k) {
      R[i, j] = (S[i] + S[j])/M[i, j]
    }
    r[i] = max(R[i, ][is.finite(R[i, ])])
  }
  DB = mean(r[is.finite(r)])
  resul <- list(DB = DB, r = r, R = R, d = M, S = S, centers = centers)
  resul
}