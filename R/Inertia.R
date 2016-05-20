#' @export Inertia
#' @title Measure of inertia
#'
#' @description
#' Calculates the within, between and total inertia of a data matrix
#'
#' @param x A numeric matrix.
#' @param y A numeric vector of data clusters.
#' @param print If \code{TRUE}, prints the inertia values
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{Name.group}}{Group membership}
#'   \item{\code{InertiaT}}{Total inertia}
#'   \item{\code{InertiaB}}{Between inertia}
#'   \item{\code{InertiaW}}{Within inertia}
#' }
#'
#' @examples
#' data("iris")
#' data = as.matrix((iris[,1:4]))
#' class = as.numeric(iris[,5])
#'
#' Res = Inertia(x = data, y = class, print = TRUE)
#'


Inertia = function(x, y, print=FALSE) {

# checks
if (! is.numeric(x)) {
  stop(deparse(substitute(x)), " is not numeric.")
}

if (! is.numeric(y)) {
  stop(deparse(substitute(y)), " is not numeric.")
}

if (! is.logical(print)) {
    stop(deparse(substitute(print)), " is not logical.")
}


# le calcul du barycentre de chaque groupe :
#############################################

nGroup=length(unique(y))

#Creation of group lists for data, y and barycenters
m=dim(x)[2]
n=dim(x)[1]

yG=c() # number of the class
nG=c() # number of elements for each group
xG=vector("list", nGroup) # data for elements of i's group
xmG=vector("list", nGroup) # barycenter for elements of i's group
j=1
for (i in unique(y)){
  yG[j]=i
  nG[j]=length(which(y==i))
  xG[[j]]=x[which(y==i),]
  xmG[[j]]=apply(xG[[j]],2,mean)
  j=j+1
}


# Raw mean
G=apply(x,2,mean)


#### Calculation of Between Inertia  (sum ng*d2(ug,u))
IB=c()
for (i in 1:length(unique(y))){
  ib = nG[i]*(sum((xmG[[i]]-G)^2))
  IB=c(IB,ib)
  InertiaB=sum(IB)
}


#### Calculation of Total Inertia
xmG=matrix(G,nrow=n,ncol=m,byrow=TRUE)
it=(x-xmG)^2
IT=apply(it,1,sum)


InertiaT = sum(IT)


iner.inter <- (InertiaB/InertiaT) # Between-group Inertia (% of tot inertia) :
iner.inter100=iner.inter*100


#### Calculation of Within Inertia
iner.intra = 1-iner.inter
InertiaW=(iner.intra * InertiaT)



iner.intra100=iner.intra*100 # Within-group Inertia (% of tot inertia) :



if (print==TRUE) {
  cat("\n Total Inertia",InertiaT)
  cat("\n Intergroup Inertia",InertiaB)
  cat("\n Intergroup Inertia (% of Itot)",iner.inter100)
  cat("\n Within Inertia ",InertiaW)
  cat("\n Within Inertia (% of Itot)",iner.intra, "\n")
}



# FIXME

#Within Inertia for each group
# dG=vector("list", nGroup)
# sumdG=vector("list", nGroup)
# IG=c()
# IGprop=c()
#
# for (i in 1:length(unique(y))){
#
# xmG[[i]]=matrix(data=as.vector(xmG[[i]]),nrow=nG[i],ncol=m,byrow=TRUE)
#
# dG[[i]]=(xG[[i]]-xmG[[i]])^2
#
# #print(dim(dG1))
# sumdG[[i]]=apply(dG[[i]],1,sum)
# IG[i]=(sum(sumdG[[i]]))
# IGprop[i]=(IG[i]/InertiaT)*100
# }
#
# cat("\n Within Inertia in each group",IG)
# cat("\n InWithin Inertia in each group (%)",IGprop)
#
# TWI=sum(IG)
# cat("\n Total Within Inertia",TWI)


PInertiaT = InertiaT*100/InertiaT
PInertiaW = InertiaW*100/InertiaT
PInertiaB = InertiaB*100/InertiaT


res=matrix(data = round(c(InertiaB, InertiaW, InertiaT, PInertiaB, PInertiaW, PInertiaT),2), ncol = 3,
           dimnames = list(c("Value", "Percentage"), c("BI", "WI", "TI")), byrow = TRUE)

return(res)


}
